//! `ngs-dose count`: one pass (scan) or targeted retrieval (fetch) over a BAM/CRAM, producing
//! the per-sample sufficient statistics that every downstream estimate is computed from.
//!
//! Counting rules, identical for controls and targets so that they cancel in the ratio:
//!   * primary alignments only (secondary, supplementary and QC-fail records are skipped);
//!   * the duplicate flag is IGNORED. MarkDuplicates keys on mapped position, and in a collapsed
//!     multi-copy locus mates scatter over paralogs, so duplicates are under-flagged there
//!     relative to single-copy sequence. Excluding flagged reads biases the ratio by a
//!     sample-specific amount; counting everything does not.
//!   * no MAPQ filter. Reads of a multi-copy class are MAPQ 0 by construction.
//!   * the unit of counting is the fragment 5' end (strand-aware, soft clips included), which
//!     is independent of read length.

use crate::controls::{Controls, GC_BINS};
use crate::kmer::{unpack_bam_seq, KmerIter};
use crate::panel::{Entry, Kind, Panel};
use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::path::PathBuf;

pub const INSERT_MAX: usize = 1500;
pub const READLEN_MAX: usize = 400;

#[derive(Clone)]
pub struct Params {
    pub min_hits: usize,
    pub min_frac: f64,
    pub bin: usize,
    pub l_grid: Vec<usize>,
    pub threads: usize,
    pub place_bin: i64,
    pub unmapped: bool,
    pub pad: i64,
    /// fetch mode: attempts per interval before giving up (transient network errors)
    pub retries: usize,
}

#[derive(Clone)]
pub struct ClassAcc {
    pub reads: u64,
    pub bases: u64,
    pub dup_flagged: u64,
    pub out_of_range: u64,
    pub fwd: Vec<u32>,
    pub rev: Vec<u32>,
    pub gc_read: Vec<u64>,
    pub hit_frac: Vec<u64>,
}

#[derive(Clone)]
pub struct Acc {
    pub classes: Vec<ClassAcc>,
    pub placements: FxHashMap<(u8, i32, i32), u32>,
    pub ambiguous: u64,
    pub below_threshold: u64,
}

/// Ticks whenever any thread reads records or finishes an interval: what the watchdog watches.
pub static PROGRESS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

#[inline]
fn tick() {
    PROGRESS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

/// A dead HTTPS connection (or a hung network filesystem) does not fail, it waits - for ever.
/// Errors are retried per interval; a stall cannot be, because the thread that would retry is
/// the one that is blocked. So a watchdog ends the process with EX_TEMPFAIL (75) when nothing
/// has been read for `limit_s` seconds, and whoever launched it tries again. Nothing is lost:
/// output is only written, under a temporary name, once counting is complete.
pub fn spawn_watchdog(limit_s: u64, what: String) {
    use std::sync::atomic::Ordering::Relaxed;
    std::thread::spawn(move || {
        let mut last = PROGRESS.load(Relaxed);
        let mut since = std::time::Instant::now();
        loop {
            std::thread::sleep(std::time::Duration::from_secs(limit_s.clamp(1, 5)));
            let now = PROGRESS.load(Relaxed);
            if now != last {
                last = now;
                since = std::time::Instant::now();
            } else if since.elapsed().as_secs() >= limit_s {
                eprintln!(
                    "[count] {}: nothing read for {} s (stalled connection or filesystem); exiting with status 75 so that the caller can retry",
                    what, limit_s
                );
                std::process::exit(75);
            }
        }
    });
}

#[derive(Clone, Default)]
pub struct ReadStats {
    pub records: u64,
    pub primary: u64,
    pub primary_dup: u64,
    pub unmapped: u64,
    pub ctrl: u64,
    pub ctrl_dup: u64,
    pub ctrl_mapq0: u64,
    pub ctrl_mapq_lt20: u64,
    pub insert: Vec<u64>,
    pub readlen: Vec<u64>,
    /// mapped primary reads per header contig; complete in scan mode only
    pub contig_reads: Vec<u64>,
}

impl ReadStats {
    pub fn new() -> Self {
        ReadStats { insert: vec![0; INSERT_MAX + 1], readlen: vec![0; READLEN_MAX + 1], ..Default::default() }
    }
    pub fn merge(&mut self, o: &ReadStats) {
        self.records += o.records;
        self.primary += o.primary;
        self.primary_dup += o.primary_dup;
        self.unmapped += o.unmapped;
        self.ctrl += o.ctrl;
        self.ctrl_dup += o.ctrl_dup;
        self.ctrl_mapq0 += o.ctrl_mapq0;
        self.ctrl_mapq_lt20 += o.ctrl_mapq_lt20;
        for (a, b) in self.insert.iter_mut().zip(&o.insert) {
            *a += b;
        }
        for (a, b) in self.readlen.iter_mut().zip(&o.readlen) {
            *a += b;
        }
        if self.contig_reads.len() < o.contig_reads.len() {
            self.contig_reads.resize(o.contig_reads.len(), 0);
        }
        for (a, b) in self.contig_reads.iter_mut().zip(&o.contig_reads) {
            *a += b;
        }
    }
}

impl Acc {
    pub fn new(panel: &Panel, bin: usize) -> Acc {
        let classes = panel
            .classes
            .iter()
            .map(|c| {
                let nb = if c.kind == Kind::Positional { c.length.div_ceil(bin) } else { 0 };
                ClassAcc {
                    reads: 0,
                    bases: 0,
                    dup_flagged: 0,
                    out_of_range: 0,
                    fwd: vec![0; nb],
                    rev: vec![0; nb],
                    gc_read: vec![0; GC_BINS],
                    hit_frac: vec![0; 11],
                }
            })
            .collect();
        Acc { classes, placements: FxHashMap::default(), ambiguous: 0, below_threshold: 0 }
    }
    pub fn merge(&mut self, o: &Acc) {
        for (a, b) in self.classes.iter_mut().zip(&o.classes) {
            a.reads += b.reads;
            a.bases += b.bases;
            a.dup_flagged += b.dup_flagged;
            a.out_of_range += b.out_of_range;
            for (x, y) in a.fwd.iter_mut().zip(&b.fwd) {
                *x += y;
            }
            for (x, y) in a.rev.iter_mut().zip(&b.rev) {
                *x += y;
            }
            for (x, y) in a.gc_read.iter_mut().zip(&b.gc_read) {
                *x += y;
            }
            for (x, y) in a.hit_frac.iter_mut().zip(&b.hit_frac) {
                *x += y;
            }
        }
        for (k, v) in &o.placements {
            *self.placements.entry(*k).or_insert(0) += v;
        }
        self.ambiguous += o.ambiguous;
        self.below_threshold += o.below_threshold;
    }
}

/// Result of classifying one *stored* sequence against the panel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Assignment {
    pub class: usize,
    /// unit coordinate of the first stored base (may fall outside [0, len) before wrapping)
    pub pos5: i64,
    /// true if the stored sequence is the reverse complement of the unit
    pub reverse: bool,
    pub hits: usize,
    pub valid: usize,
}

impl Assignment {
    /// Fragment 5' end and strand on the unit for the read *as sequenced*.
    ///
    /// BAM stores a reverse-strand alignment reverse-complemented, so for a record with the
    /// reverse flag the sequenced read is the reverse complement of the stored sequence and
    /// its 5' end is the stored sequence's last base.
    #[inline]
    pub fn sequenced(&self, read_len: usize, record_reverse: bool) -> (i64, bool) {
        let span = read_len as i64 - 1;
        match (self.reverse, record_reverse) {
            (false, false) => (self.pos5, false),
            (false, true) => (self.pos5 + span, true),
            (true, false) => (self.pos5, true),
            (true, true) => (self.pos5 - span, false),
        }
    }
}

pub struct Scratch {
    codes: Vec<u8>,
    hits: Vec<(u32, Entry, bool)>,
    diags: Vec<i64>,
    class_hits: Vec<u32>,
    out: Vec<Assignment>,
}

impl Scratch {
    pub fn new(n_classes: usize) -> Scratch {
        Scratch { codes: Vec::new(), hits: Vec::new(), diags: Vec::new(), class_hits: vec![0; n_classes], out: Vec::new() }
    }
}

/// What happened to a read that carried panel k-mers but produced no assignment (diagnostics).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Unassigned {
    pub below: bool,
    pub ambiguous: bool,
}

/// Classify a read given its 2-bit codes (already in `scratch.codes`); assignments are appended
/// to `out`.
///
/// Every k-mer is tested: an earlier version probed every 4th k-mer first, which dropped
/// 0.04-0.06% of class reads (hits in short, broken runs) while being described as lossless.
///
/// A class is judged on its own hits, so its count does not depend on which other classes are
/// loaded: a distal-junction read that also carries beta-satellite k-mers is a DJ fragment end
/// whether or not a satellite panel is present. Positional classes are therefore independent
/// (a read can mark a fragment end in one and add mass to a family). Only compositional
/// families compete with one another, winner takes all, because a read's bases can belong to
/// one family only; a read split more evenly than 5:1 between two families is left unassigned.
pub fn classify(panel: &Panel, p: &Params, s: &mut Scratch, out: &mut Vec<Assignment>) -> Unassigned {
    let k = panel.k;
    out.clear();
    s.hits.clear();
    let mut valid = 0usize;
    for km in KmerIter::new(&s.codes, k) {
        valid += 1;
        if let Some(e) = panel.get(km.canon) {
            s.hits.push((km.offset as u32, e, km.fwd_is_canon));
        }
    }
    let mut un = Unassigned::default();
    if s.hits.is_empty() {
        return un;
    }
    for h in s.class_hits.iter_mut() {
        *h = 0;
    }
    for &(_, e, _) in &s.hits {
        s.class_hits[e.class()] += 1;
    }
    let passes = |h: u32| h as usize >= p.min_hits && (h as f64) >= p.min_frac * valid as f64;
    let (mut comp_best, mut comp_best_hits, mut comp_second) = (usize::MAX, 0u32, 0u32);
    for (class, &h) in s.class_hits.iter().enumerate() {
        if h == 0 {
            continue;
        }
        match panel.classes[class].kind {
            Kind::Compositional => {
                if h > comp_best_hits {
                    comp_second = comp_best_hits;
                    comp_best_hits = h;
                    comp_best = class;
                } else if h > comp_second {
                    comp_second = h;
                }
            }
            Kind::Positional => {
                if !passes(h) {
                    un.below = true;
                    continue;
                }
                // orientation by majority, position by the median diagonal of the agreeing hits
                let same = s.hits.iter().filter(|&&(_, e, rf)| e.class() == class && rf == e.fwd()).count();
                let reverse = h as usize - same > same;
                s.diags.clear();
                for &(off, e, rf) in &s.hits {
                    if e.class() == class && ((rf == e.fwd()) != reverse) {
                        let d = if reverse { e.pos() as i64 + k as i64 - 1 + off as i64 } else { e.pos() as i64 - off as i64 };
                        s.diags.push(d);
                    }
                }
                let mid = s.diags.len() / 2;
                let (_, m, _) = s.diags.select_nth_unstable(mid);
                out.push(Assignment { class, pos5: *m, reverse, hits: h as usize, valid });
            }
        }
    }
    if comp_best != usize::MAX {
        if !passes(comp_best_hits) {
            un.below = true;
        } else if comp_second * 5 > comp_best_hits {
            un.ambiguous = true;
        } else {
            out.push(Assignment { class: comp_best, pos5: 0, reverse: false, hits: comp_best_hits as usize, valid });
        }
    }
    un
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn record_assignment(
    panel: &Panel,
    p: &Params,
    acc: &mut Acc,
    a: &Assignment,
    read_len: usize,
    gc: u32,
    dup: bool,
    rec_reverse: bool,
    tid: i32,
    pos: i64,
) {
    let cdef = &panel.classes[a.class];
    let ca = &mut acc.classes[a.class];
    ca.reads += 1;
    ca.bases += read_len as u64;
    ca.dup_flagged += dup as u64;
    if read_len > 0 {
        ca.gc_read[((gc as usize * 100) + read_len / 2) / read_len] += 1;
    }
    ca.hit_frac[(a.hits * 10 / a.valid.max(1)).min(10)] += 1;
    if cdef.kind == Kind::Positional {
        let len = cdef.length as i64;
        let (mut x, reverse) = a.sequenced(read_len, rec_reverse);
        if cdef.circular {
            x = x.rem_euclid(len);
        }
        if x < 0 || x >= len {
            ca.out_of_range += 1;
        } else {
            let b = x as usize / p.bin;
            if reverse {
                ca.rev[b] += 1
            } else {
                ca.fwd[b] += 1
            }
        }
    }
    let pb = if tid < 0 { -1 } else { (pos / p.place_bin) as i32 };
    *acc.placements.entry((a.class as u8, tid, pb)).or_insert(0) += 1;
}

/// Fragment 5' end on the reference (inclusive coordinate), soft clips included.
#[inline]
fn five_prime(rec: &bam::Record) -> i64 {
    let cig = rec.raw_cigar();
    if rec.is_reverse() {
        let mut end = rec.pos();
        for &c in cig {
            match c & 0xf {
                0 | 2 | 3 | 7 | 8 => end += (c >> 4) as i64,
                _ => {}
            }
        }
        let clip = cig.last().filter(|&&c| c & 0xf == 4).map(|&c| (c >> 4) as i64).unwrap_or(0);
        end + clip - 1
    } else {
        let clip = cig.first().filter(|&&c| c & 0xf == 4).map(|&c| (c >> 4) as i64).unwrap_or(0);
        rec.pos() - clip
    }
}

const SKIP_FLAGS: u16 = 0x100 | 0x200 | 0x800;

/// Control bookkeeping for one record. Returns false if the record is to be skipped entirely.
/// With `pending`, control increments are deferred so that a failed interval can be retried
/// without double counting; without it they are applied at once.
#[inline]
fn handle_record_meta(rec: &bam::Record, controls: &Controls, st: &mut ReadStats, pending: Option<&mut Vec<(usize, usize, bool)>>) -> bool {
    st.records += 1;
    if st.records & 0x3ff == 0 {
        tick();
    }
    let flags = rec.flags();
    if flags & SKIP_FLAGS != 0 {
        return false;
    }
    st.primary += 1;
    let dup = flags & 0x400 != 0;
    st.primary_dup += dup as u64;
    if flags & 0x4 != 0 {
        st.unmapped += 1;
        return true;
    }
    let tid = rec.tid();
    if tid >= 0 {
        if st.contig_reads.len() <= tid as usize {
            st.contig_reads.resize(tid as usize + 1, 0);
        }
        st.contig_reads[tid as usize] += 1;
    }
    let p5 = five_prime(rec);
    let Some((idx, off, is_control)) = controls.locate(rec.tid(), p5) else {
        return true;
    };
    match pending {
        Some(v) => v.push((idx, off, rec.is_reverse())),
        None => controls.bump(idx, off, rec.is_reverse()),
    }
    if is_control {
        st.ctrl += 1;
        st.ctrl_dup += dup as u64;
        let q = rec.mapq();
        st.ctrl_mapq0 += (q == 0) as u64;
        st.ctrl_mapq_lt20 += (q < 20) as u64;
        let t = rec.insert_size();
        if flags & 0x2 != 0 && t > 0 {
            st.insert[(t as usize).min(INSERT_MAX)] += 1;
        }
        st.readlen[rec.seq_len().min(READLEN_MAX)] += 1;
    }
    true
}

fn set_cram_fields<R: Read>(rd: &mut R) {
    let fields = htslib::sam_fields_SAM_FLAG
        | htslib::sam_fields_SAM_RNAME
        | htslib::sam_fields_SAM_POS
        | htslib::sam_fields_SAM_MAPQ
        | htslib::sam_fields_SAM_CIGAR
        | htslib::sam_fields_SAM_RNEXT
        | htslib::sam_fields_SAM_PNEXT
        | htslib::sam_fields_SAM_TLEN
        | htslib::sam_fields_SAM_SEQ;
    // harmless no-op for BAM
    let _ = rd.set_cram_options(htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields);
}

pub struct Input {
    pub path: String,
    pub index: Option<String>,
    pub reference: Option<PathBuf>,
}

impl Input {
    fn is_url(&self) -> bool {
        self.path.contains("://")
    }
    fn open_indexed(&self) -> Result<bam::IndexedReader> {
        let mut rd = if self.is_url() {
            let full = match &self.index {
                Some(i) => format!("{}##idx##{}", self.path, i),
                None => self.path.clone(),
            };
            let url = url::Url::parse(&full).with_context(|| format!("bad URL {}", full))?;
            bam::IndexedReader::from_url(&url)
        } else {
            match &self.index {
                Some(i) => bam::IndexedReader::from_path_and_index(&self.path, i),
                None => bam::IndexedReader::from_path(&self.path),
            }
        }
        .map_err(|e| anyhow::anyhow!("cannot open {} with its index: {}", self.path, e))?;
        if let Some(r) = &self.reference {
            rd.set_reference(r).map_err(|e| anyhow::anyhow!("set_reference: {}", e))?;
        }
        set_cram_fields(&mut rd);
        Ok(rd)
    }
    fn open_stream(&self, threads: usize) -> Result<bam::Reader> {
        let mut rd = if self.is_url() { bam::Reader::from_url(&url::Url::parse(&self.path)?) } else { bam::Reader::from_path(&self.path) }
            .map_err(|e| anyhow::anyhow!("cannot open {}: {}", self.path, e))?;
        if let Some(r) = &self.reference {
            rd.set_reference(r).map_err(|e| anyhow::anyhow!("set_reference: {}", e))?;
        }
        set_cram_fields(&mut rd);
        if threads > 1 {
            rd.set_threads(threads).map_err(|e| anyhow::anyhow!("set_threads: {}", e))?;
        }
        Ok(rd)
    }
}

struct Batch {
    packed: Vec<u8>,
    items: Vec<(u32, u32, i32, i64, bool, bool)>, // offset, len, tid, pos, dup, record is reverse-strand
}

impl Batch {
    fn new() -> Batch {
        Batch { packed: Vec::with_capacity(1 << 21), items: Vec::with_capacity(1 << 14) }
    }
}

fn classify_batch(panel: &Panel, p: &Params, b: &Batch, acc: &mut Acc, s: &mut Scratch) {
    for &(off, len, tid, pos, dup, rev) in &b.items {
        let (off, len) = (off as usize, len as usize);
        let gc = unpack_bam_seq(&b.packed[off..off + len.div_ceil(2)], len, &mut s.codes);
        let mut out = std::mem::take(&mut s.out);
        let un = classify(panel, p, s, &mut out);
        for a in &out {
            record_assignment(panel, p, acc, a, len, gc, dup, rev, tid, pos);
        }
        s.out = out;
        acc.ambiguous += un.ambiguous as u64;
        acc.below_threshold += un.below as u64;
    }
}

/// Whole-file pass: every primary read is screened against the panel wherever it was aligned.
pub fn scan(input: &Input, panel: &Panel, controls: &Controls, p: &Params) -> Result<(Acc, ReadStats)> {
    let workers = p.threads.max(1);
    let decomp = p.threads.clamp(1, 8);
    let mut rd = input.open_stream(decomp)?;
    let (tx, rx) = crossbeam_channel::bounded::<Batch>(workers * 2);
    let mut stats = ReadStats::new();
    let accs: Vec<Acc> = std::thread::scope(|scope| -> Result<Vec<Acc>> {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                let rx = rx.clone();
                scope.spawn(move || {
                    let mut acc = Acc::new(panel, p.bin);
                    let mut s = Scratch::new(panel.classes.len());
                    for b in rx.iter() {
                        classify_batch(panel, p, &b, &mut acc, &mut s);
                    }
                    acc
                })
            })
            .collect();
        drop(rx);
        let mut rec = bam::Record::new();
        let mut batch = Batch::new();
        while let Some(r) = rd.read(&mut rec) {
            r.map_err(|e| anyhow::anyhow!("read error: {}", e))?;
            if !handle_record_meta(&rec, controls, &mut stats, None) {
                continue;
            }
            let len = rec.seq_len();
            if len < panel.k {
                continue;
            }
            let enc = rec.seq().encoded;
            let rev = rec.is_reverse() && !rec.is_unmapped();
            batch.items.push((batch.packed.len() as u32, len as u32, rec.tid(), rec.pos(), rec.is_duplicate(), rev));
            batch.packed.extend_from_slice(&enc[..len.div_ceil(2)]);
            if batch.items.len() >= 16384 {
                tx.send(std::mem::replace(&mut batch, Batch::new())).ok();
            }
        }
        tx.send(batch).ok();
        drop(tx);
        Ok(handles.into_iter().map(|h| h.join().expect("worker panicked")).collect())
    })?;
    let mut acc = Acc::new(panel, p.bin);
    for a in &accs {
        acc.merge(a);
    }
    Ok((acc, stats))
}

/// Targeted pass: only reads placed in the control regions and in the class sinks are
/// retrieved (plus, optionally, the unmapped bin).
pub fn fetch(input: &Input, panel: &Panel, controls: &Controls, sinks: &[(String, i64, i64)], p: &Params) -> Result<(Acc, ReadStats)> {
    // merged, non-overlapping plan; each record is counted in the interval holding its start
    let mut by_chrom: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
    for (c, s, e) in controls.fetch_intervals(p.pad).into_iter().chain(sinks.iter().cloned()) {
        by_chrom.entry(c).or_default().push((s, e));
    }
    let probe = input.open_indexed()?;
    let header = probe.header().clone();
    drop(probe);
    let mut plan: Vec<(i32, i64, i64)> = Vec::new();
    for (c, mut v) in by_chrom {
        let tid = match header.tid(c.as_bytes()) {
            Some(t) => t as i32,
            None => bail!("fetch interval contig '{}' is not in the alignment header", c),
        };
        v.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::new();
        for (s, e) in v {
            match merged.last_mut() {
                Some(l) if s <= l.1 => l.1 = l.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        plan.extend(merged.into_iter().map(|(s, e)| (tid, s, e)));
    }
    plan.sort_unstable();
    let workers = p.threads.max(1).min(plan.len().max(1));
    let (tx, rx) = crossbeam_channel::unbounded::<Option<(i32, i64, i64)>>();
    for iv in &plan {
        tx.send(Some(*iv)).ok();
    }
    if p.unmapped {
        tx.send(None).ok();
    }
    drop(tx);
    let results: Vec<Result<(Acc, ReadStats)>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                let rx = rx.clone();
                scope.spawn(move || -> Result<(Acc, ReadStats)> {
                    let mut rd = input.open_indexed()?;
                    let mut acc = Acc::new(panel, p.bin);
                    let mut st = ReadStats::new();
                    let mut s = Scratch::new(panel.classes.len());
                    let mut rec = bam::Record::new();
                    let mut pending: Vec<(usize, usize, bool)> = Vec::new();
                    for job in rx.iter() {
                        let mut attempt = 0usize;
                        loop {
                            // an interval is committed only if it was read to its end
                            let mut iv_acc = Acc::new(panel, p.bin);
                            let mut iv_st = ReadStats::new();
                            pending.clear();
                            let res = read_interval(&mut rd, job, panel, controls, p, &mut iv_acc, &mut iv_st, &mut pending, &mut s, &mut rec);
                            match res {
                                Ok(()) => {
                                    tick();
                                    for &(idx, off, rev) in &pending {
                                        controls.bump(idx, off, rev);
                                    }
                                    acc.merge(&iv_acc);
                                    st.merge(&iv_st);
                                    break;
                                }
                                Err(e) => {
                                    tick();
                                    attempt += 1;
                                    if attempt >= p.retries.max(1) {
                                        return Err(e.context(format!("interval {:?} failed {} times", job, attempt)));
                                    }
                                    eprintln!("[count] retry {}/{} for interval {:?}: {}", attempt, p.retries, job, e);
                                    std::thread::sleep(std::time::Duration::from_millis(500 << attempt.min(5)));
                                    rd = input.open_indexed()?;
                                }
                            }
                        }
                    }
                    Ok((acc, st))
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("worker panicked")).collect()
    });
    let mut acc = Acc::new(panel, p.bin);
    let mut stats = ReadStats::new();
    for r in results {
        let (a, s) = r?;
        acc.merge(&a);
        stats.merge(&s);
    }
    Ok((acc, stats))
}

#[allow(clippy::too_many_arguments)]
fn read_interval(
    rd: &mut bam::IndexedReader,
    job: Option<(i32, i64, i64)>,
    panel: &Panel,
    controls: &Controls,
    p: &Params,
    acc: &mut Acc,
    st: &mut ReadStats,
    pending: &mut Vec<(usize, usize, bool)>,
    s: &mut Scratch,
    rec: &mut bam::Record,
) -> Result<()> {
    match job {
        Some((tid, a, b)) => rd.fetch((tid, a, b)),
        None => rd.fetch(bam::FetchDefinition::Unmapped),
    }
    .map_err(|e| anyhow::anyhow!("fetch failed: {}", e))?;
    while let Some(r) = rd.read(rec) {
        r.map_err(|e| anyhow::anyhow!("read error: {}", e))?;
        if let Some((_, a, b)) = job {
            if rec.pos() < a || rec.pos() >= b {
                continue; // belongs to a neighbouring interval, or starts before this one
            }
        }
        if !handle_record_meta(rec, controls, st, Some(pending)) {
            continue;
        }
        let len = rec.seq_len();
        if len < panel.k {
            continue;
        }
        let gc = unpack_bam_seq(rec.seq().encoded, len, &mut s.codes);
        let mut out = std::mem::take(&mut s.out);
        let un = classify(panel, p, s, &mut out);
        let rev = rec.is_reverse() && !rec.is_unmapped();
        for asg in &out {
            record_assignment(panel, p, acc, asg, len, gc, rec.is_duplicate(), rev, rec.tid(), rec.pos());
        }
        s.out = out;
        acc.ambiguous += un.ambiguous as u64;
        acc.below_threshold += un.below as u64;
    }
    Ok(())
}

// ---------------------------------------------------------------- output

#[derive(Serialize)]
pub struct GcTable {
    pub l: usize,
    pub n: Vec<u64>,
    pub o: Vec<u64>,
}

#[derive(Serialize)]
pub struct RegionOut {
    pub name: String,
    pub role: String,
    pub label: String,
    pub len: usize,
    pub gc: f64,
    pub obs: u64,
    /// the region's contig is not in the alignment header
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub absent: bool,
}

#[derive(Serialize)]
pub struct ClassOut {
    pub name: String,
    pub kind: String,
    pub length: usize,
    pub circular: bool,
    pub panel_kmers: usize,
    pub reads: u64,
    pub bases: u64,
    pub dup_flagged: u64,
    pub out_of_range: u64,
    pub bin: usize,
    pub fwd: Vec<u32>,
    pub rev: Vec<u32>,
    pub gc_read: Vec<u64>,
    pub hit_frac: Vec<u64>,
}

#[derive(Serialize)]
pub struct Placement {
    pub class: String,
    pub contig: String,
    pub start: i64,
    pub reads: u32,
}

/// A contig of the alignment header. Lengths let the estimator refuse a file aligned to another
/// reference build; in scan mode `reads` (mapped primary reads) gives the dosage of whatever
/// else the reference carries: chrM, chrEBV, chrY, decoys.
#[derive(Serialize)]
pub struct ContigOut {
    pub name: String,
    pub len: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reads: Option<u64>,
}

#[derive(Serialize)]
pub struct Output {
    pub format: &'static str,
    pub engine_version: &'static str,
    pub sample: String,
    pub input: String,
    pub mode: String,
    pub k: usize,
    pub min_hits: usize,
    pub min_frac: f64,
    pub panel: String,
    pub panel_sha256: Vec<String>,
    pub controls_sha256: String,
    pub sinks_sha256: Option<String>,
    pub controls: String,
    pub sinks: Option<String>,
    pub unmapped_fetched: bool,
    pub records: u64,
    pub primary: u64,
    pub primary_dup_flagged: u64,
    pub unmapped: u64,
    pub ctrl_reads: u64,
    pub ctrl_dup_flagged: u64,
    pub ctrl_mapq0: u64,
    pub ctrl_mapq_lt20: u64,
    pub ctrl_positions: u64,
    pub read_length_mode: usize,
    pub insert_median: usize,
    pub insert_hist_10bp: Vec<u64>,
    pub readlen_hist: Vec<(usize, u64)>,
    pub gc_tables: Vec<GcTable>,
    pub regions: Vec<RegionOut>,
    pub classes: Vec<ClassOut>,
    pub ambiguous_reads: u64,
    pub below_threshold_reads: u64,
    pub placement_bin: i64,
    pub placements: Vec<Placement>,
    pub eof_marker: &'static str,
    pub contigs: Vec<ContigOut>,
    pub elapsed_sec: f64,
}

fn median_from_hist(h: &[u64]) -> usize {
    let total: u64 = h.iter().sum();
    if total == 0 {
        return 0;
    }
    let mut c = 0u64;
    for (i, &v) in h.iter().enumerate() {
        c += v;
        if c * 2 >= total {
            return i;
        }
    }
    h.len() - 1
}

#[allow(clippy::too_many_arguments)]
pub fn make_output(
    sample: String,
    input: &Input,
    mode: &str,
    panel: &Panel,
    panel_path: String,
    hashes: (Vec<String>, String, Option<String>),
    controls: &Controls,
    controls_path: String,
    sinks_path: Option<String>,
    p: &Params,
    acc: Acc,
    st: ReadStats,
    header_contigs: &[(String, u64)],
    sink_contigs: &[String],
    eof_marker: &'static str,
    elapsed: f64,
) -> Output {
    let read_mode = st.readlen.iter().enumerate().max_by_key(|(_, &v)| v).map(|(i, _)| i).unwrap_or(0);
    let insert_median = median_from_hist(&st.insert);
    let mut grid = p.l_grid.clone();
    if read_mode > 0 && !grid.contains(&read_mode) {
        grid.push(read_mode);
    }
    grid.sort_unstable();
    grid.dedup();
    let max_l = controls.regions.iter().map(|r| r.flank).min().unwrap_or(0);
    let gc_tables = grid
        .into_iter()
        .filter(|&l| l > 0 && l <= max_l)
        .map(|l| {
            let (n, o) = controls.gc_table(l);
            GcTable { l, n, o }
        })
        .collect();
    let totals = controls.region_totals();
    let regions = controls
        .regions
        .iter()
        .zip(totals)
        .map(|(r, obs)| RegionOut {
            name: format!("{}:{}-{}", r.chrom, r.start, r.end),
            role: r.role.clone(),
            label: r.label.clone(),
            len: r.len(),
            gc: (r.mean_gc() * 1e4).round() / 1e4,
            obs,
            absent: r.absent,
        })
        .collect();
    let classes = panel
        .classes
        .iter()
        .zip(acc.classes)
        .map(|(c, a)| ClassOut {
            name: c.name.clone(),
            kind: c.kind.as_str().to_string(),
            length: c.length,
            circular: c.circular,
            panel_kmers: c.n_kmers_kept,
            reads: a.reads,
            bases: a.bases,
            dup_flagged: a.dup_flagged,
            out_of_range: a.out_of_range,
            bin: p.bin,
            fwd: a.fwd,
            rev: a.rev,
            gc_read: a.gc_read,
            hit_frac: a.hit_frac,
        })
        .collect();
    let mut placements: Vec<Placement> = acc
        .placements
        .iter()
        .map(|(&(class, tid, b), &reads)| Placement {
            class: panel.classes[class as usize].name.clone(),
            contig: if tid < 0 { "*".to_string() } else { header_contigs[tid as usize].0.clone() },
            start: if tid < 0 { 0 } else { b as i64 * p.place_bin },
            reads,
        })
        .collect();
    // total order (ties broken by contig and position): the file is reproducible byte for byte
    placements
        .sort_by(|a, b| (&a.class, std::cmp::Reverse(a.reads), &a.contig, a.start).cmp(&(&b.class, std::cmp::Reverse(b.reads), &b.contig, b.start)));
    let mut ih = vec![0u64; INSERT_MAX / 10 + 1];
    for (i, &v) in st.insert.iter().enumerate() {
        ih[i / 10] += v;
    }
    // scan: every contig that holds a read; always: the contigs the controls and sinks sit on
    let scanned = mode == "scan";
    let contigs = header_contigs
        .iter()
        .enumerate()
        .filter_map(|(tid, (name, len))| {
            let reads = st.contig_reads.get(tid).copied().unwrap_or(0);
            let used = controls.regions.iter().any(|r| &r.chrom == name) || sink_contigs.iter().any(|c| c == name);
            ((scanned && reads > 0) || used).then(|| ContigOut { name: name.clone(), len: *len, reads: scanned.then_some(reads) })
        })
        .collect();
    Output {
        format: "ngs-dose-counts/1",
        engine_version: env!("CARGO_PKG_VERSION"),
        sample,
        input: input.path.clone(),
        mode: mode.to_string(),
        k: panel.k,
        min_hits: p.min_hits,
        min_frac: p.min_frac,
        panel: panel_path,
        panel_sha256: hashes.0,
        controls_sha256: hashes.1,
        sinks_sha256: hashes.2,
        controls: controls_path,
        sinks: sinks_path,
        unmapped_fetched: p.unmapped,
        records: st.records,
        primary: st.primary,
        primary_dup_flagged: st.primary_dup,
        unmapped: st.unmapped,
        ctrl_reads: st.ctrl,
        ctrl_dup_flagged: st.ctrl_dup,
        ctrl_mapq0: st.ctrl_mapq0,
        ctrl_mapq_lt20: st.ctrl_mapq_lt20,
        ctrl_positions: controls.regions.iter().filter(|r| r.role == "control").map(|r| r.len() as u64).sum(),
        read_length_mode: read_mode,
        insert_median,
        insert_hist_10bp: ih,
        readlen_hist: st.readlen.iter().enumerate().filter(|(_, &v)| v > 0).map(|(i, &v)| (i, v)).collect(),
        gc_tables,
        regions,
        classes,
        ambiguous_reads: acc.ambiguous,
        below_threshold_reads: acc.below_threshold,
        placement_bin: p.place_bin,
        placements,
        eof_marker,
        contigs,
        elapsed_sec: (elapsed * 100.0).round() / 100.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::{encode_ascii, revcomp_ascii};
    use crate::panel::ClassDef;

    fn toy_panel(unit: &[u8], k: usize) -> Panel {
        let mut map = FxHashMap::default();
        for km in KmerIter::new(&encode_ascii(unit), k) {
            map.insert(km.canon, Entry::new(0, km.offset, km.fwd_is_canon));
        }
        Panel::from_parts(
            k,
            vec![ClassDef {
                name: "u".into(),
                kind: Kind::Positional,
                length: unit.len(),
                circular: false,
                source: "t".into(),
                n_kmers_input: map.len(),
                n_kmers_kept: map.len(),
            }],
            map,
            vec![],
        )
    }

    fn params() -> Params {
        Params { min_hits: 4, min_frac: 0.0, bin: 10, l_grid: vec![], threads: 1, place_bin: 10000, unmapped: false, pad: 0, retries: 1 }
    }

    // pseudo-random but fixed unit so that all k-mers are unique
    fn unit() -> Vec<u8> {
        let mut x = 0x2545F4914F6CDD1Du64;
        (0..400)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    }

    fn one(panel: &Panel, s: &mut Scratch) -> Assignment {
        let mut out = Vec::new();
        classify(panel, &params(), s, &mut out);
        assert_eq!(out.len(), 1, "expected exactly one assignment");
        out[0]
    }

    #[test]
    fn forward_and_reverse_reads_recover_their_five_prime_end() {
        let u = unit();
        let panel = toy_panel(&u, 15);
        let mut s = Scratch::new(1);
        // forward read starting at 120
        s.codes = encode_ascii(&u[120..220]);
        let a = one(&panel, &mut s);
        assert_eq!((a.pos5, a.reverse, a.hits), (120, false, 100 - 15 + 1));
        // reverse-strand read covering 120..220: its 5' end is unit coordinate 219
        s.codes = encode_ascii(&revcomp_ascii(&u[120..220]));
        let a = one(&panel, &mut s);
        assert_eq!((a.pos5, a.reverse), (219, true));
    }

    #[test]
    fn mismatch_and_softclip_do_not_move_the_five_prime_end() {
        let u = unit();
        let panel = toy_panel(&u, 15);
        let mut s = Scratch::new(1);
        let mut read = u[50..150].to_vec();
        read[40] = if read[40] == b'A' { b'C' } else { b'A' }; // one SNV
        for b in read.iter_mut().take(8) {
            *b = b'N'; // unalignable 5' tail: the diagonal still points at unit position 50
        }
        s.codes = encode_ascii(&read);
        let a = one(&panel, &mut s);
        assert_eq!((a.pos5, a.reverse), (50, false));
    }

    #[test]
    fn hits_in_short_broken_runs_are_not_missed() {
        // two runs of three panel k-mers, neither containing an offset divisible by 4: the
        // strided screen of the first engine dropped this read although it holds 6 >= 4 hits
        let u = unit();
        let k = 15;
        let mut map = FxHashMap::default();
        for km in KmerIter::new(&encode_ascii(&u), k) {
            if matches!(km.offset, 101..=103 | 105..=107) {
                map.insert(km.canon, Entry::new(0, km.offset, km.fwd_is_canon));
            }
        }
        let cd = ClassDef {
            name: "u".into(),
            kind: Kind::Positional,
            length: u.len(),
            circular: false,
            source: "t".into(),
            n_kmers_input: 6,
            n_kmers_kept: 6,
        };
        let panel = Panel::from_parts(k, vec![cd], map, vec![]);
        let mut s = Scratch::new(1);
        s.codes = encode_ascii(&u[100..200]);
        let a = one(&panel, &mut s);
        assert_eq!((a.pos5, a.hits), (100, 6));
    }

    #[test]
    fn a_class_is_judged_on_its_own_hits() {
        // unit A (positional) followed by family B (compositional): a read across the junction is a
        // fragment end of A *and* mass of B, and A's verdict is the same with or without B loaded
        let u = unit();
        let k = 15;
        let (mut both, mut only_a) = (FxHashMap::default(), FxHashMap::default());
        for km in KmerIter::new(&encode_ascii(&u), k) {
            if km.offset < 200 {
                both.insert(km.canon, Entry::new(0, km.offset, km.fwd_is_canon));
                only_a.insert(km.canon, Entry::new(0, km.offset, km.fwd_is_canon));
            } else if km.offset >= 215 {
                both.insert(km.canon, Entry::new(1, 0, km.fwd_is_canon));
            }
        }
        let cd =
            |n: &str, kind| ClassDef { name: n.into(), kind, length: 200, circular: false, source: "t".into(), n_kmers_input: 0, n_kmers_kept: 0 };
        let p2 = Panel::from_parts(k, vec![cd("A", Kind::Positional), cd("B", Kind::Compositional)], both, vec![]);
        let p1 = Panel::from_parts(k, vec![cd("A", Kind::Positional)], only_a, vec![]);
        let read = encode_ascii(&u[150..250]);
        let (mut s2, mut s1) = (Scratch::new(2), Scratch::new(1));
        s2.codes = read.clone();
        s1.codes = read;
        let (mut o2, mut o1) = (Vec::new(), Vec::new());
        classify(&p2, &params(), &mut s2, &mut o2);
        classify(&p1, &params(), &mut s1, &mut o1);
        assert_eq!(o1.len(), 1);
        assert_eq!(o2.len(), 2, "fragment end of A and mass of B");
        let a2 = o2.iter().find(|a| a.class == 0).unwrap();
        assert_eq!((a2.pos5, a2.reverse, a2.hits), (o1[0].pos5, o1[0].reverse, o1[0].hits));
    }

    #[test]
    fn stored_orientation_and_record_strand_combine() {
        // a 100-base read covering unit [120, 220): forward 5' end = 120, reverse 5' end = 219
        let same = Assignment { class: 0, pos5: 120, reverse: false, hits: 50, valid: 86 };
        let rc = Assignment { class: 0, pos5: 219, reverse: true, hits: 50, valid: 86 };
        assert_eq!(same.sequenced(100, false), (120, false)); // aligned forward, stored as sequenced
        assert_eq!(same.sequenced(100, true), (219, true)); // aligned reverse: stored = revcomp(read)
        assert_eq!(rc.sequenced(100, false), (219, true)); // unmapped / sink copy in opposite orientation
        assert_eq!(rc.sequenced(100, true), (120, false));
    }

    #[test]
    fn unrelated_read_is_not_classified() {
        let u = unit();
        let panel = toy_panel(&u, 15);
        let mut s = Scratch::new(1);
        s.codes = encode_ascii(b"ACACACACACACACACACACACACACACACACACACACACACACACACACAC");
        let mut out = Vec::new();
        let un = classify(&panel, &params(), &mut s, &mut out);
        assert!(out.is_empty() && un == Unassigned::default());
    }
}
