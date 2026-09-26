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
use crate::fasta::{self, BedRec};
use crate::kmer::{unpack_bam_seq, KmerIter};
use crate::panel::{Entry, Kind, Panel};
use anyhow::{bail, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

pub const INSERT_MAX: usize = 1500;
pub const READLEN_MAX: usize = 400;
/// Padding (bp) of the control regions in a fetch, shared by `count` and `plan`, which must agree
/// for a cut along the plan to reproduce the fetch. At least READLEN_MAX: see `fetch_plan`.
pub const DEFAULT_PAD: i64 = 600;
pub const DEFAULT_RETRIES: usize = 5;

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
    /// attempts per interval (fetch) and per open of a remote input before giving up
    pub retries: usize,
}

impl Params {
    /// The ranges the command line enforces, checked again for callers that build Params
    /// themselves: a zero bin divides by zero, a negative placement bin never ends.
    pub fn check(&self) -> Result<()> {
        if self.bin == 0
            || !(1..=MAX_PLACE_BIN).contains(&self.place_bin)
            || self.min_hits == 0
            || self.retries == 0
            || !(0.0..=1.0).contains(&self.min_frac)
        {
            bail!(
                "parameters out of range: bin {} (>= 1), place_bin {} (1-1000000000), min_hits {} (>= 1), retries {} (>= 1), min_frac {} (0-1)",
                self.bin,
                self.place_bin,
                self.min_hits,
                self.retries,
                self.min_frac
            );
        }
        Ok(())
    }
}

/// An error another attempt cannot cure (a missing file, a refused credential, a reference that
/// does not match the CRAM): never retried, and never exit status 75.
#[derive(Debug)]
pub struct Permanent(pub String);

impl std::fmt::Display for Permanent {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Permanent {}

/// Context on an error of a remote input that outlasted every retry: main exits with EX_TEMPFAIL
/// (75), as the stall watchdog does, so that the caller tries the sample again later.
#[derive(Debug)]
pub struct TempFail;

impl std::fmt::Display for TempFail {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("remote input still failing after every retry (exit status 75: try again later)")
    }
}

fn is_permanent(e: &anyhow::Error) -> bool {
    e.downcast_ref::<Permanent>().is_some()
}

/// Wait before attempt `attempt + 1`: 1, 2, 4, 8, 16, then 16 s.
fn backoff(attempt: usize) {
    std::thread::sleep(std::time::Duration::from_millis(500 << attempt.min(5)));
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

/// Mapped primary reads per placement bin of the alignment, by leftmost position: all of them and
/// those flagged duplicate. A depth tool (mosdepth) sees a bin's total; the class placements say
/// how much of that total is the class. Dense in scan mode, where every read of the file is
/// tallied by one thread; sparse in fetch mode, where only the fetched intervals are seen.
#[derive(Clone)]
pub enum BinTally {
    Dense(Vec<Vec<[u32; 2]>>),
    Sparse(FxHashMap<(i32, i32), [u32; 2]>),
}

impl Default for BinTally {
    fn default() -> Self {
        BinTally::Sparse(FxHashMap::default())
    }
}

impl BinTally {
    #[inline]
    fn bump(&mut self, tid: i32, bin: usize, dup: bool) {
        let c = match self {
            BinTally::Dense(v) => {
                if v.len() <= tid as usize {
                    v.resize(tid as usize + 1, Vec::new());
                }
                let t = &mut v[tid as usize];
                if t.len() <= bin {
                    t.resize(bin + 4096, [0, 0]);
                }
                &mut t[bin]
            }
            BinTally::Sparse(m) => m.entry((tid, bin as i32)).or_insert([0, 0]),
        };
        c[0] += 1;
        c[1] += dup as u32;
    }
    fn get(&self, tid: i32, bin: usize) -> [u32; 2] {
        match self {
            BinTally::Dense(v) => v.get(tid as usize).and_then(|t| t.get(bin)).copied().unwrap_or([0, 0]),
            BinTally::Sparse(m) => m.get(&(tid, bin as i32)).copied().unwrap_or([0, 0]),
        }
    }
    fn merge(&mut self, o: &BinTally) {
        let mut add = |tid: i32, bin: usize, c: [u32; 2]| match self {
            BinTally::Dense(v) => {
                if v.len() <= tid as usize {
                    v.resize(tid as usize + 1, Vec::new());
                }
                let t = &mut v[tid as usize];
                if t.len() <= bin {
                    t.resize(bin + 1, [0, 0]);
                }
                t[bin][0] += c[0];
                t[bin][1] += c[1];
            }
            BinTally::Sparse(m) => {
                let e = m.entry((tid, bin as i32)).or_insert([0, 0]);
                e[0] += c[0];
                e[1] += c[1];
            }
        };
        match o {
            BinTally::Dense(v) => {
                for (tid, t) in v.iter().enumerate() {
                    for (bin, c) in t.iter().enumerate().filter(|(_, c)| c[0] > 0) {
                        add(tid as i32, bin, *c);
                    }
                }
            }
            BinTally::Sparse(m) => {
                for (&(tid, bin), c) in m {
                    add(tid, bin as usize, *c);
                }
            }
        }
    }
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
    /// mapped primary reads per placement bin (width `place_bin`)
    pub bins: BinTally,
    pub place_bin: i64,
    /// scan mode: whether the stream, read to its end, ended on an end-of-file marker (None when
    /// the format has none, or in fetch mode, which does not read to the end)
    pub end_marker: Option<bool>,
}

impl ReadStats {
    pub fn new(place_bin: i64, dense: bool) -> Self {
        ReadStats {
            insert: vec![0; INSERT_MAX + 1],
            readlen: vec![0; READLEN_MAX + 1],
            bins: if dense { BinTally::Dense(Vec::new()) } else { BinTally::default() },
            place_bin,
            ..Default::default()
        }
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
        self.bins.merge(&o.bins);
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

/// Width of the bins in which the alignment positions of a class's reads are recorded. A
/// satellite family covers tens of megabases, so it gets ten times the width of a positional class.
pub const COMPOSITIONAL_BIN_FACTOR: i64 = 10;

/// Largest placement bin accepted (1 Gb, longer than any chromosome): a wider one would only
/// overflow once multiplied by COMPOSITIONAL_BIN_FACTOR.
pub const MAX_PLACE_BIN: i64 = 1_000_000_000;

#[inline]
fn place_width(p: &Params, kind: Kind) -> i64 {
    if kind == Kind::Positional {
        p.place_bin
    } else {
        p.place_bin * COMPOSITIONAL_BIN_FACTOR
    }
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
    // percent GC of the read, rounded; a read of length zero has none
    if let Some(pct) = ((gc as usize * 100) + read_len / 2).checked_div(read_len) {
        ca.gc_read[pct] += 1;
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
    let pb = if tid < 0 { -1 } else { (pos / place_width(p, cdef.kind)) as i32 };
    *acc.placements.entry((a.class as u8, tid, pb)).or_insert(0) += 1;
}

/// Fragment 5' end on the reference (inclusive coordinate), soft clips included. A hard clip
/// outside the soft clip is stepped over but not added: class reads are placed from the stored
/// sequence, which holds soft-clipped bases and not hard-clipped ones, and the two must agree.
#[inline]
fn five_prime(rec: &bam::Record) -> i64 {
    five_prime_of(rec.pos(), rec.is_reverse(), rec.raw_cigar())
}

#[inline]
fn five_prime_of(pos: i64, reverse: bool, cig: &[u32]) -> i64 {
    if reverse {
        let mut end = pos;
        for &c in cig {
            match c & 0xf {
                0 | 2 | 3 | 7 | 8 => end += (c >> 4) as i64,
                _ => {}
            }
        }
        end + soft_clip(cig.iter().rev()) - 1
    } else {
        pos - soft_clip(cig.iter())
    }
}

/// Length of the soft clip at the end of the CIGAR the ops start from, inside any hard clip.
#[inline]
fn soft_clip<'a>(mut ops: impl Iterator<Item = &'a u32>) -> i64 {
    match ops.find(|&&c| c & 0xf != 5) {
        Some(&c) if c & 0xf == 4 => (c >> 4) as i64,
        _ => 0,
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
        st.bins.bump(tid, (rec.pos().max(0) / st.place_bin) as usize, dup);
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
    /// the input is a CRAM (set from the probe)
    pub cram: bool,
}

/// What opening the input once tells: its header, its end-of-file marker, and its format.
pub struct Probe {
    pub header: bam::HeaderView,
    /// "present", "absent", or "unchecked" (an unseekable stream, or a format without one)
    pub eof_marker: &'static str,
    pub cram: bool,
    /// reads without a coordinate, from the index of a BAM (a CRAM index does not record them)
    pub unplaced_unmapped: Option<u64>,
}

/// The input as it may be written down: a URL loses its query, fragment and user info, which
/// is where signed URLs carry their credentials.
pub fn redact(s: &str) -> String {
    if !s.contains("://") {
        return s.to_string();
    }
    match url::Url::parse(s) {
        Ok(mut u) => {
            let query = u.query().is_some();
            if !query && u.fragment().is_none() && u.username().is_empty() && u.password().is_none() {
                return s.to_string();
            }
            u.set_query(None);
            u.set_fragment(None);
            let _ = u.set_username("");
            let _ = u.set_password(None);
            let mut t = u.to_string();
            if query {
                t.push_str("?<redacted>");
            }
            t
        }
        Err(_) => {
            let t = s.split(['?', '#']).next().unwrap_or("");
            match t.split_once("://") {
                Some((scheme, rest)) => {
                    let (auth, path) = rest.split_at(rest.find('/').unwrap_or(rest.len()));
                    format!("{}://{}{}", scheme, auth.rsplit_once('@').map(|(_, h)| h).unwrap_or(auth), path)
                }
                None => t.to_string(),
            }
        }
    }
}

/// Every URL in a message redacted as `redact` does, whatever form it takes there: htslib reports
/// the URL as the url crate normalised it (host in lower case, spaces escaped, ...), and joins an
/// input and its index with "##idx##".
pub fn scrub_urls(msg: &str) -> String {
    let mut out = String::with_capacity(msg.len());
    let mut rest = msg;
    while let Some(i) = rest.find("://") {
        let start = rest[..i].char_indices().rev().take_while(|&(_, c)| c.is_ascii_alphanumeric() || "+.-".contains(c)).last().map(|(j, _)| j);
        let end = rest[i..].find(char::is_whitespace).map_or(rest.len(), |e| i + e);
        match start {
            Some(s) => {
                out.push_str(&rest[..s]);
                let parts: Vec<String> = rest[s..end].split("##idx##").map(redact_token).collect();
                out.push_str(&parts.join("##idx##"));
                rest = &rest[end..];
            }
            None => {
                out.push_str(&rest[..i + 3]);
                rest = &rest[i + 3..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// One URL as it appears in a message: cut at its query or fragment, without user info. A closing
/// punctuation mark after it is kept.
fn redact_token(t: &str) -> String {
    let body = t.trim_end_matches([';', ':', ',', ')', '\'', '"']);
    let tail = &t[body.len()..];
    let Some((scheme, rest)) = body.split_once("://") else {
        return t.to_string();
    };
    let cut = rest.find(['?', '#']).unwrap_or(rest.len());
    let (addr, query) = (&rest[..cut], rest[cut..].contains('?'));
    let (auth, path) = addr.split_at(addr.find('/').unwrap_or(addr.len()));
    let host = auth.rsplit_once('@').map_or(auth, |(_, h)| h);
    format!("{}://{}{}{}{}", scheme, host, path, if query { "?<redacted>" } else { "" }, tail)
}

/// Whether an errno from opening a file or URL says it is not there or is refused (htslib maps
/// HTTP 404 and 410 to ENOENT, 401 and 407 to EPERM, 403 to EACCES).
fn gone_errno(n: i32) -> bool {
    matches!(std::io::Error::from_raw_os_error(n).kind(), std::io::ErrorKind::NotFound | std::io::ErrorKind::PermissionDenied)
}

/// Ask the server for `url` once: true when it answers that the file is not there or is refused.
fn remote_gone(url: &str) -> bool {
    // as htslib is given it (the url crate's form)
    let url = url::Url::parse(url).map_or(url.to_string(), |u| u.to_string());
    let Ok(c) = std::ffi::CString::new(url) else {
        return true;
    };
    clear_errno();
    let fp = unsafe { htslib::hopen(c.as_ptr(), c"r".as_ptr()) };
    if fp.is_null() {
        return std::io::Error::last_os_error().raw_os_error().is_some_and(gone_errno);
    }
    unsafe { htslib::hclose(fp) };
    false
}

#[cfg(any(target_os = "macos", target_os = "ios", target_os = "freebsd"))]
extern "C" {
    #[link_name = "__error"]
    fn errno_location() -> *mut std::os::raw::c_int;
}
#[cfg(any(target_os = "linux", target_os = "android"))]
extern "C" {
    #[link_name = "__errno_location"]
    fn errno_location() -> *mut std::os::raw::c_int;
}

/// errno = 0, so that what htslib leaves in it after a failed call is its own.
pub fn clear_errno() {
    #[cfg(any(target_os = "macos", target_os = "ios", target_os = "freebsd", target_os = "linux", target_os = "android"))]
    unsafe {
        *errno_location() = 0;
    }
}

impl Input {
    pub fn new(path: String, index: Option<String>, reference: Option<PathBuf>) -> Input {
        Input { path, index, reference, cram: false }
    }
    pub fn is_url(&self) -> bool {
        self.path.contains("://")
    }
    /// The input for messages and the counts file (see `redact`).
    pub fn display(&self) -> String {
        redact(&self.path)
    }
    /// A failed open. `errno` is read right after htslib gave up, with errno cleared before the call:
    /// it says why the input itself could not be opened (from the server's answer for a URL). An
    /// index that failed to load is judged apart, by asking the server for it again (`index_gone`),
    /// since htslib's index search leaves an errno of its own.
    fn open_failed(&self, what: &str, e: rust_htslib::errors::Error, errno: Option<i32>) -> anyhow::Error {
        let msg = format!("cannot open {}{}: {}", self.display(), what, scrub_urls(&e.to_string()));
        let gone = match e {
            rust_htslib::errors::Error::BamInvalidIndex { .. } => self.index_gone(),
            _ => errno.filter(|&n| gone_errno(n)).map(|_| String::new()),
        };
        match gone {
            // a missing file or a refused credential stays so; anything else on a network may pass
            Some(why) => anyhow::Error::new(Permanent(msg + &why)),
            None => anyhow::anyhow!(msg),
        }
    }
    /// After the input opened but its index did not: Some(why) when no attempt can cure that (the
    /// index is not there, or refused), None when the server may answer next time.
    fn index_gone(&self) -> Option<String> {
        match &self.index {
            Some(i) if i.contains("://") => remote_gone(i).then(|| format!(" ({} is not there, or refused)", redact(i))),
            Some(i) => Some(format!(" (the local index {} {})", i, if Path::new(i).exists() { "cannot be read" } else { "does not exist" })),
            None if self.is_url() => {
                // htslib looks for X.csi, X.bai or X.crai, with or without X's extension
                let (base, query) = match self.path.split_once('?') {
                    Some((b, q)) => (b, format!("?{}", q)),
                    None => (self.path.as_str(), String::new()),
                };
                let file = base.rsplit('/').next().unwrap_or("");
                let stem = file.rfind('.').map(|d| &base[..base.len() - file.len() + d]);
                let tried = std::iter::once(base).chain(stem).flat_map(|b| [".csi", ".bai", ".crai"].map(|ext| format!("{}{}{}", b, ext, query)));
                tried.collect::<Vec<_>>().iter().all(|u| remote_gone(u)).then(|| {
                    " (no index beside it: .csi, .bai and .crai, after the file name or in place of its extension, are not there, or refused)"
                        .to_string()
                })
            }
            None => None,
        }
    }
    /// s3:// and gs:// need htslib's S3/GCS plugins, which this build does not include.
    fn check_scheme(&self) -> Result<()> {
        let scheme = self.path.split("://").next().unwrap_or("");
        if ["s3", "s3+http", "s3+https", "gs", "gs+http", "gs+https"].contains(&scheme) {
            let c = std::ffi::CString::new(self.path.as_bytes())?;
            if unsafe { htslib::hisremote(c.as_ptr()) } == 0 {
                return Err(anyhow::Error::new(Permanent(format!(
                    "cannot open {}: this build of ngs-dose reads no {}:// URLs (htslib without its S3/GCS plugins). Use an https:// (signed) URL, a mounted path, or `ngs-dose plan` and a samtools cut",
                    self.display(),
                    scheme
                ))));
            }
        }
        Ok(())
    }
    /// One attempt at opening the input with its index. htslib does not check the header read of
    /// an indexed open, so a failed one (a 503 from a server) leaves a reader without a header:
    /// that is an open failure too, and may pass on another attempt.
    pub fn open_indexed(&self) -> Result<bam::IndexedReader> {
        self.check_scheme()?;
        let rd = if self.is_url() {
            let full = match &self.index {
                Some(i) => format!("{}##idx##{}", self.path, i),
                None => self.path.clone(),
            };
            let url = url::Url::parse(&full).map_err(|e| anyhow::Error::new(Permanent(format!("bad URL {}: {}", redact(&full), e))))?;
            clear_errno();
            bam::IndexedReader::from_url(&url)
        } else {
            clear_errno();
            match &self.index {
                Some(i) => bam::IndexedReader::from_path_and_index(&self.path, i),
                None => bam::IndexedReader::from_path(&self.path),
            }
        };
        let mut rd = rd.map_err(|e| {
            let errno = std::io::Error::last_os_error().raw_os_error();
            self.open_failed(" with its index", e, errno)
        })?;
        if rd.header().inner_ptr().is_null() {
            bail!("cannot read the header of {}", self.display());
        }
        if let Some(r) = &self.reference {
            rd.set_reference(r).map_err(|e| anyhow::anyhow!("set_reference: {}", e))?;
        }
        set_cram_fields(&mut rd);
        Ok(rd)
    }
    pub fn open_stream(&self, threads: usize) -> Result<bam::Reader> {
        self.check_scheme()?;
        let rd = if self.is_url() {
            let url = url::Url::parse(&self.path).map_err(|e| anyhow::Error::new(Permanent(format!("bad URL {}: {}", self.display(), e))))?;
            clear_errno();
            bam::Reader::from_url(&url)
        } else {
            clear_errno();
            bam::Reader::from_path(&self.path)
        };
        let mut rd = rd.map_err(|e| {
            let errno = std::io::Error::last_os_error().raw_os_error();
            self.open_failed("", e, errno)
        })?;
        if let Some(r) = &self.reference {
            rd.set_reference(r).map_err(|e| anyhow::anyhow!("set_reference: {}", e))?;
        }
        set_cram_fields(&mut rd);
        if threads > 1 {
            rd.set_threads(threads).map_err(|e| anyhow::anyhow!("set_threads: {}", e))?;
        }
        Ok(rd)
    }
    /// Run `f` up to `retries` times for a remote input (once for a local one), backing off in
    /// between: the probe's and the scan's opens. (The fetch workers open inside their per-interval
    /// retries, so that a failed reopen costs one attempt, not the run.) `f` is told whether this is the last attempt. A remote input still failing at
    /// the end gives up with TempFail, unless the error is Permanent.
    fn retry<T>(&self, what: &str, retries: usize, mut f: impl FnMut(bool) -> Result<T>) -> Result<T> {
        let tries = if self.is_url() { retries.max(1) } else { 1 };
        let mut attempt = 0;
        loop {
            match f(attempt + 1 >= tries) {
                Ok(v) => return Ok(v),
                Err(e) if is_permanent(&e) => return Err(e),
                Err(e) => {
                    attempt += 1;
                    if attempt >= tries {
                        return Err(if self.is_url() { e.context(TempFail) } else { e });
                    }
                    eprintln!("[count] retry {}/{} to {}: {}", attempt, tries, what, scrub_urls(&e.to_string()));
                    backoff(attempt);
                }
            }
        }
    }
    /// Read the header and check the end-of-file marker, opening with the index (fetch, plan) or
    /// as a stream (scan). For a URL, a failed end-of-file check is retried like a failed open,
    /// so that one lost range request does not turn the truncation guard off (a scan checks the
    /// marker again at the end of the stream). A URL to fetch from (`seek`) must answer range
    /// requests: htslib reports a lost one and a server without them alike ("cannot seek"), so two
    /// in a row (or one on the last try) are taken for a server without them, which no retry cures.
    pub fn probe(&self, indexed: bool, seek: bool, retries: usize) -> Result<Probe> {
        let mut unseekable = 0;
        self.retry("open the input", retries, |last| {
            let (header, eof, cram, unplaced_unmapped) = if indexed {
                let rd = self.open_indexed()?;
                let n = if is_cram(&rd) { None } else { Some(unsafe { htslib::hts_idx_get_n_no_coor(rd.index().inner_ptr()) }) };
                (rd.header().clone(), eof_code(&rd), is_cram(&rd), n)
            } else {
                let rd = self.open_stream(1)?;
                (rd.header().clone(), eof_code(&rd), is_cram(&rd), None)
            };
            if self.is_url() && seek && eof == 2 {
                unseekable += 1;
                if unseekable >= 2 || last {
                    return Err(anyhow::Error::new(Permanent(format!(
                        "{} does not answer range requests (twice in a row): a fetch needs them. Scan the file instead, or copy it \
                         (or a samtools cut along `ngs-dose plan`) and fetch from the copy",
                        self.display()
                    ))));
                }
                bail!("the end-of-file check of {} failed (a range request went unanswered)", self.display());
            }
            unseekable = 0;
            if self.is_url() && !last && eof < 0 {
                bail!("the end-of-file check of {} failed", self.display());
            }
            let eof_marker = match eof {
                1 => "present",
                0 => "absent",
                _ => "unchecked",
            };
            Ok(Probe { header, eof_marker, cram, unplaced_unmapped })
        })
    }
    /// Why a read of `contig` may have failed for good: a CRAM decoded against a -T FASTA whose
    /// sequence for it differs from the one the CRAM was made with (htslib reports that only as
    /// a truncated record). None when there is nothing to compare or the two agree. `missing`:
    /// whether a FASTA without the contig counts as a mismatch.
    fn reference_mismatch(&self, contig: &str, want_m5: Option<&String>, missing: bool) -> Option<String> {
        let (Some(fa), Some(want)) = (self.reference.as_ref().filter(|_| self.cram), want_m5) else {
            return None;
        };
        let got = match fasta_m5(fa, contig) {
            Ok(Some(m)) => m,
            Ok(None) => return missing.then(|| format!("the reference {} has no contig {}, which the CRAM uses", fa.display(), contig)),
            Err(_) => return None,
        };
        (!got.eq_ignore_ascii_case(want)).then(|| {
            format!(
                "the reference {} does not match the CRAM on {} (@SQ M5 {}, FASTA M5 {}): pass -T the FASTA the CRAM was made with",
                fa.display(),
                contig,
                want,
                got
            )
        })
    }
}

fn is_cram<R: Read>(rd: &R) -> bool {
    unsafe { (*rd.htsfile()).format.format == htslib::htsExactFormat_cram }
}

/// htslib's end-of-file check: 1 present, 0 absent, 2 unseekable, 3 no marker in this format,
/// negative when the check itself failed. Must be asked before any record is read.
fn eof_code<R: Read>(rd: &R) -> i32 {
    unsafe { htslib::hts_check_EOF(rd.htsfile()) }
}

/// The error htslib recorded on the reader's file handle since the last call, cleared. A CRAM reader
/// ends an interval as if at its end when a request of a remote file fails mid-read (no record,
/// no error), so every interval, and a scan, is checked for one.
fn take_io_error<R: Read>(rd: &R) -> Option<i32> {
    unsafe {
        let fp = &*rd.htsfile();
        let h = if fp.format.format == htslib::htsExactFormat_cram {
            htslib::cram_fd_get_fp(fp.fp.cram)
        } else if matches!(fp.format.compression, htslib::htsCompression_bgzf | htslib::htsCompression_gzip) && !fp.fp.bgzf.is_null() {
            (*fp.fp.bgzf).fp
        } else {
            return None;
        };
        if h.is_null() || (*h).has_errno == 0 {
            return None;
        }
        let n = (*h).has_errno;
        (*h).has_errno = 0;
        Some(n)
    }
}

fn io_error(n: i32) -> anyhow::Error {
    anyhow::anyhow!("read error: the input failed mid-read ({})", std::io::Error::from_raw_os_error(n))
}

/// After a stream was read to its end: whether it ended on an end-of-file marker. For an input
/// whose marker could not be checked up front (a server without range requests) this is
/// the only evidence that nothing was cut off. A CRAM decoded with several threads does not keep
/// it (htslib marks every end alike), so there it stays unknown.
fn end_marker<R: Read>(rd: &R, threaded: bool) -> Option<bool> {
    unsafe {
        let fp = &*rd.htsfile();
        if fp.format.format == htslib::htsExactFormat_cram {
            if threaded {
                return None;
            }
            match htslib::cram_eof(fp.fp.cram) {
                1 => Some(true),
                2 => Some(false),
                _ => None,
            }
        } else if fp.format.compression == htslib::htsCompression_bgzf && !fp.fp.bgzf.is_null() {
            Some((*fp.fp.bgzf).no_eof_block() == 0)
        } else {
            None
        }
    }
}

/// The M5 tags of the header's @SQ lines, by contig name.
pub fn sq_m5(header: &bam::HeaderView) -> FxHashMap<String, String> {
    let text = String::from_utf8_lossy(header.as_bytes()).to_string();
    let mut m = FxHashMap::default();
    for line in text.lines().filter(|l| l.starts_with("@SQ\t")) {
        let tag = |t: &str| line.split('\t').find_map(|f| f.strip_prefix(t)).map(str::to_string);
        if let (Some(sn), Some(m5)) = (tag("SN:"), tag("M5:")) {
            m.insert(sn, m5);
        }
    }
    m
}

/// MD5 of a contig of an indexed FASTA as the SAM specification defines it (upper case, no
/// white space), which is what a CRAM's @SQ M5 records. None when the FASTA lacks the contig.
fn fasta_m5(fa: &Path, contig: &str) -> Result<Option<String>> {
    let rd = rust_htslib::faidx::Reader::from_path(fa).map_err(|e| anyhow::anyhow!("{}", e))?;
    let len = rd.fetch_seq_len(contig);
    if len == 0 || len > i64::MAX as u64 {
        return Ok(None);
    }
    unsafe {
        let ctx = htslib::hts_md5_init();
        if ctx.is_null() {
            bail!("hts_md5_init");
        }
        let step = 1 << 24;
        let mut at = 0usize;
        while (at as u64) < len {
            let end = (at + step).min(len as usize) - 1;
            let mut s = rd.fetch_seq(contig, at, end).map_err(|e| anyhow::anyhow!("{}", e))?;
            s.retain(|b| (b'!'..=b'~').contains(b));
            s.make_ascii_uppercase();
            htslib::hts_md5_update(ctx, s.as_ptr().cast(), s.len() as _);
            at = end + 1;
        }
        let mut digest = [0u8; 16];
        let mut hex = [0 as std::os::raw::c_char; 33];
        htslib::hts_md5_final(digest.as_mut_ptr(), ctx);
        htslib::hts_md5_hex(hex.as_mut_ptr(), digest.as_ptr());
        htslib::hts_md5_destroy(ctx);
        Ok(Some(std::ffi::CStr::from_ptr(hex.as_ptr()).to_string_lossy().to_string()))
    }
}

/// A record the reader returned as read although htslib failed inside its fixed-size core (an I/O
/// error mid-record gives an empty record that rust-htslib passes on as Ok).
#[inline]
fn check_record(rec: &bam::Record) -> Result<()> {
    if rec.inner().l_data == 0 {
        bail!("read error: empty record (an I/O error inside the record)");
    }
    Ok(())
}

fn panicked(p: Box<dyn std::any::Any + Send>) -> anyhow::Error {
    let what = p.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| p.downcast_ref::<String>().cloned()).unwrap_or_default();
    anyhow::anyhow!("a worker thread panicked: {}", what)
}

/// What the class side needs of a record besides its sequence.
#[derive(Clone, Copy)]
struct ReadKey {
    tid: i32,
    pos: i64,
    dup: bool,
    /// the stored sequence is the reverse complement of the read as sequenced
    rev: bool,
    len: usize,
}

/// The one place where scan and fetch read a record for classification: they must not differ.
#[inline]
fn read_key(rec: &bam::Record) -> ReadKey {
    ReadKey { tid: rec.tid(), pos: rec.pos(), dup: rec.is_duplicate(), rev: rec.is_reverse() && !rec.is_unmapped(), len: rec.seq_len() }
}

/// Classify one read (4-bit BAM encoding) and record its assignments.
#[inline]
fn process_read(panel: &Panel, p: &Params, acc: &mut Acc, s: &mut Scratch, enc: &[u8], key: &ReadKey) {
    let gc = unpack_bam_seq(enc, key.len, &mut s.codes);
    let mut out = std::mem::take(&mut s.out);
    let un = classify(panel, p, s, &mut out);
    for a in &out {
        record_assignment(panel, p, acc, a, key.len, gc, key.dup, key.rev, key.tid, key.pos);
    }
    s.out = out;
    acc.ambiguous += un.ambiguous as u64;
    acc.below_threshold += un.below as u64;
}

struct Batch {
    packed: Vec<u8>,
    items: Vec<(u32, ReadKey)>, // offset into packed
}

impl Batch {
    fn new() -> Batch {
        Batch { packed: Vec::with_capacity(1 << 21), items: Vec::with_capacity(1 << 14) }
    }
}

fn classify_batch(panel: &Panel, p: &Params, b: &Batch, acc: &mut Acc, s: &mut Scratch) {
    for (off, key) in &b.items {
        let off = *off as usize;
        process_read(panel, p, acc, s, &b.packed[off..off + key.len.div_ceil(2)], key);
    }
}

/// Whole-file pass: every primary read is screened against the panel wherever it was aligned.
pub fn scan(input: &Input, panel: &Panel, controls: &Controls, p: &Params) -> Result<(Acc, ReadStats)> {
    p.check()?;
    let workers = p.threads.max(1);
    let decomp = p.threads.clamp(1, 8);
    let mut rd = input.retry("open the input", p.retries, |_| input.open_stream(decomp))?;
    let (tx, rx) = crossbeam_channel::bounded::<Batch>(workers * 2);
    let mut stats = ReadStats::new(p.place_bin, true);
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
        let mut last = (-1i32, -1i64);
        while let Some(r) = rd.read(&mut rec) {
            if let Err(e) = r.map_err(|e| anyhow::anyhow!("read error: {}", e)).and_then(|()| check_record(&rec)) {
                return Err(scan_failure(input, rd.header(), e, stats.records, last));
            }
            last = (rec.tid(), rec.pos());
            if !handle_record_meta(&rec, controls, &mut stats, None) {
                continue;
            }
            let key = read_key(&rec);
            if key.len < panel.k {
                continue;
            }
            batch.items.push((batch.packed.len() as u32, key));
            batch.packed.extend_from_slice(&rec.seq().encoded[..key.len.div_ceil(2)]);
            if batch.items.len() >= 16384 {
                tx.send(std::mem::replace(&mut batch, Batch::new())).ok();
            }
        }
        tx.send(batch).ok();
        drop(tx);
        let accs: Result<Vec<Acc>> = handles.into_iter().map(|h| h.join().map_err(panicked)).collect();
        if let Some(n) = take_io_error(&rd) {
            return Err(scan_failure(input, rd.header(), io_error(n), stats.records, last));
        }
        accs
    })?;
    stats.end_marker = end_marker(&rd, decomp > 1);
    let mut acc = Acc::new(panel, p.bin);
    for a in &accs {
        acc.merge(a);
    }
    Ok((acc, stats))
}

/// A scan read error, located (records read so far, the last good record) so that the bad block
/// of a large file can be found, and explained when the CRAM's reference does not match.
fn scan_failure(input: &Input, header: &bam::HeaderView, e: anyhow::Error, records: u64, last: (i32, i64)) -> anyhow::Error {
    let name = |tid: i32| String::from_utf8_lossy(header.tid2name(tid as u32)).to_string();
    let at = if last.0 >= 0 { format!("{}:{}", name(last.0), last.1 + 1) } else { "none, or unmapped".to_string() };
    let e = e.context(format!("reading {} failed after {} records (last good record: {})", input.display(), records, at));
    if input.cram && input.reference.is_some() {
        // the contig of the last good record and the next; before any, every contig (in order)
        let m5 = sq_m5(header);
        let n = header.target_count() as i32;
        let tids: Vec<i32> = if last.0 >= 0 { (last.0..n.min(last.0 + 2)).collect() } else { (0..n).collect() };
        for tid in tids {
            if let Some(why) = input.reference_mismatch(&name(tid), m5.get(&name(tid)), last.0 >= 0) {
                return e.context(Permanent(why));
            }
        }
    }
    e
}

/// Sink intervals whose contig the alignment header lacks, per class: they cannot be fetched.
#[derive(Serialize, Default, Debug, Clone, PartialEq)]
pub struct SkippedSinks {
    pub intervals: u64,
    pub bp: u64,
}

/// A sinks BED filtered against the alignment header.
pub struct Sinks {
    pub kept: Vec<BedRec>,
    pub skipped: BTreeMap<String, SkippedSinks>,
    /// the BED carries a class column
    pub named: bool,
}

impl Sinks {
    pub fn intervals(&self) -> Vec<(String, i64, i64)> {
        self.kept.iter().map(|r| (r.chrom.clone(), r.start, r.end)).collect()
    }
    /// "DJ 2 (1.5 kb), rDNA45S 1 (40.0 kb)"
    pub fn skipped_summary(&self) -> String {
        self.skipped
            .iter()
            .map(|(c, s)| format!("{} {} ({:.1} kb)", if c.is_empty() { "(no class)" } else { c }, s.intervals, s.bp as f64 / 1e3))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

/// Read a sinks BED and keep the intervals on contigs `tid_of` knows, recording the rest per class.
pub fn load_sinks(path: &Path, tid_of: &dyn Fn(&str) -> Option<i32>) -> Result<Sinks> {
    let recs = fasta::read_bed(path)?;
    let named = recs.iter().any(|r| !r.name.is_empty());
    let (mut kept, mut skipped) = (Vec::new(), BTreeMap::<String, SkippedSinks>::new());
    for r in recs {
        if tid_of(&r.chrom).is_some() {
            kept.push(r);
        } else {
            let s = skipped.entry(r.name.clone()).or_default();
            s.intervals += 1;
            s.bp += (r.end - r.start) as u64;
        }
    }
    Ok(Sinks { kept, skipped, named })
}

/// The intervals a fetch reads: the control regions padded by `pad` on both sides (a fetch finds a
/// read by its leftmost aligned base but counts it at its fragment 5' end, and only when that 5'
/// end lies inside a region: a reverse-strand read starts up to its aligned reference span plus
/// trailing soft clip upstream of its 5' end, and a forward read with a leading soft clip starts
/// that clip's length downstream of it, so `pad` must be at least the longest read's reference
/// span including clips; reads found only through the extra pad are not counted), plus the class
/// sinks, merged where they touch or overlap, sorted by contig name and start. `ngs-dose plan`
/// writes this as BED for sites that must cut the reads out of a CRAM with samtools first.
pub fn fetch_plan(controls: &Controls, sinks: &[(String, i64, i64)], pad: i64) -> Vec<(String, i64, i64)> {
    let mut by_chrom: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
    for (c, s, e) in controls.fetch_intervals(pad).into_iter().chain(sinks.iter().cloned()) {
        by_chrom.entry(c).or_default().push((s, e));
    }
    let mut plan: Vec<(String, i64, i64)> = Vec::new();
    for (c, mut v) in by_chrom {
        v.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::new();
        for (s, e) in v {
            match merged.last_mut() {
                Some(l) if s <= l.1 => l.1 = l.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        plan.extend(merged.into_iter().map(|(s, e)| (c.clone(), s, e)));
    }
    plan.sort_unstable();
    plan
}

/// Targeted pass: only reads placed in the control regions and in the class sinks are
/// retrieved (plus, optionally, the unmapped bin). `header` is the input's, from the probe.
pub fn fetch(
    input: &Input,
    header: &bam::HeaderView,
    panel: &Panel,
    controls: &Controls,
    sinks: &[(String, i64, i64)],
    p: &Params,
) -> Result<(Acc, ReadStats)> {
    p.check()?;
    if p.pad < READLEN_MAX as i64 {
        bail!("pad {} is below {} bp: reverse-strand 5' ends at the start of control regions would be missed", p.pad, READLEN_MAX);
    }
    // merged, non-overlapping plan; each record is counted in the interval holding its start
    let mut plan: Vec<(i32, i64, i64, String)> = Vec::new();
    for (c, s, e) in fetch_plan(controls, sinks, p.pad) {
        let tid = match header.tid(c.as_bytes()) {
            Some(t) => t as i32,
            None => bail!("fetch interval contig '{}' is not in the alignment header", c),
        };
        plan.push((tid, s, e, c));
    }
    plan.sort_unstable();
    let m5 = sq_m5(header);
    let label = |job: Option<usize>| match job {
        Some(i) => format!("{}:{}-{}", plan[i].3, plan[i].1, plan[i].2),
        None => "the unmapped bin".to_string(),
    };
    let workers = p.threads.max(1).min(plan.len().max(1));
    let (tx, rx) = crossbeam_channel::unbounded::<Option<usize>>();
    for i in 0..plan.len() {
        tx.send(Some(i)).ok();
    }
    if p.unmapped {
        tx.send(None).ok();
    }
    drop(tx);
    // set by a worker that gives up, so that the others stop taking intervals
    let abort = AtomicBool::new(false);
    let results: Vec<Result<(Acc, ReadStats)>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                let rx = rx.clone();
                let (plan, m5, label, abort) = (&plan, &m5, &label, &abort);
                scope.spawn(move || -> Result<(Acc, ReadStats)> {
                    // opened on first use and again after any failure, each open one attempt
                    let mut rd: Option<bam::IndexedReader> = None;
                    let mut acc = Acc::new(panel, p.bin);
                    let mut st = ReadStats::new(p.place_bin, false);
                    let mut s = Scratch::new(panel.classes.len());
                    let mut rec = bam::Record::new();
                    let mut pending: Vec<(usize, usize, bool)> = Vec::new();
                    for job in rx.iter() {
                        if abort.load(Relaxed) {
                            break;
                        }
                        let iv = job.map(|i| (plan[i].0, plan[i].1, plan[i].2));
                        let mut attempt = 0usize;
                        loop {
                            // an interval is committed only if it was read to its end
                            let mut iv_acc = Acc::new(panel, p.bin);
                            let mut iv_st = ReadStats::new(p.place_bin, false);
                            pending.clear();
                            let res = (|| -> Result<()> {
                                if rd.is_none() {
                                    rd = Some(input.open_indexed()?);
                                }
                                let r = rd.as_mut().expect("opened above");
                                read_interval(r, iv, panel, controls, p, &mut iv_acc, &mut iv_st, &mut pending, &mut s, &mut rec)
                            })();
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
                                    rd = None;
                                    attempt += 1;
                                    let why = match job {
                                        Some(i) if attempt == 1 && !is_permanent(&e) => {
                                            input.reference_mismatch(&plan[i].3, m5.get(&plan[i].3), true)
                                        }
                                        _ => None,
                                    };
                                    if let Some(why) = why {
                                        abort.store(true, Relaxed);
                                        return Err(e.context(format!("interval {} failed", label(job))).context(Permanent(why)));
                                    }
                                    if is_permanent(&e) || attempt >= p.retries.max(1) {
                                        abort.store(true, Relaxed);
                                        let e = e.context(format!("interval {} failed {} times", label(job), attempt));
                                        return Err(if input.is_url() && !is_permanent(&e) { e.context(TempFail) } else { e });
                                    }
                                    if abort.load(Relaxed) {
                                        return Ok((acc, st)); // another worker gave up; its error is the one reported
                                    }
                                    eprintln!("[count] retry {}/{} for interval {}: {}", attempt, p.retries, label(job), scrub_urls(&e.to_string()));
                                    backoff(attempt);
                                }
                            }
                        }
                    }
                    Ok((acc, st))
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap_or_else(|x| Err(panicked(x)))).collect()
    });
    let mut acc = Acc::new(panel, p.bin);
    let mut stats = ReadStats::new(p.place_bin, false);
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
    take_io_error(rd);
    match job {
        Some((tid, a, b)) => rd.fetch((tid, a, b)),
        None => rd.fetch(bam::FetchDefinition::Unmapped),
    }
    .map_err(|e| anyhow::anyhow!("fetch failed: {}", e))?;
    while let Some(r) = rd.read(rec) {
        r.map_err(|e| anyhow::anyhow!("read error: {}", e))?;
        check_record(rec)?;
        if let Some((_, a, b)) = job {
            if rec.pos() < a || rec.pos() >= b {
                continue; // belongs to a neighbouring interval, or starts before this one
            }
        }
        if !handle_record_meta(rec, controls, st, Some(pending)) {
            continue;
        }
        let key = read_key(rec);
        if key.len < panel.k {
            continue;
        }
        process_read(panel, p, acc, s, rec.seq().encoded, &key);
    }
    match take_io_error(rd) {
        Some(n) => Err(io_error(n)),
        None => Ok(()),
    }
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
    /// reads of the class placed in this bin (by leftmost position; unmapped reads that sit at
    /// their mate's position are included, as an index query would return them)
    pub reads: u32,
    /// every mapped primary read in the bin, and those of them flagged duplicate: what a depth
    /// tool sees there, with and without its duplicate filter
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dup: Option<u32>,
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
    /// the commit the engine was built from (NGSDOSE_BUILD at compile time; "dev" for a local build)
    pub engine_build: &'static str,
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
    /// fetch mode: loaded panel classes for which the sinks BED had no interval, counted only where
    /// their reads fall inside other intervals (only with --allow-missing-sinks; otherwise refused)
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sinks_missing_classes: Vec<String>,
    /// fetch mode: sink intervals left out, per class, because their contig is not in the alignment
    /// header. A class whose every interval was left out is in sinks_missing_classes as well.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub sinks_skipped: BTreeMap<String, SkippedSinks>,
    /// fetch mode: bp added to both sides of every control region
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pad: Option<i64>,
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
    pub placement_bin_compositional: i64,
    pub placements: Vec<Placement>,
    pub eof_marker: &'static str,
    pub contigs: Vec<ContigOut>,
    /// which pipeline made the input: sinks are specific to an aligner and a reference
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pipeline: Option<Pipeline>,
    pub elapsed_sec: f64,
}

/// One @PG line of the alignment header.
#[derive(Serialize)]
pub struct PgLine {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pn: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vn: Option<String>,
}

/// The programs and the reference behind an input, from its header. Two inputs made against the
/// same reference sequences have the same `sq_sha256`, whatever their contig order; a masked
/// reference with the same contig lengths differs in it wherever the header carries M5 tags.
#[derive(Serialize)]
pub struct Pipeline {
    /// the @PG lines in header order
    pub pg: Vec<PgLine>,
    /// sha256 of the @SQ lines as "name<TAB>length<TAB>M5" (M5 empty where absent), sorted, each
    /// ending in a newline
    pub sq_sha256: String,
    /// @SQ lines, and how many of them carry M5
    pub sq_n: usize,
    pub sq_m5: usize,
}

pub fn pipeline(header: &bam::HeaderView) -> Option<Pipeline> {
    use sha2::{Digest, Sha256};
    let text = String::from_utf8_lossy(header.as_bytes()).to_string();
    let (mut pg, mut sq) = (Vec::new(), Vec::new());
    let mut sq_m5 = 0;
    for line in text.lines() {
        let tag = |t: &str| line.split('\t').skip(1).find_map(|f| f.strip_prefix(t)).map(str::to_string);
        if line.starts_with("@PG\t") {
            pg.push(PgLine { id: tag("ID:"), pn: tag("PN:"), vn: tag("VN:") });
        } else if line.starts_with("@SQ\t") {
            let m5 = tag("M5:");
            sq_m5 += m5.is_some() as usize;
            sq.push(format!("{}\t{}\t{}\n", tag("SN:").unwrap_or_default(), tag("LN:").unwrap_or_default(), m5.unwrap_or_default()));
        }
    }
    if pg.is_empty() && sq.is_empty() {
        return None;
    }
    sq.sort_unstable();
    let mut h = Sha256::new();
    for l in &sq {
        h.update(l.as_bytes());
    }
    let sq_sha256 = h.finalize().iter().map(|b| format!("{:02x}", b)).collect();
    Some(Pipeline { pg, sq_sha256, sq_n: sq.len(), sq_m5 })
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

/// The most frequent read length (the shorter on a tie); 0 when no control read was seen.
fn read_length_mode(readlen: &[u64]) -> usize {
    // the last of equal maxima, as before; only an all-empty histogram changes (0, not the top bin)
    readlen.iter().enumerate().filter(|(_, &v)| v > 0).max_by_key(|&(_, &v)| v).map(|(i, _)| i).unwrap_or(0)
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
    let read_mode = read_length_mode(&st.readlen);
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
        .map(|(&(class, tid, b), &reads)| {
            let width = place_width(p, panel.classes[class as usize].kind);
            // the tally is kept at the positional width; a wider bin is the sum of its parts
            let parts = (width / p.place_bin) as usize;
            let census = (tid >= 0)
                .then(|| (0..parts).map(|k| st.bins.get(tid, b as usize * parts + k)).fold([0u32, 0u32], |a, c| [a[0] + c[0], a[1] + c[1]]));
            Placement {
                class: panel.classes[class as usize].name.clone(),
                contig: if tid < 0 { "*".to_string() } else { header_contigs[tid as usize].0.clone() },
                start: if tid < 0 { 0 } else { b as i64 * width },
                reads,
                all: census.map(|c| c[0]),
                dup: census.map(|c| c[1]),
            }
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
        engine_build: option_env!("NGSDOSE_BUILD").unwrap_or("dev"),
        sample,
        input: input.display(),
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
        sinks_missing_classes: Vec::new(),
        sinks_skipped: BTreeMap::new(),
        pad: (mode == "fetch").then_some(p.pad),
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
        placement_bin_compositional: p.place_bin * COMPOSITIONAL_BIN_FACTOR,
        placements,
        eof_marker,
        contigs,
        pipeline: None,
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

    fn cigar(ops: &str) -> Vec<u32> {
        let mut v = Vec::new();
        let mut n = 0u32;
        for c in ops.chars() {
            match c.to_digit(10) {
                Some(d) => n = n * 10 + d,
                None => {
                    v.push(n << 4 | "MIDNSHP=X".find(c).unwrap() as u32);
                    n = 0;
                }
            }
        }
        v
    }

    #[test]
    fn five_prime_ends_on_the_reference() {
        // (CIGAR, reverse, 5' end) for a record whose leftmost aligned base is 1000
        for (ops, rev, want) in [
            ("100M", false, 1000),
            ("100M", true, 1099),
            ("10S90M", false, 990),
            ("10S90M", true, 1089),
            ("90M10S", true, 1099),
            ("40M5D55M", false, 1000),
            ("40M5D55M", true, 1099),
            ("40M5I55M", true, 1094),
            ("40M100N30=1X29M", false, 1000),
            ("40M100N30=1X29M", true, 1199),
            // a hard clip outside the soft clip is stepped over, and only the soft clip counts
            ("5H10S85M", false, 990),
            ("85M10S5H", true, 1094),
            ("5H95M", false, 1000),
            ("95M5H", true, 1094),
        ] {
            assert_eq!(five_prime_of(1000, rev, &cigar(ops)), want, "{} reverse={}", ops, rev);
        }
    }

    #[test]
    fn params_out_of_range_are_refused() {
        assert!(params().check().is_ok());
        assert!(Params { place_bin: MAX_PLACE_BIN, ..params() }.check().is_ok());
        let bad = [
            Params { bin: 0, ..params() },
            Params { place_bin: 0, ..params() },
            Params { place_bin: -1000, ..params() },
            Params { place_bin: MAX_PLACE_BIN + 1, ..params() },
            Params { place_bin: i64::MAX, ..params() },
            Params { min_frac: 2.0, ..params() },
            Params { min_frac: f64::NAN, ..params() },
            Params { min_hits: 0, ..params() },
            Params { retries: 0, ..params() },
        ];
        for p in bad {
            assert!(p.check().is_err());
        }
    }

    #[test]
    fn signed_urls_are_redacted() {
        assert_eq!(redact("/data/x.cram"), "/data/x.cram");
        assert_eq!(redact("https://host/x.cram"), "https://host/x.cram");
        assert_eq!(redact("https://b.s3.amazonaws.com/x.cram?X-Amz-Signature=SECRET&a=b"), "https://b.s3.amazonaws.com/x.cram?<redacted>");
        assert_eq!(redact("https://user:pw@host/x.cram#frag"), "https://host/x.cram");
        assert_eq!(redact("https://u:p@bad host/a@b.cram?sig=S"), "https://bad host/a@b.cram"); // not a valid URL: cut by hand
        let i = Input::new("https://host/x.bam?sig=SECRET".into(), Some("https://host/x.bam.csi?sig=SECRET".into()), None);
        assert!(!i.display().contains("SECRET"));
        // htslib reports the URL as the url crate normalised it, not as given
        for given in ["http://LOCALHOST:8080/x.bam?sig=SECRET1", "http://127.0.0.1:80/x.bam?sig=SE%20CRET2&x=a b", "https://u:SECRET3@HOST/x.bam#f"] {
            let norm = url::Url::parse(given).unwrap().to_string();
            let msg = format!("unable to open SAM/BAM/CRAM index for {}##idx##{}.csi; please create an index", norm, norm.replace("/x.bam", "/y"));
            let got = scrub_urls(&msg);
            assert!(!got.contains("SECRET") && got.ends_with("; please create an index"), "{}", got);
        }
        assert_eq!(
            scrub_urls("open http://h/m.bam?t=S: file at http://localhost:1/m.bam?sig=S&x=a%20b##idx##HTTPS://h/m.csi?t=S; please"),
            "open http://h/m.bam?<redacted>: file at http://localhost:1/m.bam?<redacted>##idx##HTTPS://h/m.csi?<redacted>; please"
        );
        assert_eq!(scrub_urls("no URL: a://, ://x and /data/x.bam?y"), "no URL: a://, ://x and /data/x.bam?y");
        assert_eq!(scrub_urls("(https://h/x.bam?s=1) and ftp://u@h/y"), "(https://h/x.bam?<redacted>) and ftp://h/y");
    }

    #[test]
    fn the_read_length_mode_is_zero_without_control_reads() {
        let mut h = vec![0u64; READLEN_MAX + 1];
        assert_eq!(read_length_mode(&h), 0);
        h[151] = 9;
        h[100] = 5;
        assert_eq!(read_length_mode(&h), 151);
        h[100] = 9;
        assert_eq!(read_length_mode(&h), 151); // a tie goes to the longer length, as in earlier builds
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
