//! Single-copy control regions: fragment 5'-end counters and fragment-GC tables.
//!
//! The controls resource is a FASTA whose records are named `chrom:start-end` (0-based,
//! half-open) and carry `flank=F` in the description; each sequence is the region plus F
//! flanking bases on both sides, so the GC of a fragment-length window anchored at any
//! position of the region can be computed without the reference genome.

use crate::fasta;
use anyhow::{bail, Context, Result};
use rustc_hash::FxHashMap;
use std::io::Write;
use std::path::Path;
use std::sync::Mutex;

pub const GC_BINS: usize = 101;

pub struct Region {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub flank: usize,
    /// `control` regions fit the GC curve and set the denominator; any other role (e.g. `test`)
    /// is tallied and reported per region only - known-copy-number truth regions (chrX, held-out
    /// autosomal sequence) that the estimator is scored on
    pub role: String,
    pub label: String,
    /// the contig is not in this file's header (a reference without chrEBV, say); allowed for
    /// every role but `control`, reported, and never confused with "no reads"
    pub absent: bool,
    /// prefix sums over region+flanks: GC bases and non-ACGT bases
    gc_cum: Vec<u32>,
    n_cum: Vec<u32>,
}

#[derive(Default)]
pub struct RegionCounts {
    pub fwd: Vec<u16>,
    pub rev: Vec<u16>,
}

pub struct Controls {
    pub regions: Vec<Region>,
    pub counts: Vec<Mutex<RegionCounts>>,
    /// tid -> sorted (start, end, region index)
    by_tid: FxHashMap<i32, Vec<(i64, i64, usize)>>,
}

#[inline]
pub fn gc_bin(gc: u32, l: usize) -> usize {
    ((gc as usize * 100) + l / 2) / l
}

impl Region {
    pub fn len(&self) -> usize {
        (self.end - self.start) as usize
    }
    /// GC count of the window of length `l` starting at region offset `idx` (forward 5' end),
    /// or None if the window leaves the stored sequence or contains a non-ACGT base.
    #[inline]
    pub fn window_fwd(&self, idx: usize, l: usize) -> Option<u32> {
        let a = self.flank + idx;
        let b = a + l;
        if b >= self.gc_cum.len() || self.n_cum[b] != self.n_cum[a] {
            return None;
        }
        Some(self.gc_cum[b] - self.gc_cum[a])
    }
    /// Same for a reverse-strand 5' end at region offset `idx` (inclusive): window [idx-l+1, idx].
    #[inline]
    pub fn window_rev(&self, idx: usize, l: usize) -> Option<u32> {
        let b = self.flank + idx + 1;
        if b < l {
            return None;
        }
        let a = b - l;
        if self.n_cum[b] != self.n_cum[a] {
            return None;
        }
        Some(self.gc_cum[b] - self.gc_cum[a])
    }
    pub fn mean_gc(&self) -> f64 {
        let a = self.flank;
        let b = a + self.len();
        (self.gc_cum[b] - self.gc_cum[a]) as f64 / self.len() as f64
    }
}

fn parse_name(name: &str) -> Result<(String, i64, i64)> {
    let (chrom, range) = name.rsplit_once(':').with_context(|| format!("control record name '{}' is not chrom:start-end", name))?;
    let (s, e) = range.split_once('-').with_context(|| format!("control record name '{}' is not chrom:start-end", name))?;
    Ok((chrom.to_string(), s.parse()?, e.parse()?))
}

impl Controls {
    pub fn load(path: &Path, tid_of: &dyn Fn(&str) -> Option<i32>) -> Result<Controls> {
        let mut regions = Vec::new();
        let mut by_tid: FxHashMap<i32, Vec<(i64, i64, usize)>> = FxHashMap::default();
        let mut missing = 0usize;
        fasta::for_each_record(path, |name, desc, seq| {
            let (chrom, start, end) = parse_name(name)?;
            let flank: usize = desc
                .split_whitespace()
                .find_map(|t| t.strip_prefix("flank="))
                .with_context(|| format!("control record {} lacks flank=", name))?
                .parse()?;
            let tag = |key: &str| desc.split_whitespace().find_map(|t| t.strip_prefix(key)).map(|v| v.to_string());
            let role = tag("role=").unwrap_or_else(|| "control".to_string());
            let label = tag("label=").unwrap_or_default();
            if seq.len() != (end - start) as usize + 2 * flank {
                bail!("control record {}: sequence length {} != region + 2*flank", name, seq.len());
            }
            let mut gc_cum = Vec::with_capacity(seq.len() + 1);
            let mut n_cum = Vec::with_capacity(seq.len() + 1);
            let (mut g, mut n) = (0u32, 0u32);
            gc_cum.push(0);
            n_cum.push(0);
            for &b in seq.iter() {
                match b {
                    b'G' | b'C' | b'g' | b'c' => g += 1,
                    b'A' | b'T' | b'a' | b't' => {}
                    _ => n += 1,
                }
                gc_cum.push(g);
                n_cum.push(n);
            }
            let absent = match tid_of(&chrom) {
                Some(tid) => {
                    by_tid.entry(tid).or_default().push((start, end, regions.len()));
                    false
                }
                None => {
                    missing += (role == "control") as usize;
                    true
                }
            };
            regions.push(Region { chrom, start, end, flank, role, label, absent, gc_cum, n_cum });
            Ok(())
        })?;
        if regions.is_empty() {
            bail!("no control regions in {}", path.display());
        }
        if missing > 0 {
            bail!("{} of {} control regions are on contigs absent from the alignment header - wrong reference build?", missing, regions.len());
        }
        for v in by_tid.values_mut() {
            v.sort_unstable();
            for w in v.windows(2) {
                if w[1].0 < w[0].1 {
                    bail!("control regions overlap");
                }
            }
        }
        let counts = regions.iter().map(|r| Mutex::new(RegionCounts { fwd: vec![0; r.len()], rev: vec![0; r.len()] })).collect();
        Ok(Controls { regions, counts, by_tid })
    }

    /// Region containing reference position `pos` on `tid`, if any.
    #[inline]
    pub fn find(&self, tid: i32, pos: i64) -> Option<usize> {
        let v = self.by_tid.get(&tid)?;
        let i = v.partition_point(|&(s, _, _)| s <= pos);
        if i == 0 {
            return None;
        }
        let (s, e, idx) = v[i - 1];
        (pos >= s && pos < e).then_some(idx)
    }

    /// Locate a fragment 5' end (reference coordinate, inclusive): (region, offset, is_control).
    #[inline]
    pub fn locate(&self, tid: i32, five_prime: i64) -> Option<(usize, usize, bool)> {
        let idx = self.find(tid, five_prime)?;
        Some((idx, (five_prime - self.regions[idx].start) as usize, self.regions[idx].role == "control"))
    }

    /// Record a located fragment 5' end on the given strand.
    #[inline]
    pub fn bump(&self, idx: usize, off: usize, reverse: bool) {
        let mut c = self.counts[idx].lock().unwrap();
        let slot = if reverse { &mut c.rev[off] } else { &mut c.fwd[off] };
        *slot = slot.saturating_add(1);
    }

    /// Pooled tables for window length `l`: (positions N[g], observed 5' ends O[g]).
    pub fn gc_table(&self, l: usize) -> (Vec<u64>, Vec<u64>) {
        let mut n = vec![0u64; GC_BINS];
        let mut o = vec![0u64; GC_BINS];
        for (r, c) in self.regions.iter().zip(self.counts.iter()) {
            if r.role != "control" {
                continue;
            }
            let c = c.lock().unwrap();
            for idx in 0..r.len() {
                if let Some(g) = r.window_fwd(idx, l) {
                    let b = gc_bin(g, l);
                    n[b] += 1;
                    o[b] += c.fwd[idx] as u64;
                }
                if let Some(g) = r.window_rev(idx, l) {
                    let b = gc_bin(g, l);
                    n[b] += 1;
                    o[b] += c.rev[idx] as u64;
                }
            }
        }
        (n, o)
    }

    pub fn region_totals(&self) -> Vec<u64> {
        self.counts
            .iter()
            .map(|c| {
                let c = c.lock().unwrap();
                c.fwd.iter().map(|&x| x as u64).sum::<u64>() + c.rev.iter().map(|&x| x as u64).sum::<u64>()
            })
            .collect()
    }

    /// Padded, merged fetch intervals per contig name.
    pub fn fetch_intervals(&self, pad: i64) -> Vec<(String, i64, i64)> {
        self.regions.iter().filter(|r| !r.absent).map(|r| (r.chrom.clone(), (r.start - pad).max(0), r.end + pad)).collect()
    }
}

/// Build the controls FASTA from a BED and an indexed reference FASTA.
pub fn build(bed: &Path, reference: &Path, flank: usize, out: &Path) -> Result<usize> {
    let recs = fasta::read_bed(bed)?;
    let rd = rust_htslib::faidx::Reader::from_path(reference).map_err(|e| anyhow::anyhow!("cannot open reference {}: {}", reference.display(), e))?;
    let f = std::fs::File::create(out)?;
    let mut w = std::io::BufWriter::new(flate2::write::GzEncoder::new(f, flate2::Compression::new(9)));
    let mut n = 0usize;
    for r in &recs {
        if r.start < flank as i64 {
            bail!("control region {}:{}-{} is closer than flank={} to the contig start", r.chrom, r.start, r.end, flank);
        }
        let (s, e) = (r.start as usize - flank, r.end as usize + flank);
        let seq = rd.fetch_seq_string(&r.chrom, s, e - 1).map_err(|e| anyhow::anyhow!("fetch {}:{}-{} failed: {}", r.chrom, r.start, r.end, e))?;
        if seq.len() != e - s {
            bail!("control region {}:{}-{} runs off the contig end", r.chrom, r.start, r.end);
        }
        // BED name column: empty/`control`, or `role:label` such as `test:chrX`
        let (role, label) = match r.name.split_once(':') {
            Some((a, b)) => (a.to_string(), b.to_string()),
            None if r.name.is_empty() => ("control".to_string(), String::new()),
            None => (r.name.clone(), String::new()),
        };
        if label.is_empty() {
            writeln!(w, ">{}:{}-{} flank={} role={}", r.chrom, r.start, r.end, flank, role)?;
        } else {
            writeln!(w, ">{}:{}-{} flank={} role={} label={}", r.chrom, r.start, r.end, flank, role, label)?;
        }
        for chunk in seq.as_bytes().chunks(80) {
            w.write_all(&chunk.to_ascii_uppercase())?;
            w.write_all(b"\n")?;
        }
        n += 1;
    }
    w.into_inner()?.finish()?;
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region_from(seq: &[u8], flank: usize) -> Region {
        let mut gc_cum = vec![0u32];
        let mut n_cum = vec![0u32];
        let (mut g, mut n) = (0, 0);
        for &b in seq {
            match b {
                b'G' | b'C' => g += 1,
                b'A' | b'T' => {}
                _ => n += 1,
            }
            gc_cum.push(g);
            n_cum.push(n);
        }
        Region {
            chrom: "c".into(),
            start: 100,
            end: 100 + (seq.len() - 2 * flank) as i64,
            flank,
            role: "control".into(),
            label: String::new(),
            absent: false,
            gc_cum,
            n_cum,
        }
    }

    #[test]
    fn windows_are_strand_symmetric() {
        //            flank      region       flank
        let seq = b"AAAAAAAAAA GGGGGCCCCC AAAAAAAAAA".iter().filter(|b| **b != b' ').cloned().collect::<Vec<u8>>();
        let r = region_from(&seq, 10);
        assert_eq!(r.len(), 10);
        assert_eq!(r.window_fwd(0, 10), Some(10)); // whole region
        assert_eq!(r.window_rev(9, 10), Some(10)); // same window seen from the reverse 5' end
        assert_eq!(r.window_fwd(5, 10), Some(5));
        assert_eq!(r.window_rev(4, 10), Some(5));
        assert_eq!(gc_bin(5, 10), 50);
        assert_eq!(gc_bin(434, 434), 100);
    }

    #[test]
    fn n_in_window_is_skipped() {
        let seq = b"AAAAANAAAAGGGGGCCCCCAAAAAAAAAA".to_vec();
        let r = region_from(&seq, 10);
        assert_eq!(r.window_rev(0, 10), None); // window reaches back over the N
        assert!(r.window_fwd(0, 10).is_some());
    }
}
