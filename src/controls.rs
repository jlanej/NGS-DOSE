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
    /// `control` regions fit the GC curve and set the denominator; `test` and `dosage` regions
    /// (which need a label) are tallied and reported per region only - known-copy-number truth
    /// regions (chrX, held-out autosomal sequence) that the estimator is scored on, and the
    /// dosage of chrM and chrEBV. No other role is accepted.
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

/// Region roles: `control`, or `test` / `dosage` with a non-empty label. ngsdose/resources.py
/// (`Bundle.regions`) applies the same rule.
pub const ROLES: [&str; 3] = ["control", "test", "dosage"];

fn check_role(role: &str, label: &str) -> Result<()> {
    if !ROLES.contains(&role) {
        bail!("role={} is not one of control, test, dosage", role);
    }
    if role != "control" && label.is_empty() {
        bail!("role={} needs a label (label=...)", role);
    }
    Ok(())
}

fn parse_name(name: &str) -> Result<(String, i64, i64)> {
    let bad = || format!("control record name '{}' is not chrom:start-end", name);
    let (chrom, range) = name.rsplit_once(':').with_context(bad)?;
    let (s, e) = range.split_once('-').with_context(bad)?;
    let s: i64 = s.parse().with_context(|| format!("control record name '{}': bad start '{}'", name, s))?;
    let e: i64 = e.parse().with_context(|| format!("control record name '{}': bad end '{}'", name, e))?;
    if chrom.is_empty() || s < 0 || e <= s {
        bail!("{}: needs a contig and 0 <= start < end", bad());
    }
    Ok((chrom.to_string(), s, e))
}

/// A hint for contigs the header lacks: the same name with or without `chr`, or chrM / MT.
fn naming_hint(absent: &[String], tid_of: &dyn Fn(&str) -> Option<i32>) -> String {
    for c in absent {
        let alt = match c.as_str() {
            "chrM" => "MT".to_string(),
            "MT" => "chrM".to_string(),
            _ => match c.strip_prefix("chr") {
                Some(bare) => bare.to_string(),
                None => format!("chr{}", c),
            },
        };
        if tid_of(&alt).is_some() {
            return format!(" (the header has {} where the controls say {}: the controls and the file name their contigs differently)", alt, c);
        }
    }
    String::new()
}

impl Controls {
    pub fn load(path: &Path, tid_of: &dyn Fn(&str) -> Option<i32>) -> Result<Controls> {
        let mut regions = Vec::new();
        let mut by_tid: FxHashMap<i32, Vec<(i64, i64, usize)>> = FxHashMap::default();
        let mut missing: Vec<String> = Vec::new();
        let at = |name: &str| format!("{}: control record {}", path.display(), name);
        fasta::for_each_record(path, |name, desc, seq| {
            let (chrom, start, end) = parse_name(name).with_context(|| path.display().to_string())?;
            let flank = desc.split_whitespace().find_map(|t| t.strip_prefix("flank=")).with_context(|| format!("{} lacks flank=", at(name)))?;
            let flank: usize = flank.parse().with_context(|| format!("{}: bad flank={}", at(name), flank))?;
            let tag = |key: &str| desc.split_whitespace().find_map(|t| t.strip_prefix(key)).map(|v| v.to_string());
            let role = tag("role=").unwrap_or_else(|| "control".to_string());
            let label = tag("label=").unwrap_or_default();
            check_role(&role, &label).with_context(|| at(name))?;
            if seq.len() != (end - start) as usize + 2 * flank {
                bail!("{}: sequence length {} != region + 2*flank", at(name), seq.len());
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
                    if role == "control" && !missing.contains(&chrom) {
                        missing.push(chrom.clone());
                    }
                    true
                }
            };
            regions.push(Region { chrom, start, end, flank, role, label, absent, gc_cum, n_cum });
            Ok(())
        })?;
        if regions.is_empty() {
            bail!("no control regions in {}", path.display());
        }
        let n_control = regions.iter().filter(|r| r.role == "control").count();
        if n_control == 0 {
            bail!("{}: no region with role=control, so there is nothing to fit the GC curve on", path.display());
        }
        if !missing.is_empty() {
            let n = regions.iter().filter(|r| r.role == "control" && r.absent).count();
            let shown: Vec<&str> = missing.iter().take(5).map(|s| s.as_str()).collect();
            bail!(
                "{} of {} control regions of {} are on contigs absent from the alignment header ({}{}){} - wrong reference build?",
                n,
                n_control,
                path.display(),
                shown.join(", "),
                if missing.len() > 5 { format!(" and {} more", missing.len() - 5) } else { String::new() },
                naming_hint(&missing, tid_of)
            );
        }
        for v in by_tid.values_mut() {
            v.sort_unstable();
            for w in v.windows(2) {
                if w[1].0 < w[0].1 {
                    let (a, b) = (&regions[w[0].2], &regions[w[1].2]);
                    bail!(
                        "{}: control regions overlap: {}:{}-{} (role={}) and {}:{}-{} (role={})",
                        path.display(),
                        a.chrom,
                        a.start,
                        a.end,
                        a.role,
                        b.chrom,
                        b.start,
                        b.end,
                        b.role
                    );
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

    /// One interval per region whose contig is in the alignment header, padded by `pad` on both
    /// sides (start clamped at 0); not merged - count::fetch_plan merges them with the sinks.
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
        // BED name column: empty/`control`, or `test:<label>` / `dosage:<label>` such as `test:chrX`
        let (role, label) = match r.name.split_once(':') {
            Some((a, b)) => (a.to_string(), b.to_string()),
            None if r.name.is_empty() => ("control".to_string(), String::new()),
            None => (r.name.clone(), String::new()),
        };
        let checked = check_role(&role, &label).and_then(|_| {
            if role == "control" && !label.is_empty() {
                bail!("a control region takes no label");
            }
            if label.chars().any(char::is_whitespace) {
                bail!("the label must not contain whitespace");
            }
            Ok(())
        });
        checked.with_context(|| {
            format!(
                "{}: region {}:{}-{} named '{}' (the name column must be empty, control, test:<label> or dosage:<label>)",
                bed.display(),
                r.chrom,
                r.start,
                r.end,
                r.name
            )
        })?;
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

    fn tmpfile(name: &str, content: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("ngsdose_controls_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::write(&p, content).unwrap();
        p
    }

    /// A controls FASTA of 100-bp regions with flank 10 from `(name, description tail)` pairs.
    fn controls_fa(file: &str, recs: &[(&str, &str)]) -> std::path::PathBuf {
        let body: String = recs.iter().map(|(n, d)| format!(">{} flank=10 {}\n{}\n", n, d, "ACGT".repeat(30))).collect();
        tmpfile(file, &body)
    }

    fn load_err(path: &Path, contigs: &[&str]) -> String {
        let tid_of = |c: &str| contigs.iter().position(|x| *x == c).map(|i| i as i32);
        format!("{:#}", Controls::load(path, &tid_of).err().expect("load must fail"))
    }

    #[test]
    fn roles_are_checked() {
        let tid_of = |c: &str| (c == "chr1").then_some(0);
        let ok = controls_fa(
            "ok.fa",
            &[
                ("chr1:100-200", "role=control"),
                ("chr1:300-400", "role=test label=auto"),
                ("chr1:500-600", "role=dosage label=chrM"),
                ("chr1:700-800", ""),
            ],
        );
        let c = Controls::load(&ok, &tid_of).unwrap();
        assert_eq!(c.regions.iter().map(|r| r.role.as_str()).collect::<Vec<_>>(), vec!["control", "test", "dosage", "control"]);
        for (file, desc, needle) in [
            ("typo.fa", "role=contol", "role=contol is not one of"),
            ("case.fa", "role=Control", "role=Control is not one of"),
            ("nolabel.fa", "role=test", "role=test needs a label"),
        ] {
            let p = controls_fa(file, &[("chr1:100-200", "role=control"), ("chr1:300-400", desc)]);
            let e = load_err(&p, &["chr1"]);
            assert!(e.contains(needle) && e.contains("chr1:300-400") && e.contains(file), "{}", e);
        }
        let p = controls_fa("nocontrol.fa", &[("chr1:300-400", "role=test label=x")]);
        assert!(load_err(&p, &["chr1"]).contains("no region with role=control"));
    }

    #[test]
    fn load_errors_name_the_record() {
        let p = controls_fa("overlap.fa", &[("chr1:1000-1100", ""), ("chr1:1050-1150", "role=test label=x")]);
        let e = load_err(&p, &["chr1"]);
        assert!(e.contains("chr1:1000-1100 (role=control) and chr1:1050-1150 (role=test)"), "{}", e);
        let p = tmpfile("flank.fa", &format!(">chr1:100-200 flank=1O\n{}\n", "A".repeat(120)));
        let e = load_err(&p, &["chr1"]);
        assert!(e.contains("chr1:100-200") && e.contains("flank=1O") && e.contains("flank.fa"), "{}", e);
        let p = controls_fa("name.fa", &[("chr1:100-2O0", "")]);
        assert!(load_err(&p, &["chr1"]).contains("bad end '2O0'"));
        let p = controls_fa(
            "absent.fa",
            &[("chr1:100-200", ""), ("chr2:100-200", ""), ("chrM:100-200", ""), ("chrEBV:100-200", "role=dosage label=chrEBV")],
        );
        let e = load_err(&p, &["1", "2", "MT"]);
        assert!(
            e.contains("3 of 3 control regions") && e.contains("(chr1, chr2, chrM)") && e.contains("the header has 1 where the controls say chr1"),
            "{}",
            e
        );
        // an absent truth or dosage region is allowed
        let tid_of = |c: &str| (c != "chrEBV").then_some(0);
        let p = controls_fa("ebv.fa", &[("chr1:100-200", ""), ("chrEBV:100-200", "role=dosage label=chrEBV")]);
        assert!(Controls::load(&p, &tid_of).unwrap().regions[1].absent);
    }

    #[test]
    fn build_checks_the_name_column() {
        let reference = tmpfile("ref.fa", &format!(">chr1\n{}\n", "ACGTTGCA".repeat(100)));
        let out = std::env::temp_dir().join(format!("ngsdose_controls_test_{}", std::process::id())).join("out.fa.gz");
        let bed = tmpfile("ok.bed", "chr1\t100\t200\nchr1\t300\t400\tcontrol\nchr1\t500\t600\ttest:chrX\n");
        assert_eq!(build(&bed, &reference, 50, &out).unwrap(), 3);
        let tid_of = |c: &str| (c == "chr1").then_some(0);
        let c = Controls::load(&out, &tid_of).unwrap();
        assert_eq!(c.regions[2].role, "test");
        assert_eq!(c.regions[2].label, "chrX");
        for name in ["contols", "Control", "test", "test:", "dosage", "control:x", "region1"] {
            let bed = tmpfile("bad.bed", &format!("chr1\t100\t200\nchr1\t300\t400\t{}\n", name));
            let e = format!("{:#}", build(&bed, &reference, 50, &out).unwrap_err());
            assert!(e.contains(&format!("chr1:300-400 named '{}'", name)), "{}", e);
        }
    }

    #[test]
    fn n_in_window_is_skipped() {
        let seq = b"AAAAANAAAAGGGGGCCCCCAAAAAAAAAA".to_vec();
        let r = region_from(&seq, 10);
        assert_eq!(r.window_rev(0, 10), None); // window reaches back over the N
        assert!(r.window_fwd(0, 10).is_some());
    }
}
