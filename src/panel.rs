//! k-mer panels: class-diagnostic k-mers with (class, unit position, strand).
//!
//! A *positional* class has a unit consensus (rDNA 45S, 5S, DJ ...): every retained k-mer occurs
//! exactly once in the unit, so a hit places the read on the unit. A *compositional* class is a
//! family without a stable unit (HSat, alpha satellite ...): k-mers only say "this read is
//! class X". In both cases a k-mer is retained only if it never occurs in the background
//! genomes outside the intervals where the class legitimately lives.

use crate::fasta;
use crate::kmer::{encode_ascii, kmer_to_string, string_to_kmer, KmerIter};
use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rustc_hash::FxHashMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Positional,
    Compositional,
}

impl Kind {
    pub fn as_str(&self) -> &'static str {
        match self {
            Kind::Positional => "positional",
            Kind::Compositional => "compositional",
        }
    }
    pub fn parse(s: &str) -> Result<Kind> {
        match s {
            "positional" => Ok(Kind::Positional),
            "compositional" => Ok(Kind::Compositional),
            _ => bail!("unknown class kind '{}' (positional|compositional)", s),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ClassDef {
    pub name: String,
    pub kind: Kind,
    /// unit length for positional classes; total input bp for compositional ones
    pub length: usize,
    pub circular: bool,
    pub source: String,
    pub n_kmers_input: usize,
    pub n_kmers_kept: usize,
}

/// Packed panel entry: class (8 bits) | fwd flag (1 bit) | position (23 bits).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Entry(pub u32);

impl Entry {
    pub const MAX_POS: usize = (1 << 23) - 1;
    pub fn new(class: usize, pos: usize, fwd: bool) -> Entry {
        debug_assert!(class < 256 && pos <= Self::MAX_POS);
        Entry(((class as u32) << 24) | ((fwd as u32) << 23) | pos as u32)
    }
    #[inline]
    pub fn class(&self) -> usize {
        (self.0 >> 24) as usize
    }
    /// true if the canonical k-mer is the unit's forward-strand k-mer at `pos`
    #[inline]
    pub fn fwd(&self) -> bool {
        (self.0 >> 23) & 1 == 1
    }
    #[inline]
    pub fn pos(&self) -> usize {
        (self.0 & 0x7f_ffff) as usize
    }
}

pub struct Panel {
    pub k: usize,
    pub classes: Vec<ClassDef>,
    pub map: FxHashMap<u64, Entry>,
    pub header: Vec<String>,
    /// One-hash bitset over the panel k-mers, sized for a ~1.5% false-positive rate. Every k-mer
    /// of every read is tested against it, and only positives reach the hash map, so that
    /// classification is exact (no k-mer is skipped) at a few nanoseconds per k-mer.
    filter: Vec<u64>,
    filter_shift: u32,
}

pub struct ClassInput {
    pub name: String,
    pub kind: Kind,
    pub fasta: PathBuf,
    pub circular: bool,
    /// compositional classes: keep k-mers seen at least this many times in the class FASTA
    pub min_count: u32,
    /// positional classes: optional BED (unit coordinates) of k-mer start positions to retain,
    /// e.g. the core of a paralogous family that is shared by every copy
    pub keep: Option<PathBuf>,
}

pub struct Background {
    pub fasta: PathBuf,
    /// intervals of this background genome that belong to panel classes and must not count
    pub mask: Option<PathBuf>,
}

/// Parse a class manifest: tab-separated `name kind fasta circular(0/1) [min_count] [keep.bed]`.
pub fn read_manifest(path: &Path) -> Result<Vec<ClassInput>> {
    let rd = fasta::open_maybe_gz(path)?;
    let base = path.parent().unwrap_or(Path::new("."));
    let mut v = Vec::new();
    for line in rd.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 {
            bail!("manifest line needs: name<TAB>kind<TAB>fasta<TAB>circular[<TAB>min_count]: {}", line);
        }
        let fa = PathBuf::from(p[2]);
        v.push(ClassInput {
            name: p[0].to_string(),
            kind: Kind::parse(p[1])?,
            fasta: if fa.is_absolute() { fa } else { base.join(fa) },
            circular: p[3] == "1" || p[3].eq_ignore_ascii_case("true"),
            min_count: p.get(4).and_then(|s| s.parse().ok()).unwrap_or(1),
            keep: p.get(5).filter(|s| !s.is_empty() && **s != ".").map(|s| {
                let k = PathBuf::from(s);
                if k.is_absolute() {
                    k
                } else {
                    base.join(k)
                }
            }),
        });
    }
    if v.is_empty() {
        bail!("no classes in manifest {}", path.display());
    }
    if v.len() > 255 {
        bail!("at most 255 classes per panel");
    }
    Ok(v)
}

pub struct BuildStats {
    pub per_class: Vec<(String, usize, usize, usize, usize, usize)>, // name, input, multi, shared, background, kept
}

pub fn build(k: usize, inputs: &[ClassInput], backgrounds: &[Background], max_bg: u32, report: Option<&Path>) -> Result<(Panel, BuildStats)> {
    // 1. enumerate class k-mers
    struct Cand {
        class: usize,
        pos: usize,
        fwd: bool,
        count: u32,
        shared: bool,
        bg: u32,
    }
    let mut cand: FxHashMap<u64, Cand> = FxHashMap::default();
    let mut classes = Vec::new();
    for (ci, inp) in inputs.iter().enumerate() {
        let recs = fasta::read_all(&inp.fasta)?;
        if recs.is_empty() {
            bail!("{}: no sequences", inp.fasta.display());
        }
        if inp.kind == Kind::Positional && recs.len() != 1 {
            bail!("positional class {} needs exactly one FASTA record (the unit), got {}", inp.name, recs.len());
        }
        let mut total = 0usize;
        let mut source = Vec::new();
        for (name, _d, seq) in &recs {
            total += seq.len();
            if source.len() < 3 {
                source.push(name.clone());
            }
            let mut s = seq.clone();
            if inp.circular {
                let ext: Vec<u8> = seq.iter().take(k - 1).cloned().collect();
                s.extend_from_slice(&ext);
            }
            let codes = encode_ascii(&s);
            for km in KmerIter::new(&codes, k) {
                let pos = if inp.kind == Kind::Positional { km.offset } else { 0 };
                match cand.get_mut(&km.canon) {
                    Some(c) => {
                        if c.class != ci {
                            c.shared = true;
                        }
                        c.count += 1;
                    }
                    None => {
                        cand.insert(km.canon, Cand { class: ci, pos, fwd: km.fwd_is_canon, count: 1, shared: false, bg: 0 });
                    }
                }
            }
        }
        if let Some(kp) = &inp.keep {
            if inp.kind != Kind::Positional {
                bail!("class {}: a keep BED only applies to positional classes", inp.name);
            }
            let iv: Vec<(i64, i64)> = fasta::merged_by_chrom(&fasta::read_bed(kp)?).into_values().flatten().collect();
            let mut keep = vec![false; total];
            for (s0, e0) in iv {
                for slot in keep.iter_mut().take((e0 as usize).min(total)).skip(s0.max(0) as usize) {
                    *slot = true;
                }
            }
            cand.retain(|_, c| c.class != ci || keep[c.pos]);
        }
        if inp.kind == Kind::Positional && total > Entry::MAX_POS {
            bail!("unit {} too long ({} bp) for the panel position field", inp.name, total);
        }
        classes.push(ClassDef {
            name: inp.name.clone(),
            kind: inp.kind,
            length: total,
            circular: inp.circular,
            source: source.join(","),
            n_kmers_input: 0,
            n_kmers_kept: 0,
        });
    }
    // 2. background occurrences
    for bg in backgrounds {
        let mask = match &bg.mask {
            Some(p) => fasta::merged_by_chrom(&fasta::read_bed(p)?),
            None => FxHashMap::default(),
        };
        eprintln!("[panel] scanning background {}", bg.fasta.display());
        fasta::for_each_record(&bg.fasta, |name, _d, seq| {
            if let Some(iv) = mask.get(name) {
                for &(s, e) in iv {
                    let (s, e) = (s.max(0) as usize, (e as usize).min(seq.len()));
                    for b in &mut seq[s..e] {
                        *b = b'N';
                    }
                }
            }
            let codes = encode_ascii(seq);
            for km in KmerIter::new(&codes, k) {
                if let Some(c) = cand.get_mut(&km.canon) {
                    c.bg = c.bg.saturating_add(1);
                }
            }
            Ok(())
        })?;
    }
    if let Some(rp) = report {
        let f = std::fs::File::create(rp).with_context(|| format!("cannot create {}", rp.display()))?;
        let mut w = std::io::BufWriter::new(GzEncoder::new(f, Compression::new(6)));
        writeln!(w, "#class\tpos\tcount_in_class\tshared\tbackground")?;
        let mut rows: Vec<(usize, usize, u32, bool, u32)> = cand.values().map(|c| (c.class, c.pos, c.count, c.shared, c.bg)).collect();
        rows.sort_unstable();
        for (class, pos, count, shared, bg) in rows {
            writeln!(w, "{}\t{}\t{}\t{}\t{}", inputs[class].name, pos, count, shared as u8, bg)?;
        }
        w.into_inner()?.finish()?;
    }
    // 3. filter
    let mut stats: Vec<(String, usize, usize, usize, usize, usize)> = classes.iter().map(|c| (c.name.clone(), 0, 0, 0, 0, 0)).collect();
    let mut map: FxHashMap<u64, Entry> = FxHashMap::default();
    for (kmer, c) in cand.iter() {
        let st = &mut stats[c.class];
        st.1 += 1;
        let inp = &inputs[c.class];
        if c.shared {
            st.3 += 1;
            continue;
        }
        if inp.kind == Kind::Positional && c.count > 1 {
            st.2 += 1;
            continue;
        }
        if inp.kind == Kind::Compositional && c.count < inp.min_count {
            st.2 += 1;
            continue;
        }
        if c.bg > max_bg {
            st.4 += 1;
            continue;
        }
        st.5 += 1;
        map.insert(*kmer, Entry::new(c.class, c.pos, c.fwd));
    }
    for (i, c) in classes.iter_mut().enumerate() {
        c.n_kmers_input = stats[i].1;
        c.n_kmers_kept = stats[i].5;
    }
    let mut header = vec![format!("##max_background_count={}", max_bg)];
    for bg in backgrounds {
        header.push(format!(
            "##background={} mask={}",
            bg.fasta.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default(),
            bg.mask.as_ref().and_then(|p| p.file_name()).map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "none".into())
        ));
    }
    Ok((Panel::from_parts(k, classes, map, header), BuildStats { per_class: stats }))
}

impl Panel {
    pub fn from_parts(k: usize, classes: Vec<ClassDef>, map: FxHashMap<u64, Entry>, header: Vec<String>) -> Panel {
        let bits = (map.len().max(1) * 64).next_power_of_two().max(1 << 20);
        let shift = 64 - bits.trailing_zeros();
        let mut filter = vec![0u64; bits / 64];
        for &kmer in map.keys() {
            let h = (kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> shift) as usize;
            filter[h >> 6] |= 1u64 << (h & 63);
        }
        Panel { k, classes, map, header, filter, filter_shift: shift }
    }

    /// False means the k-mer is certainly not in the panel.
    #[inline]
    pub fn maybe(&self, kmer: u64) -> bool {
        let h = (kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.filter_shift) as usize;
        self.filter[h >> 6] & (1u64 << (h & 63)) != 0
    }

    #[inline]
    pub fn get(&self, kmer: u64) -> Option<Entry> {
        if self.maybe(kmer) {
            self.map.get(&kmer).copied()
        } else {
            None
        }
    }

    /// Load and merge several panels (e.g. the positional bundle panel and a satellite panel).
    /// k must agree, class names must be unique, and a k-mer claimed by two panels is dropped
    /// from both: it cannot say which class a read belongs to.
    pub fn load_many(paths: &[PathBuf]) -> Result<Panel> {
        let mut it = paths.iter();
        let first = it.next().context("no panel given")?;
        let mut merged = Panel::load(first)?;
        let mut shared: Vec<u64> = Vec::new();
        for path in it {
            let other = Panel::load(path)?;
            if other.k != merged.k {
                bail!("{}: k={} but {} has k={}", path.display(), other.k, first.display(), merged.k);
            }
            let offset = merged.classes.len();
            if offset + other.classes.len() > 255 {
                bail!("more than 255 classes after merging panels");
            }
            for c in &other.classes {
                if merged.classes.iter().any(|m| m.name == c.name) {
                    bail!("class {} occurs in more than one panel", c.name);
                }
            }
            merged.classes.extend(other.classes.iter().cloned());
            merged.header.extend(other.header.iter().cloned());
            for (kmer, e) in other.map {
                match merged.map.entry(kmer) {
                    std::collections::hash_map::Entry::Occupied(_) => shared.push(kmer),
                    std::collections::hash_map::Entry::Vacant(v) => {
                        v.insert(Entry::new(e.class() + offset, e.pos(), e.fwd()));
                    }
                }
            }
        }
        for kmer in &shared {
            merged.map.remove(kmer);
        }
        if !shared.is_empty() {
            eprintln!("[panel] {} k-mers occur in more than one panel and were dropped", shared.len());
        }
        let mut kept = vec![0usize; merged.classes.len()];
        for e in merged.map.values() {
            kept[e.class()] += 1;
        }
        for (c, n) in merged.classes.iter_mut().zip(kept) {
            c.n_kmers_kept = n;
        }
        Ok(Panel::from_parts(merged.k, merged.classes, merged.map, merged.header))
    }

    pub fn write(&self, path: &Path) -> Result<()> {
        let f = std::fs::File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
        let mut w = std::io::BufWriter::new(GzEncoder::new(f, Compression::new(6)));
        writeln!(w, "##ngs-dose-panel v1")?;
        writeln!(w, "##k={}", self.k)?;
        for h in &self.header {
            writeln!(w, "{}", h)?;
        }
        for (i, c) in self.classes.iter().enumerate() {
            writeln!(
                w,
                "##class\tid={}\tname={}\tkind={}\tlength={}\tcircular={}\tsource={}\tkmers_input={}\tkmers_kept={}",
                i,
                c.name,
                c.kind.as_str(),
                c.length,
                c.circular as u8,
                c.source,
                c.n_kmers_input,
                c.n_kmers_kept
            )?;
        }
        writeln!(w, "#kmer\tclass\tpos\tstrand")?;
        let mut rows: Vec<(usize, usize, u64, bool)> = self.map.iter().map(|(k, e)| (e.class(), e.pos(), *k, e.fwd())).collect();
        rows.sort_unstable();
        for (class, pos, kmer, fwd) in rows {
            writeln!(w, "{}\t{}\t{}\t{}", kmer_to_string(kmer, self.k), class, pos, if fwd { '+' } else { '-' })?;
        }
        w.into_inner()?.finish()?;
        Ok(())
    }

    pub fn load(path: &Path) -> Result<Panel> {
        let rd = fasta::open_maybe_gz(path)?;
        let mut k = 0usize;
        let mut classes: Vec<ClassDef> = Vec::new();
        let mut map: FxHashMap<u64, Entry> = FxHashMap::default();
        let mut header = Vec::new();
        for line in rd.lines() {
            let line = line?;
            if let Some(rest) = line.strip_prefix("##k=") {
                k = rest.trim().parse().context("bad ##k= line")?;
            } else if line.starts_with("##class\t") {
                let mut kv: FxHashMap<&str, &str> = FxHashMap::default();
                for f in line.split('\t').skip(1) {
                    if let Some((a, b)) = f.split_once('=') {
                        kv.insert(a, b);
                    }
                }
                let get = |key: &str| -> Result<&str> { kv.get(key).copied().with_context(|| format!("##class line missing {}", key)) };
                classes.push(ClassDef {
                    name: get("name")?.to_string(),
                    kind: Kind::parse(get("kind")?)?,
                    length: get("length")?.parse()?,
                    circular: get("circular")? == "1",
                    source: get("source").unwrap_or("").to_string(),
                    n_kmers_input: get("kmers_input").unwrap_or("0").parse().unwrap_or(0),
                    n_kmers_kept: get("kmers_kept").unwrap_or("0").parse().unwrap_or(0),
                });
            } else if line.starts_with('#') {
                if line.starts_with("##") && !line.starts_with("##ngs-dose-panel") {
                    header.push(line);
                }
            } else if !line.is_empty() {
                let p: Vec<&str> = line.split('\t').collect();
                if p.len() < 4 {
                    bail!("bad panel row: {}", line);
                }
                if p[0].len() != k {
                    bail!("panel k-mer length {} != k {}", p[0].len(), k);
                }
                let v = string_to_kmer(p[0]).with_context(|| format!("bad k-mer {}", p[0]))?;
                let class: usize = p[1].parse()?;
                let pos: usize = p[2].parse()?;
                map.insert(v, Entry::new(class, pos, p[3] == "+"));
            }
        }
        if k == 0 || classes.is_empty() || map.is_empty() {
            bail!("{} is not a valid ngs-dose panel", path.display());
        }
        Ok(Panel::from_parts(k, classes, map, header))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpfile(name: &str, content: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("ngsdose_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::File::create(&p).unwrap().write_all(content.as_bytes()).unwrap();
        p
    }

    #[test]
    fn build_filters_background_and_roundtrips() {
        // fixed pseudo-random unit: all 11-mers distinct, none shared with its own reverse strand
        let mut x = 0x9E3779B97F4A7C15u64;
        let unit: String = (0..60)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 33) as usize & 3] as char
            })
            .collect();
        let k = 11;
        let distinct: std::collections::HashSet<u64> = KmerIter::new(&encode_ascii(unit.as_bytes()), k).map(|m| m.canon).collect();
        assert_eq!(distinct.len(), unit.len() - k + 1, "test unit must have unique k-mers");
        let fa = tmpfile("unit.fa", &format!(">u\n{}\n", unit));
        // the background holds the first 20 bases of the unit, i.e. its first 10 k-mers
        let bg = tmpfile("bg.fa", &format!(">c\nTTTTTTTT{}GGGGGGGG\n", &unit[..20]));
        let inputs = vec![ClassInput { name: "u".into(), kind: Kind::Positional, fasta: fa, circular: false, min_count: 1, keep: None }];
        let (panel, st) = build(k, &inputs, &[Background { fasta: bg, mask: None }], 0, None).unwrap();
        let n_in = unit.len() - k + 1;
        assert_eq!(st.per_class[0].1, n_in);
        assert_eq!(st.per_class[0].4, 10);
        assert_eq!(panel.map.len(), n_in - 10);
        assert!(panel.map.values().all(|e| e.pos() >= 10), "background k-mers are the first ten positions");
        let out = tmpfile("p.tsv.gz", "");
        panel.write(&out).unwrap();
        let re = Panel::load(&out).unwrap();
        assert_eq!(re.k, k);
        assert_eq!(re.map.len(), panel.map.len());
        for (kmer, e) in &panel.map {
            assert_eq!(re.map[kmer], *e);
        }
    }

    #[test]
    fn merging_panels_offsets_classes_and_drops_shared_kmers() {
        // two units that share their first 30 bases: the 20 shared 11-mers must vanish from both
        let mut x = 0x1234_5678_9ABC_DEF1u64;
        let mut rnd = |n: usize| -> String {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x >> 33) as usize & 3] as char
                })
                .collect()
        };
        let (shared, ta, tb) = (rnd(30), rnd(50), rnd(50));
        let k = 11;
        let mk = |name: &str, seq: &str| {
            let fa = tmpfile(&format!("{name}.fa"), &format!(">{name}\n{seq}\n"));
            let inputs = vec![ClassInput { name: name.into(), kind: Kind::Positional, fasta: fa, circular: false, min_count: 1, keep: None }];
            let (panel, _) = build(k, &inputs, &[], 0, None).unwrap();
            let out = tmpfile(&format!("{name}.panel.gz"), "");
            panel.write(&out).unwrap();
            (out, panel.map.len())
        };
        let (pa, na) = mk("A", &format!("{shared}{ta}"));
        let (pb, nb) = mk("B", &format!("{shared}{tb}"));
        let merged = Panel::load_many(&[pa.clone(), pb.clone()]).unwrap();
        assert_eq!(merged.classes.iter().map(|c| c.name.as_str()).collect::<Vec<_>>(), vec!["A", "B"]);
        assert_eq!(merged.map.len(), na + nb - 2 * 20);
        assert_eq!(merged.classes[0].n_kmers_kept + merged.classes[1].n_kmers_kept, merged.map.len());
        assert!(merged.map.values().any(|e| e.class() == 1), "second panel's class ids are offset");
        assert!(merged.map.keys().all(|&kmer| merged.maybe(kmer)), "the prefilter has no false negatives");
        assert!(Panel::load_many(&[pa.clone(), pa]).is_err(), "duplicate class names are refused");
    }
}
