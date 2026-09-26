//! k-mer panels: class-diagnostic k-mers with (class, unit position, strand).
//!
//! A *positional* class has a unit consensus (rDNA 45S, 5S, DJ ...): every retained k-mer occurs
//! exactly once in the unit, so a hit places the read on the unit. A *compositional* class is a
//! family without a stable unit (HSat, alpha satellite ...): k-mers only say "this read is
//! class X". In both cases a k-mer is retained only if it never occurs in the background
//! genomes outside the intervals where the class legitimately lives.

use crate::fasta;
use crate::kmer::{canonical, encode_ascii, kmer_to_string, string_to_kmer, KmerIter};
use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rustc_hash::{FxHashMap, FxHashSet};
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
        // a real assert: in a release build an out-of-range value would silently wrap the class
        // or flip the strand bit (entries are built only at panel build and load time)
        assert!(class < 256 && pos <= Self::MAX_POS, "panel entry out of range: class {} pos {}", class, pos);
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

#[derive(Debug)]
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

const MANIFEST_COLUMNS: &str = "name<TAB>kind<TAB>fasta<TAB>circular[<TAB>min_count[<TAB>keep.bed]]";

/// Parse a class manifest: tab-separated `name kind fasta circular [min_count] [keep.bed]`, with
/// kind positional|compositional, circular 0/1/true/false, min_count an integer >= 1 (empty or
/// '.' for the default 1; compositional classes only) and keep.bed a BED of unit intervals ('.'
/// or empty for none; positional classes only). Relative paths are taken from the manifest's
/// directory. Any other value is an error.
pub fn read_manifest(path: &Path) -> Result<Vec<ClassInput>> {
    let rd = fasta::open_maybe_gz(path)?;
    let base = path.parent().unwrap_or(Path::new("."));
    let resolve = |s: &str| {
        let p = PathBuf::from(s);
        if p.is_absolute() {
            p
        } else {
            base.join(p)
        }
    };
    let mut v: Vec<ClassInput> = Vec::new();
    let mut seen: FxHashMap<String, usize> = FxHashMap::default();
    for (i, line) in rd.lines().enumerate() {
        let line = line.with_context(|| format!("cannot read {}", path.display()))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let at = || format!("{}:{}", path.display(), i + 1);
        let p: Vec<&str> = line.split('\t').map(str::trim).collect();
        if p.len() < 4 || (p.len() > 6 && p[6..].iter().any(|s| !s.is_empty())) {
            bail!("{}: a manifest line is {}, got: {}", at(), MANIFEST_COLUMNS, line);
        }
        let name = p[0];
        if name.is_empty() {
            bail!("{}: empty class name", at());
        }
        if let Some(prev) = seen.insert(name.to_string(), i + 1) {
            bail!("{}: class {} is already defined on line {}", at(), name, prev);
        }
        let kind = Kind::parse(p[1]).with_context(at)?;
        if p[2].is_empty() {
            bail!("{}: class {} has no FASTA", at(), name);
        }
        let circular = match p[3].to_ascii_lowercase().as_str() {
            "0" | "false" => false,
            "1" | "true" => true,
            _ => bail!("{}: circular must be 0, 1, true or false, got '{}'", at(), p[3]),
        };
        let min_count = match p.get(4).copied().unwrap_or("") {
            "" | "." => 1,
            m => match m.parse::<u32>() {
                Ok(n) if n >= 1 => n,
                _ => bail!("{}: min_count must be an integer >= 1 (or '.' for the default 1), got '{}'", at(), m),
            },
        };
        if kind == Kind::Positional && min_count != 1 {
            bail!("{}: min_count applies to compositional classes only; a positional class keeps k-mers seen exactly once", at());
        }
        let keep = match p.get(5).copied().unwrap_or("") {
            "" | "." => None,
            k if kind == Kind::Positional => Some(resolve(k)),
            _ => bail!("{}: a keep BED applies to positional classes only", at()),
        };
        v.push(ClassInput { name: name.to_string(), kind, fasta: resolve(p[2]), circular, min_count, keep });
    }
    if v.is_empty() {
        bail!("no classes in manifest {}", path.display());
    }
    if v.len() > 255 {
        bail!("at most 255 classes per panel");
    }
    Ok(v)
}

#[derive(Debug)]
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
        /// first seen in its class outside that class's keep BED: left out of the panel, but kept
        /// here until every class is enumerated, so that a later class sharing it still marks it
        /// shared (the result does not depend on the manifest order)
        outside_keep: bool,
    }
    let mut names = FxHashSet::default();
    for inp in inputs {
        if !names.insert(inp.name.as_str()) {
            bail!("class {} is defined twice", inp.name);
        }
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
        if inp.kind == Kind::Positional && recs[0].2.len() > Entry::MAX_POS {
            bail!("unit {} too long ({} bp) for the panel position field", inp.name, recs[0].2.len());
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
                    // a k-mer that an earlier class left out through its keep BED now belongs to
                    // this class, shared: as if this class had come first in the manifest
                    Some(c) if c.outside_keep && c.class != ci => {
                        *c = Cand { class: ci, pos, fwd: km.fwd_is_canon, count: c.count + 1, shared: true, bg: 0, outside_keep: false };
                    }
                    Some(c) => {
                        if c.class != ci {
                            c.shared = true;
                        }
                        c.count += 1;
                    }
                    None => {
                        cand.insert(km.canon, Cand { class: ci, pos, fwd: km.fwd_is_canon, count: 1, shared: false, bg: 0, outside_keep: false });
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
                if s0 < 0 || e0 as usize > total {
                    bail!("{}: interval {}-{} is outside unit {} (0-{})", kp.display(), s0, e0, inp.name, total);
                }
                for slot in &mut keep[s0 as usize..e0 as usize] {
                    *slot = true;
                }
            }
            for c in cand.values_mut() {
                if c.class == ci && !keep[c.pos] {
                    c.outside_keep = true;
                }
            }
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
        let mask_path = || bg.mask.as_ref().map(|p| p.display().to_string()).unwrap_or_default();
        let mut seen: FxHashSet<String> = FxHashSet::default();
        eprintln!("[panel] scanning background {}", bg.fasta.display());
        fasta::for_each_record(&bg.fasta, |name, _d, seq| {
            if let Some(iv) = mask.get(name) {
                seen.insert(name.to_string());
                for &(s, e) in iv {
                    if s.max(0) as usize >= seq.len() {
                        bail!(
                            "mask {}: interval {}:{}-{} starts past the end of {} in {} ({} bp): a mask for another assembly?",
                            mask_path(),
                            name,
                            s,
                            e,
                            name,
                            bg.fasta.display(),
                            seq.len()
                        );
                    }
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
        if seen.len() < mask.len() {
            let mut absent: Vec<&str> = mask.keys().filter(|c| !seen.contains(*c)).map(|c| c.as_str()).collect();
            absent.sort_unstable();
            bail!(
                "mask {} names {} contig(s) that background {} does not have ({}{}): the class loci there would count as background (contig names with and without 'chr'? a mask for another assembly?)",
                mask_path(),
                absent.len(),
                bg.fasta.display(),
                absent.iter().take(5).cloned().collect::<Vec<_>>().join(", "),
                if absent.len() > 5 { ", ..." } else { "" }
            );
        }
    }
    if let Some(rp) = report {
        let f = std::fs::File::create(rp).with_context(|| format!("cannot create {}", rp.display()))?;
        let mut w = std::io::BufWriter::new(GzEncoder::new(f, Compression::new(6)));
        writeln!(w, "#class\tpos\tcount_in_class\tshared\tbackground")?;
        let mut rows: Vec<(usize, usize, u32, bool, u32)> =
            cand.values().filter(|c| !c.outside_keep).map(|c| (c.class, c.pos, c.count, c.shared, c.bg)).collect();
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
        if c.outside_keep {
            continue;
        }
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
        let (input, bg, kept) = (stats[i].1, stats[i].4, stats[i].5);
        if kept == 0 {
            eprintln!("[panel] WARNING: class {} keeps no k-mer ({} candidates, {} dropped as background)", c.name, input, bg);
        } else if 2 * bg > input {
            eprintln!(
                "[panel] WARNING: class {} loses {} of its {} k-mers to the background: check that the background masks cover the class's own loci",
                c.name, bg, input
            );
        }
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

    /// Load a panel written by `ngs-dose panel`. Every line is checked: k in 11..=32, `##class`
    /// ids in order and names unique, and each row a canonical k-mer of length k that occurs once,
    /// with a declared class, strand + or -, and pos within the unit (0 for a compositional class).
    pub fn load(path: &Path) -> Result<Panel> {
        let rd = fasta::open_maybe_gz(path)?;
        let mut k = 0usize;
        let mut classes: Vec<ClassDef> = Vec::new();
        let mut map: FxHashMap<u64, Entry> = FxHashMap::default();
        let mut header = Vec::new();
        let mut declared: Vec<Option<usize>> = Vec::new();
        for (i, line) in rd.lines().enumerate() {
            let line = line.with_context(|| format!("cannot read {}", path.display()))?;
            let at = || format!("{}:{}", path.display(), i + 1);
            if let Some(rest) = line.strip_prefix("##k=") {
                if k != 0 {
                    bail!("{}: a second ##k= line", at());
                }
                k = rest.trim().parse().with_context(|| format!("{}: bad ##k= line", at()))?;
                if !(11..=32).contains(&k) {
                    bail!("{}: k={} is outside 11..=32", at(), k);
                }
            } else if line.starts_with("##class\t") {
                let mut kv: FxHashMap<&str, &str> = FxHashMap::default();
                for f in line.split('\t').skip(1) {
                    if let Some((a, b)) = f.split_once('=') {
                        kv.insert(a, b);
                    }
                }
                let get = |key: &str| -> Result<&str> { kv.get(key).copied().with_context(|| format!("{}: ##class line missing {}", at(), key)) };
                let num = |key: &str| -> Result<Option<usize>> {
                    kv.get(key).map(|v| v.parse().with_context(|| format!("{}: bad {}={}", at(), key, v))).transpose()
                };
                let id: usize = get("id")?.parse().with_context(|| format!("{}: bad id", at()))?;
                if id != classes.len() {
                    bail!("{}: ##class id={} where id={} comes next (class lines must be in id order)", at(), id, classes.len());
                }
                let name = get("name")?;
                if name.is_empty() || classes.iter().any(|c| c.name == name) {
                    bail!("{}: class name '{}' is empty or already used", at(), name);
                }
                if classes.len() == 255 {
                    bail!("{}: more than 255 classes", at());
                }
                classes.push(ClassDef {
                    name: name.to_string(),
                    kind: Kind::parse(get("kind")?).with_context(at)?,
                    length: get("length")?.parse().with_context(|| format!("{}: bad length", at()))?,
                    circular: match get("circular")? {
                        "0" => false,
                        "1" => true,
                        c => bail!("{}: circular={} (0 or 1)", at(), c),
                    },
                    source: get("source").unwrap_or("").to_string(),
                    n_kmers_input: num("kmers_input")?.unwrap_or(0),
                    n_kmers_kept: num("kmers_kept")?.unwrap_or(0),
                });
                declared.push(num("kmers_kept")?);
            } else if line.starts_with('#') {
                if line.starts_with("##") && !line.starts_with("##ngs-dose-panel") {
                    header.push(line);
                }
            } else if !line.is_empty() {
                let p: Vec<&str> = line.split('\t').collect();
                if p.len() != 4 {
                    bail!("{}: a panel row is kmer<TAB>class<TAB>pos<TAB>strand: {}", at(), line);
                }
                if k == 0 || classes.is_empty() {
                    bail!("{}: a k-mer row before the ##k= and ##class lines", at());
                }
                if p[0].len() != k {
                    bail!("{}: panel k-mer length {} != k {}", at(), p[0].len(), k);
                }
                let v = string_to_kmer(p[0]).with_context(|| format!("{}: bad k-mer {}", at(), p[0]))?;
                if v != canonical(v, k) {
                    bail!("{}: k-mer {} is not in canonical form (the smaller of it and its reverse complement)", at(), p[0]);
                }
                let class: usize = p[1].parse().with_context(|| format!("{}: bad class {}", at(), p[1]))?;
                let pos: usize = p[2].parse().with_context(|| format!("{}: bad pos {}", at(), p[2]))?;
                let c = classes.get(class).with_context(|| format!("{}: class {} is not declared ({} ##class lines)", at(), class, classes.len()))?;
                let pos_ok = match c.kind {
                    Kind::Positional => pos < c.length && pos <= Entry::MAX_POS,
                    Kind::Compositional => pos == 0,
                };
                if !pos_ok {
                    bail!("{}: pos {} is outside class {} ({} bp, {})", at(), pos, c.name, c.length, c.kind.as_str());
                }
                let fwd = match p[3] {
                    "+" => true,
                    "-" => false,
                    s => bail!("{}: strand '{}' (+ or -)", at(), s),
                };
                if map.insert(v, Entry::new(class, pos, fwd)).is_some() {
                    bail!("{}: k-mer {} occurs twice", at(), p[0]);
                }
            }
        }
        if k == 0 || classes.is_empty() || map.is_empty() {
            bail!("{} is not a valid ngs-dose panel", path.display());
        }
        let mut n = vec![0usize; classes.len()];
        for e in map.values() {
            n[e.class()] += 1;
        }
        for ((c, n), d) in classes.iter().zip(n).zip(declared) {
            if d.is_some_and(|d| d != n) {
                eprintln!(
                    "[panel] WARNING: {}: class {} has {} k-mer rows but its ##class line says kmers_kept={}",
                    path.display(),
                    c.name,
                    n,
                    c.n_kmers_kept
                );
            }
        }
        Ok(Panel::from_parts(k, classes, map, header))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn tmpfile(name: &str, content: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("ngsdose_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::File::create(&p).unwrap().write_all(content.as_bytes()).unwrap();
        p
    }

    fn rnd_seq(seed: u64, n: usize) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 33) as usize & 3] as char
            })
            .collect()
    }

    fn positional(name: &str, fasta: PathBuf, circular: bool, keep: Option<PathBuf>) -> ClassInput {
        ClassInput { name: name.into(), kind: Kind::Positional, fasta, circular, min_count: 1, keep }
    }

    /// k-mer -> (class name, pos, fwd), which does not depend on the class ids
    fn by_name(p: &Panel) -> std::collections::BTreeMap<u64, (String, usize, bool)> {
        p.map.iter().map(|(k, e)| (*k, (p.classes[e.class()].name.clone(), e.pos(), e.fwd()))).collect()
    }

    #[test]
    fn circular_units_wrap() {
        let k = 11;
        let unit = rnd_seq(0xC0FFEE, 60);
        let fa = tmpfile("circ.fa", &format!(">c\n{}\n", unit));
        let (lin, _) = build(k, &[positional("c", fa.clone(), false, None)], &[], 0, None).unwrap();
        let (circ, _) = build(k, &[positional("c", fa, true, None)], &[], 0, None).unwrap();
        assert_eq!(lin.map.len(), 60 - k + 1);
        assert_eq!(circ.map.len(), 60, "one k-mer per unit position");
        let mut pos: Vec<usize> = circ.map.values().map(|e| e.pos()).collect();
        pos.sort_unstable();
        assert_eq!(pos, (0..60).collect::<Vec<_>>());
        // the wrap k-mers start in the last k-1 bases and continue at the unit start
        let wrapped = format!("{}{}", unit, &unit[..k - 1]);
        for p in 60 - (k - 1)..60 {
            let v = crate::kmer::canonical(string_to_kmer(&wrapped[p..p + k]).unwrap(), k);
            assert_eq!(circ.map[&v].pos(), p);
            assert!(!lin.map.contains_key(&v));
        }
    }

    #[test]
    fn a_keep_bed_gives_the_same_panel_in_any_manifest_order() {
        // A and B share a 30-bp segment R that lies outside A's keep BED: R's 11-mers are shared,
        // so they must leave B too, whichever class comes first
        let k = 11;
        let r = rnd_seq(7, 30);
        let a = format!("{}{}{}", rnd_seq(1, 40), r, rnd_seq(2, 40));
        let b = format!("{}{}{}", rnd_seq(3, 40), r, rnd_seq(4, 40));
        let fa = tmpfile("keepA.fa", &format!(">A\n{}\n", a));
        let fb = tmpfile("keepB.fa", &format!(">B\n{}\n", b));
        let keep = tmpfile("keepA.bed", "A\t0\t25\n");
        let ia = || positional("A", fa.clone(), false, Some(keep.clone()));
        let ib = || positional("B", fb.clone(), false, None);
        let (ab, st_ab) = build(k, &[ia(), ib()], &[], 0, None).unwrap();
        let (ba, st_ba) = build(k, &[ib(), ia()], &[], 0, None).unwrap();
        assert_eq!(by_name(&ab), by_name(&ba));
        let r_kmers: Vec<u64> = KmerIter::new(&encode_ascii(r.as_bytes()), k).map(|m| m.canon).collect();
        assert!(r_kmers.iter().all(|v| !ab.map.contains_key(v)), "R's k-mers are shared, so not in B");
        assert!(ab.map.values().filter(|e| e.class() == 1).count() >= 70);
        assert!(ab.map.values().filter(|e| e.class() == 0).all(|e| e.pos() < 25), "A keeps only starts in its keep BED");
        assert_eq!(ab.map.values().filter(|e| e.class() == 0).count(), 25);
        let (sa, sb) = (&st_ab.per_class, &st_ba.per_class);
        assert_eq!((&sa[0], &sa[1]), (&sb[1], &sb[0]), "per-class statistics agree too");
        // a keep interval beyond the unit is an error
        let far = tmpfile("keepFar.bed", "A\t100\t200\n");
        let e = build(k, &[positional("A", fa.clone(), false, Some(far))], &[], 0, None).err().unwrap();
        assert!(format!("{:#}", e).contains("outside unit A"), "{:#}", e);
    }

    #[test]
    fn manifests_are_parsed_strictly() {
        let dir = std::env::temp_dir().join(format!("ngsdose_test_{}", std::process::id()));
        let ok = tmpfile(
            "m_ok.tsv",
            "# comment\nU\tpositional\tU.fa\t 1\nV\tpositional\t/abs/V.fa\tFALSE\t.\tcore.bed\nC\tcompositional\tC.fa\t0\t10\nD\tcompositional\tD.fa\ttrue\t\t.\n",
        );
        let m = read_manifest(&ok).unwrap();
        assert_eq!(
            m.iter().map(|c| (c.name.as_str(), c.circular, c.min_count)).collect::<Vec<_>>(),
            vec![("U", true, 1), ("V", false, 1), ("C", false, 10), ("D", true, 1)]
        );
        assert_eq!(m[0].fasta, dir.join("U.fa"));
        assert_eq!(m[1].fasta, PathBuf::from("/abs/V.fa"));
        assert_eq!(m[1].keep, Some(dir.join("core.bed")));
        assert!(m[3].keep.is_none());
        for (i, (line, needle)) in [
            ("U\tpositional\tU.fa\tyes", "circular must be"),
            ("U\tpositional\tU.fa\t2", "circular must be"),
            ("C\tcompositional\tC.fa\t0\t5x", "min_count must be"),
            ("C\tcompositional\tC.fa\t0\t10.0", "min_count must be"),
            ("C\tcompositional\tC.fa\t0\t0", "min_count must be"),
            ("U\tpositional\tU.fa\t1\t3", "compositional classes only"),
            ("C\tcompositional\tC.fa\t0\t1\tk.bed", "positional classes only"),
            ("U\tposicional\tU.fa\t1", "unknown class kind"),
            ("U\tpositional\tU.fa", "a manifest line is"),
            ("U\tpositional\tU.fa\t1\t1\t.\textra", "a manifest line is"),
            ("U\tpositional\tU.fa\t1\nU\tcompositional\tC.fa\t0", "already defined on line 2"),
        ]
        .iter()
        .enumerate()
        {
            let p = tmpfile(&format!("m_bad{}.tsv", i), &format!("#h\n{}\n", line));
            let e = format!("{:#}", read_manifest(&p).unwrap_err());
            assert!(e.contains(needle) && e.contains(&format!("m_bad{}.tsv:", i)), "{}: {}", line, e);
        }
        // a duplicate class name given to build directly
        let fa = tmpfile("dup.fa", &format!(">u\n{}\n", rnd_seq(9, 40)));
        assert!(build(11, &[positional("u", fa.clone(), false, None), positional("u", fa, false, None)], &[], 0, None).is_err());
    }

    #[test]
    fn background_masks_must_match_the_background() {
        let k = 11;
        let unit = rnd_seq(11, 50);
        let fa = tmpfile("mU.fa", &format!(">U\n{}\n", unit));
        let bg = tmpfile("mbg.fa", &format!(">chrO\n{}\n>chrU\n{}\n", rnd_seq(12, 200), unit));
        let run = |mask: &str, file: &str| {
            build(k, &[positional("U", fa.clone(), false, None)], &[Background { fasta: bg.clone(), mask: Some(tmpfile(file, mask)) }], 0, None)
        };
        let (p, _) = run("chrU\t0\t50\n", "mask_ok.bed").unwrap();
        assert_eq!(p.map.len(), 40);
        let (p, _) = run("chrU\t20\t80\n", "mask_straddle.bed").unwrap();
        assert_eq!(p.map.len(), 40 - (20 - k + 1), "a mask end past the contig end is clipped");
        let e = format!("{:#}", run("U\t0\t50\nchrQ\t0\t9\n", "mask_names.bed").err().unwrap());
        assert!(e.contains("mask_names.bed names 2 contig(s)") && e.contains("(U, chrQ)") && e.contains("mbg.fa"), "{}", e);
        let e = format!("{:#}", run("chrU\t60\t70\n", "mask_past.bed").err().unwrap());
        assert!(e.contains("chrU:60-70 starts past the end of chrU") && e.contains("(50 bp)"), "{}", e);
    }

    #[test]
    fn malformed_panels_are_refused() {
        let k = 11;
        let unit = rnd_seq(21, 40);
        let (panel, _) = build(k, &[positional("A", tmpfile("lA.fa", &format!(">A\n{}\n", unit)), false, None)], &[], 0, None).unwrap();
        let good = tmpfile("lgood.tsv.gz", "");
        panel.write(&good).unwrap();
        let text = String::from_utf8(
            std::fs::read(&good)
                .map(|b| {
                    let mut s = Vec::new();
                    flate2::read::MultiGzDecoder::new(&b[..]).read_to_end(&mut s).unwrap();
                    s
                })
                .unwrap(),
        )
        .unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let first_row = lines.iter().position(|l| !l.starts_with('#')).unwrap();
        let row: Vec<&str> = lines[first_row].split('\t').collect();
        let with_row = |r: String| {
            let mut v: Vec<String> = lines.iter().map(|s| s.to_string()).collect();
            v[first_row] = r;
            v.join("\n") + "\n"
        };
        let rc = String::from_utf8(crate::kmer::revcomp_ascii(row[0].as_bytes())).unwrap();
        let class_b = "##class\tid=1\tname=B\tkind=compositional\tlength=10\tcircular=0\tsource=x\tkmers_input=0\tkmers_kept=0";
        let class_line = lines.iter().position(|l| l.starts_with("##class")).unwrap();
        let add_line = |at: usize, l: &str| {
            let mut v: Vec<String> = lines.iter().map(|s| s.to_string()).collect();
            v.insert(at, l.to_string());
            v.join("\n") + "\n"
        };
        let cases: Vec<(String, &str)> = vec![
            (with_row(format!("{}\t{}\t{}\t{}", rc, row[1], row[2], row[3])), "not in canonical form"),
            (with_row(format!("{}\t1\t{}\t{}", row[0], row[2], row[3])), "class 1 is not declared"),
            (with_row(format!("{}\t0\t40\t{}", row[0], row[3])), "pos 40 is outside class A"),
            (with_row(format!("{}\t0\t8388608\t{}", row[0], row[3])), "pos 8388608 is outside"),
            (with_row(format!("{}\t{}\t{}\tx", row[0], row[1], row[2])), "strand 'x'"),
            (with_row(format!("{}\t{}\t{}", row[0], row[1], row[2])), "a panel row is"),
            (add_line(first_row, lines[first_row + 1]), "occurs twice"),
            (text.replace("##k=11", "##k=33"), "outside 11..=32"),
            (add_line(class_line, &class_b.replace("id=1", "id=0").replace("name=B", "name=A")), "id=0 where id=1 comes next"),
            (add_line(class_line + 1, &class_b.replace("name=B", "name=A")), "already used"),
            (add_line(class_line + 1, &class_b.replace("circular=0", "circular=yes")), "circular=yes"),
        ];
        for (i, (content, needle)) in cases.iter().enumerate() {
            let p = tmpfile(&format!("lbad{}.tsv", i), content);
            let e = format!("{:#}", Panel::load(&p).err().unwrap_or_else(|| panic!("case {} ({}) loaded", i, needle)));
            assert!(e.contains(needle) && e.contains(&format!("lbad{}.tsv:", i)), "case {}: {}", i, e);
        }
        // the unchanged text still loads, as does a panel with a second, compositional class
        assert_eq!(Panel::load(&tmpfile("lsame.tsv", &text)).unwrap().map.len(), panel.map.len());
        let two = add_line(class_line + 1, class_b).replace("kmers_kept=0", "kmers_kept=1") + &format!("{}\t1\t0\t+\n", "A".repeat(k));
        let p = Panel::load(&tmpfile("ltwo.tsv", &two)).unwrap();
        assert_eq!(p.classes[1].name, "B");
        assert_eq!(p.map[&0].class(), 1);
    }

    #[test]
    #[should_panic(expected = "panel entry out of range")]
    fn entries_out_of_range_panic_in_release_too() {
        Entry::new(0, Entry::MAX_POS + 1, true);
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
