//! Minimal FASTA / BED readers (plain or gzip).

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut f = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    let f = File::open(path)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, f)))
    }
}

/// Stream FASTA records, calling `f(name, description, sequence)` for each.
pub fn for_each_record<F>(path: &Path, mut f: F) -> Result<()>
where
    F: FnMut(&str, &str, &mut Vec<u8>) -> Result<()>,
{
    let mut rd = open_maybe_gz(path)?;
    let mut line = Vec::with_capacity(256);
    let mut name = String::new();
    let mut desc = String::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut have = false;
    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            if have {
                f(&name, &desc, &mut seq)?;
            }
            let h = String::from_utf8_lossy(&line[1..]).to_string();
            let mut it = h.splitn(2, char::is_whitespace);
            name = it.next().unwrap_or("").to_string();
            desc = it.next().unwrap_or("").to_string();
            seq.clear();
            have = true;
        } else if have {
            seq.extend_from_slice(&line);
        } else if !line.is_empty() {
            bail!("{}: sequence data before first FASTA header", path.display());
        }
    }
    if have {
        f(&name, &desc, &mut seq)?;
    }
    Ok(())
}

pub fn read_all(path: &Path) -> Result<Vec<(String, String, Vec<u8>)>> {
    let mut v = Vec::new();
    for_each_record(path, |n, d, s| {
        v.push((n.to_string(), d.to_string(), std::mem::take(s)));
        Ok(())
    })?;
    Ok(v)
}

#[derive(Debug, Clone)]
pub struct BedRec {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
}

pub fn read_bed(path: &Path) -> Result<Vec<BedRec>> {
    let rd = open_maybe_gz(path)?;
    let mut v = Vec::new();
    for (i, line) in rd.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 3 {
            bail!("{}:{}: expected >=3 tab-separated BED columns", path.display(), i + 1);
        }
        let start: i64 = p[1].parse().with_context(|| format!("{}:{} bad start", path.display(), i + 1))?;
        let end: i64 = p[2].parse().with_context(|| format!("{}:{} bad end", path.display(), i + 1))?;
        if end <= start {
            bail!("{}:{}: end <= start", path.display(), i + 1);
        }
        v.push(BedRec { chrom: p[0].to_string(), start, end, name: p.get(3).unwrap_or(&"").to_string() });
    }
    Ok(v)
}

/// BED records grouped by contig, sorted and merged.
pub fn merged_by_chrom(recs: &[BedRec]) -> FxHashMap<String, Vec<(i64, i64)>> {
    let mut m: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
    for r in recs {
        m.entry(r.chrom.clone()).or_default().push((r.start, r.end));
    }
    for v in m.values_mut() {
        v.sort_unstable();
        let mut out: Vec<(i64, i64)> = Vec::with_capacity(v.len());
        for &(s, e) in v.iter() {
            match out.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => out.push((s, e)),
            }
        }
        *v = out;
    }
    m
}
