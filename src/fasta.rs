//! Minimal FASTA / BED readers (plain or gzip).

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::Path;

/// Open a plain or gzip file. The file is opened once and its first two bytes are put back in
/// front of the stream, so a pipe, FIFO or `<(...)` reads correctly too.
pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    maybe_gz(f).with_context(|| format!("cannot read {}", path.display()))
}

fn maybe_gz<R: Read + 'static>(mut r: R) -> std::io::Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let mut n = 0;
    while n < 2 {
        match r.read(&mut magic[n..]) {
            Ok(0) => break,
            Ok(m) => n += m,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    let r = Cursor::new(magic[..n].to_vec()).chain(r);
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(r))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, r)))
    }
}

/// Bytes allowed on a FASTA sequence line: the IUPAC nucleotide codes in either case (1) and
/// ASCII whitespace, which is dropped (2); anything else ('-', '*', '.', digits) is an error (0).
const SEQ_BYTE: [u8; 256] = {
    let mut t = [0u8; 256];
    let codes = b"ACGTUNRYSWKMBDHVacgtunryswkmbdhv";
    let mut i = 0;
    while i < codes.len() {
        t[codes[i] as usize] = 1;
        i += 1;
    }
    let ws = b" \t\r\n\x0b\x0c";
    let mut i = 0;
    while i < ws.len() {
        t[ws[i] as usize] = 2;
        i += 1;
    }
    t
};

/// Stream FASTA records, calling `f(name, description, sequence)` for each. Whitespace inside
/// sequence lines is dropped; a byte that is not an IUPAC nucleotide code is an error.
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
    let mut ln = 0usize;
    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line).with_context(|| format!("cannot read {}", path.display()))?;
        if n == 0 {
            break;
        }
        ln += 1;
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
            if name.is_empty() {
                bail!("{}:{}: FASTA header without a name", path.display(), ln);
            }
            seq.clear();
            have = true;
        } else if have {
            let mut clean = true;
            for &b in &line {
                match SEQ_BYTE[b as usize] {
                    1 => {}
                    2 => clean = false,
                    _ => bail!(
                        "{}:{}: byte '{}' in the sequence of record {} is not an IUPAC nucleotide code",
                        path.display(),
                        ln,
                        b.escape_ascii(),
                        name
                    ),
                }
            }
            if clean {
                seq.extend_from_slice(&line);
            } else {
                seq.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
            }
        } else if !line.iter().all(|b| b.is_ascii_whitespace()) {
            bail!("{}:{}: sequence data before first FASTA header", path.display(), ln);
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
        let line = line.with_context(|| format!("cannot read {}", path.display()))?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    fn tmpfile(name: &str, content: &[u8]) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("ngsdose_fasta_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::File::create(&p).unwrap().write_all(content).unwrap();
        p
    }

    fn gz(data: &[u8]) -> Vec<u8> {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    #[test]
    fn whitespace_in_sequence_lines_is_dropped() {
        let clean = read_all(&tmpfile("clean.fa", b">u desc here\nACGTACGTAC\nGTACGT\n>v\nNNRY\n")).unwrap();
        let spaced = read_all(&tmpfile("spaced.fa", b">u desc here\r\nACGTA CGTAC \r\n\tGTACGT\t\r\n\n>v\r\nNN RY\r\n")).unwrap();
        assert_eq!(clean, spaced);
        assert_eq!(clean[0], ("u".to_string(), "desc here".to_string(), b"ACGTACGTACGTACGT".to_vec()));
        assert_eq!(clean[1].2, b"NNRY");
    }

    #[test]
    fn gap_characters_and_headerless_data_are_errors() {
        for (name, content, needle) in [
            ("gap.fa", &b">u\nACGT\nAC--GT\n"[..], "gap.fa:3"),
            ("star.fa", b">u\nACGT*\n", "not an IUPAC"),
            ("semi.fa", b">u\n;comment\nACGT\n", "not an IUPAC"),
            ("nohead.fa", b"ACGT\n>u\nACGT\n", "before first FASTA header"),
            ("noname.fa", b"> u\nACGT\n", "without a name"),
        ] {
            let e = read_all(&tmpfile(name, content)).unwrap_err();
            let msg = format!("{:#}", e);
            assert!(msg.contains(needle), "{}: {}", name, msg);
        }
        // an empty record and lowercase IUPAC codes are fine
        let v = read_all(&tmpfile("empty.fa", b">a\n>b\nacgtnmrwsykvhdb\n")).unwrap();
        assert_eq!(v.len(), 2);
        assert!(v[0].2.is_empty());
    }

    #[test]
    fn gzip_and_plain_files_read_alike() {
        let text = b">u\nACGT\nTTGA\n";
        let a = read_all(&tmpfile("p.fa", text)).unwrap();
        let b = read_all(&tmpfile("p.fa.gz", &gz(text))).unwrap();
        assert_eq!(a, b);
        // files shorter than the gzip magic
        assert!(read_bed(&tmpfile("zero.bed", b"")).unwrap().is_empty());
        assert!(read_bed(&tmpfile("one.bed", b"\n")).unwrap().is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn a_pipe_loses_no_bytes() {
        use std::os::fd::AsRawFd;
        for (data, want) in [(b"chrU\t0\t500\tU\n".to_vec(), "chrU"), (gz(b"chrV\t10\t20\n"), "chrV")] {
            let (rd, mut wr) = std::io::pipe().unwrap();
            wr.write_all(&data).unwrap();
            drop(wr);
            let path = PathBuf::from(format!("/dev/fd/{}", rd.as_raw_fd()));
            let recs = read_bed(&path).unwrap();
            assert_eq!(recs.len(), 1);
            assert_eq!(recs[0].chrom, want);
        }
    }
}
