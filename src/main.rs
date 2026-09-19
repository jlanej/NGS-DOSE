//! ngs-dose: sequence-class dosage from short-read WGS. This binary is the counting engine;
//! estimation and cohort modelling live in the `ngsdose` Python package.

mod controls;
mod count;
mod fasta;
mod kmer;
mod panel;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use rust_htslib::bam::Read;
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "ngs-dose", version, about = "Sequence-class dosage from short-read WGS: counting engine")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Clone, Copy, ValueEnum, PartialEq)]
enum Mode {
    /// read every record; placement-independent, the reference mode
    Scan,
    /// retrieve only control regions and class sinks through the index
    Fetch,
}

#[derive(Subcommand)]
enum Cmd {
    /// Count control and class fragment ends in one BAM/CRAM -> per-sample counts JSON
    Count {
        /// BAM/CRAM path or URL
        #[arg(short, long)]
        input: String,
        /// index path (needed for URLs whose index should not be downloaded to the cwd)
        #[arg(long)]
        index: Option<String>,
        /// reference FASTA for CRAM decoding (otherwise REF_PATH / REF_CACHE are used)
        #[arg(short = 'T', long)]
        reference: Option<PathBuf>,
        /// k-mer panel; repeat to merge several (e.g. the bundle panel and a satellite panel)
        #[arg(short, long, required = true)]
        panel: Vec<PathBuf>,
        #[arg(short, long)]
        controls: PathBuf,
        /// BED of class sink intervals; required for --mode fetch
        #[arg(long)]
        sinks: Option<PathBuf>,
        #[arg(short, long, value_enum, default_value = "fetch")]
        mode: Mode,
        /// also retrieve the unmapped bin in fetch mode
        #[arg(long)]
        unmapped: bool,
        #[arg(short, long)]
        sample: Option<String>,
        /// output JSON (.gz for gzip); '-' for stdout
        #[arg(short, long, default_value = "-")]
        out: String,
        #[arg(short = '@', long, default_value_t = 4)]
        threads: usize,
        /// minimum panel k-mers for a read to be assigned to a class
        #[arg(long, default_value_t = 4)]
        min_hits: usize,
        /// minimum fraction of the read's k-mers that must be panel k-mers
        #[arg(long, default_value_t = 0.0)]
        min_frac: f64,
        /// bin width (bp) of positional class profiles
        #[arg(long, default_value_t = 50)]
        bin: usize,
        /// fragment-GC window lengths to tabulate (the modal read length is always added)
        #[arg(long, value_delimiter = ',', default_value = "100,150,200,250,300,350,400,450,500,550,600")]
        l_grid: Vec<usize>,
        /// bin width (bp) of the placement histogram: where class reads were aligned, and how many
        /// reads of any kind each of those bins holds. 1000 is the grid of `mosdepth --by 1000`.
        /// Compositional classes (satellite families) are recorded at ten times this width.
        #[arg(long, default_value_t = 1000)]
        place_bin: i64,
        /// fetch mode: attempts per interval (remote inputs fail transiently)
        #[arg(long, default_value_t = 5)]
        retries: usize,
        /// count a BAM/CRAM that lacks its end-of-file marker instead of refusing it
        #[arg(long)]
        allow_truncated: bool,
        /// give up with exit status 75 when nothing has been read for this many seconds (a dead
        /// HTTPS connection waits for ever rather than failing); 0 disables
        #[arg(long, default_value_t = 300)]
        stall_timeout: u64,
    },
    /// Build a k-mer panel from class FASTAs, filtered against background genomes
    Panel {
        /// TSV: name, kind (positional|compositional), fasta, circular (0/1), [min_count]
        #[arg(short, long)]
        manifest: PathBuf,
        #[arg(short, long, default_value_t = 31)]
        k: usize,
        /// background genome FASTA, optionally `fasta:mask.bed` to exempt the class's own loci
        #[arg(short, long)]
        background: Vec<String>,
        /// drop k-mers seen more than this many times in the unmasked background
        #[arg(long, default_value_t = 0)]
        max_bg: u32,
        /// write every candidate k-mer with its class and background counts (diagnostics)
        #[arg(long)]
        report: Option<PathBuf>,
        #[arg(short, long)]
        out: PathBuf,
    },
    /// Build the controls FASTA (regions + flanks) from a BED and an indexed reference
    Controls {
        #[arg(short, long)]
        bed: PathBuf,
        #[arg(short = 'T', long)]
        reference: PathBuf,
        #[arg(long, default_value_t = 1000)]
        flank: usize,
        #[arg(short, long)]
        out: PathBuf,
    },
}

fn sample_from_header(header: &rust_htslib::bam::HeaderView) -> Option<String> {
    let text = String::from_utf8_lossy(header.as_bytes()).to_string();
    for line in text.lines() {
        if line.starts_with("@RG") {
            for f in line.split('\t') {
                if let Some(sm) = f.strip_prefix("SM:") {
                    return Some(sm.to_string());
                }
            }
        }
    }
    None
}

/// The bundled libcurl has no compiled-in CA store; point it at the system one so that
/// https:// and s3:// inputs work out of the box. An explicit CURL_CA_BUNDLE always wins.
fn ensure_ca_bundle() {
    if std::env::var_os("CURL_CA_BUNDLE").is_some() {
        return;
    }
    for p in [
        "/etc/ssl/certs/ca-certificates.crt",
        "/etc/pki/tls/certs/ca-bundle.crt",
        "/etc/ssl/ca-bundle.pem",
        "/etc/ssl/cert.pem",
        "/usr/local/share/certs/ca-root-nss.crt",
    ] {
        if Path::new(p).exists() {
            std::env::set_var("CURL_CA_BUNDLE", p);
            return;
        }
    }
}

fn main() -> Result<()> {
    ensure_ca_bundle();
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Panel { manifest, k, background, max_bg, report, out } => {
            if !(11..=32).contains(&k) {
                bail!("k must be in 11..=32");
            }
            let inputs = panel::read_manifest(&manifest)?;
            let bgs: Vec<panel::Background> = background
                .iter()
                .map(|b| match b.rsplit_once(':') {
                    Some((fa, mask)) if Path::new(mask).exists() => panel::Background { fasta: fa.into(), mask: Some(mask.into()) },
                    _ => panel::Background { fasta: b.into(), mask: None },
                })
                .collect();
            let (p, st) = panel::build(k, &inputs, &bgs, max_bg, report.as_deref())?;
            p.write(&out)?;
            eprintln!("class\tkmers_input\tdropped_multicopy_or_rare\tdropped_shared\tdropped_background\tkept");
            for (n, a, b, c, d, e) in st.per_class {
                eprintln!("{}\t{}\t{}\t{}\t{}\t{}", n, a, b, c, d, e);
            }
            eprintln!("wrote {} ({} k-mers, k={})", out.display(), p.map.len(), k);
        }
        Cmd::Controls { bed, reference, flank, out } => {
            let n = controls::build(&bed, &reference, flank, &out)?;
            eprintln!("wrote {} control regions to {}", n, out.display());
        }
        Cmd::Count {
            input,
            index,
            reference,
            panel: panel_path,
            controls: controls_path,
            sinks,
            mode,
            unmapped,
            sample,
            out,
            threads,
            min_hits,
            min_frac,
            bin,
            l_grid,
            place_bin,
            retries,
            allow_truncated,
            stall_timeout,
        } => {
            let t0 = std::time::Instant::now();
            let panel = panel::Panel::load_many(&panel_path)?;
            if stall_timeout > 0 {
                count::spawn_watchdog(stall_timeout, input.clone());
            }
            let inp = count::Input { path: input.clone(), index, reference };
            // header: contig ids and the sample name
            let (header, eof_marker) = if mode == Mode::Fetch {
                let sinks_given = sinks.is_some();
                if !sinks_given {
                    bail!("--mode fetch needs --sinks (class sink intervals for this reference build)");
                }
                probe_header(&inp, true)?
            } else {
                probe_header(&inp, false)?
            };
            // a truncated file still decodes; what is lost is whatever sorted last - for a
            // coordinate-sorted GRCh38 file that is chr21, chr22 and the unplaced contigs: the rDNA
            if eof_marker == "absent" {
                if allow_truncated {
                    eprintln!("[count] WARNING: {} has no end-of-file marker (truncated?); continuing as asked", inp.path);
                } else {
                    bail!("{} has no end-of-file marker: the file is truncated (an interrupted copy?). --allow-truncated overrides", inp.path);
                }
            }
            let tid_of = |name: &str| header.tid(name.as_bytes()).map(|t| t as i32);
            let ctrl = controls::Controls::load(&controls_path, &tid_of)?;
            let params = count::Params { min_hits, min_frac, bin, l_grid, threads, place_bin, unmapped, pad: 600, retries };
            let mut sink_contigs: Vec<String> = Vec::new();
            let (acc, st) = match mode {
                Mode::Scan => count::scan(&inp, &panel, &ctrl, &params)?,
                Mode::Fetch => {
                    let recs = fasta::read_bed(sinks.as_ref().unwrap())?;
                    let mut iv = Vec::new();
                    let mut skipped = 0;
                    for r in recs {
                        if tid_of(&r.chrom).is_some() {
                            iv.push((r.chrom, r.start, r.end));
                        } else {
                            skipped += 1;
                        }
                    }
                    if skipped > 0 {
                        eprintln!("[count] note: {} sink intervals are on contigs absent from this file's header", skipped);
                    }
                    if iv.is_empty() {
                        bail!("no sink interval matches the alignment header - wrong reference build?");
                    }
                    sink_contigs = iv.iter().map(|(c, _, _)| c.clone()).collect();
                    sink_contigs.sort();
                    sink_contigs.dedup();
                    count::fetch(&inp, &panel, &ctrl, &iv, &params)?
                }
            };
            let sample = sample.or_else(|| sample_from_header(&header)).unwrap_or_else(|| "unknown".into());
            let header_contigs: Vec<(String, u64)> = (0..header.target_count())
                .map(|t| (String::from_utf8_lossy(header.tid2name(t)).to_string(), header.target_len(t).unwrap_or(0)))
                .collect();
            let o = count::make_output(
                sample,
                &inp,
                if mode == Mode::Scan { "scan" } else { "fetch" },
                &panel,
                panel_path.iter().map(|p| p.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default()).collect::<Vec<_>>().join(","),
                (
                    panel_path.iter().map(|p| sha256_file(p)).collect::<Result<Vec<_>>>()?,
                    sha256_file(&controls_path)?,
                    sinks.as_ref().map(|p| sha256_file(p)).transpose()?,
                ),
                &ctrl,
                controls_path.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default(),
                sinks.as_ref().and_then(|p| p.file_name()).map(|s| s.to_string_lossy().to_string()),
                &params,
                acc,
                st,
                &header_contigs,
                &sink_contigs,
                eof_marker,
                t0.elapsed().as_secs_f64(),
            );
            let json = serde_json::to_vec(&o)?;
            if out == "-" {
                std::io::stdout().write_all(&json)?;
            } else {
                // written under another name and renamed: a job killed mid-write leaves no file
                // that a resumed run would mistake for a finished sample
                let part = format!("{}.partial", out);
                let f = std::fs::File::create(&part).with_context(|| format!("cannot create {}", part))?;
                if out.ends_with(".gz") {
                    let mut w = flate2::write::GzEncoder::new(f, flate2::Compression::new(6));
                    w.write_all(&json)?;
                    w.finish()?.sync_all()?;
                } else {
                    let mut f = f;
                    f.write_all(&json)?;
                    f.sync_all()?;
                }
                std::fs::rename(&part, &out).with_context(|| format!("cannot move {} into place", part))?;
            }
            eprintln!(
                "[count] {}: {} primary reads, {} control 5' ends, classes: {}; {:.1}s",
                o.sample,
                o.primary,
                o.ctrl_reads,
                o.classes.iter().map(|c| format!("{}={}", c.name, c.reads)).collect::<Vec<_>>().join(" "),
                o.elapsed_sec
            );
        }
    }
    Ok(())
}

/// Content hash of a resource file, recorded in every counts file so that a cohort cannot be
/// assembled from runs against different bundles without it showing.
fn sha256_file(path: &Path) -> Result<String> {
    use sha2::{Digest, Sha256};
    let mut f = std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut h = Sha256::new();
    std::io::copy(&mut f, &mut h)?;
    Ok(h.finalize().iter().map(|b| format!("{:02x}", b)).collect())
}

/// End-of-file marker of a BAM/CRAM: "present", "absent", or "unchecked" (unseekable stream, or a
/// format without one). Must be asked before any record is read.
fn eof_status<R: rust_htslib::bam::Read>(rd: &R) -> &'static str {
    match unsafe { rust_htslib::htslib::hts_check_EOF(rd.htsfile()) } {
        1 => "present",
        0 => "absent",
        _ => "unchecked",
    }
}

fn probe_header(inp: &count::Input, indexed: bool) -> Result<(rust_htslib::bam::HeaderView, &'static str)> {
    // opening costs one small read; the engine re-opens per worker
    if indexed {
        let rd = if inp.path.contains("://") {
            let full = match &inp.index {
                Some(i) => format!("{}##idx##{}", inp.path, i),
                None => inp.path.clone(),
            };
            rust_htslib::bam::IndexedReader::from_url(&url::Url::parse(&full)?)
        } else {
            match &inp.index {
                Some(i) => rust_htslib::bam::IndexedReader::from_path_and_index(&inp.path, i),
                None => rust_htslib::bam::IndexedReader::from_path(&inp.path),
            }
        }
        .map_err(|e| anyhow::anyhow!("cannot open {} (index present?): {}", inp.path, e))?;
        Ok((rd.header().clone(), eof_status(&rd)))
    } else {
        let rd = if inp.path.contains("://") {
            rust_htslib::bam::Reader::from_url(&url::Url::parse(&inp.path)?)
        } else {
            rust_htslib::bam::Reader::from_path(&inp.path)
        }
        .map_err(|e| anyhow::anyhow!("cannot open {}: {}", inp.path, e))?;
        Ok((rd.header().clone(), eof_status(&rd)))
    }
}
