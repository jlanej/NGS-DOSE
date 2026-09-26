//! ngs-dose: sequence-class dosage from short-read WGS. This binary is the counting engine;
//! estimation and cohort modelling live in the `ngsdose` Python package.

mod controls;
mod count;
mod fasta;
mod kmer;
mod panel;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
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
        /// BAM/CRAM path or URL (http, https, ftp). A URL is written to the counts file and to
        /// messages without its query, fragment and user info, where signed URLs keep their keys.
        #[arg(short, long)]
        input: String,
        /// index path or URL. A local path is safest: without --index, and for a BAM index given as a
        /// URL, htslib saves the index in the working directory and reuses a same-named file there
        /// unchecked (a CRAM index URL is read into memory)
        #[arg(long)]
        index: Option<String>,
        /// reference FASTA the CRAM was made with (for DRAGEN CRAMs, DRAGEN's own hg38). A CRAM is
        /// refused without -T or a non-empty REF_PATH: htslib would otherwise download the reference from the
        /// EBI server (REF_CACHE alone does not stop that on a cache miss). REF_PATH set to any
        /// local directory, even an empty one, lets htslib use a populated REF_CACHE and the @SQ UR
        /// path, and is all a CRAM that needs no reference (made with no_ref or embed_ref) asks for.
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
        /// minimum panel k-mers for a read to be assigned to a class (>= 1)
        #[arg(long, default_value_t = 4, value_parser = at_least_one)]
        min_hits: usize,
        /// minimum fraction of the read's k-mers that must be panel k-mers (0-1)
        #[arg(long, default_value_t = 0.0, value_parser = fraction)]
        min_frac: f64,
        /// bin width (bp) of positional class profiles (>= 1)
        #[arg(long, default_value_t = 50, value_parser = at_least_one)]
        bin: usize,
        /// fragment-GC window lengths to tabulate (the modal read length is always added)
        #[arg(long, value_delimiter = ',', default_value = "100,150,200,250,300,350,400,450,500,550,600", value_parser = at_least_one)]
        l_grid: Vec<usize>,
        /// bin width (bp, 1 to 1,000,000,000) of the placement histogram: where class reads were aligned, and how
        /// many mapped primary reads (not secondary, supplementary or QC-fail; also how many of them
        /// are duplicate-flagged) start in each of those bins, by leftmost position. In fetch mode
        /// only reads that start inside a fetched interval are counted, so a bin that is only partly
        /// fetched is under-counted. 1000 is the grid of `mosdepth --by 1000`.
        /// Compositional classes (satellite families) are recorded at ten times this width.
        #[arg(long, default_value_t = 1000, value_parser = clap::value_parser!(i64).range(1..=count::MAX_PLACE_BIN))]
        place_bin: i64,
        /// attempts (>= 1) per fetched interval, and per open of a remote input. A remote input
        /// still failing after the last one ends with exit status 75, for the caller to retry later
        #[arg(long, default_value_t = count::DEFAULT_RETRIES, value_parser = at_least_one)]
        retries: usize,
        /// count a BAM/CRAM that lacks its end-of-file marker instead of refusing it - also a scanned
        /// URL whose marker could not be checked before reading (a server without range requests)
        /// and that turns out, read to its end, to have none. (A pipe cannot be counted: the input
        /// is opened more than once.)
        #[arg(long)]
        allow_truncated: bool,
        /// fetch mode: padding (bp) added to both sides of every control region when it is retrieved;
        /// at least 400, the longest read span it has to cover
        #[arg(long, default_value_t = count::DEFAULT_PAD, value_parser = clap::value_parser!(i64).range(count::READLEN_MAX as i64..))]
        pad: i64,
        /// fetch mode: go on when a loaded panel class has no interval in the sinks BED, or none on a
        /// contig of this file's header. Its reads are then counted only where they fall inside other
        /// intervals - an undercount - and the class is listed in the output's sinks_missing_classes.
        /// Without this flag the run refuses. (Sink intervals on contigs the header lacks are always
        /// left out and recorded per class in sinks_skipped.)
        #[arg(long)]
        allow_missing_sinks: bool,
        /// give up with exit status 75 when nothing has been read for this many seconds (a dead
        /// HTTPS connection waits for ever rather than failing); 0 disables
        #[arg(long, default_value_t = 300)]
        stall_timeout: u64,
    },
    /// Build a k-mer panel from class FASTAs, filtered against background genomes
    Panel {
        /// TSV: name, kind (positional|compositional), fasta, circular (0/1 or true), [min_count:
        /// compositional classes keep k-mers seen at least this often, default 1], [keep.bed:
        /// positional classes keep only k-mers starting in these unit intervals; '.' for none]
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
    /// Write the intervals a fetch reads (control regions padded, plus the sinks, merged) as BED.
    /// For sites that cannot let the engine read CRAMs directly: cut the reads out first with
    /// `samtools view -M -L plan.bed`, index the cut, and count it in fetch mode with the same
    /// bundle, sinks and padding - that reproduces the fetch of the whole file exactly, except for
    /// the unmapped bin, which the plan leaves out and `-L` never outputs. To reproduce `count
    /// --unmapped`, add the reads that have no coordinate to the cut (`samtools view -u in.cram '*'`,
    /// joined to the cut with `samtools cat` or `samtools merge` before indexing) and pass --unmapped
    /// when counting the cut; `plan --unmapped` prints these steps. Counted in scan mode the cut
    /// would pass for a whole-file scan, which it is not (`ngsdose sinks` refuses such a file).
    Plan {
        #[arg(short, long)]
        controls: PathBuf,
        /// BED of class sink intervals (the bundle's sinks.bed)
        #[arg(long)]
        sinks: PathBuf,
        /// a BAM/CRAM whose header decides which contigs the plan keeps (samtools ignores a region
        /// on a contig a CRAM lacks and exits 0, so an unfiltered plan can lose intervals silently).
        /// The sink intervals left out are reported per class, and, for a BAM, the reads without a
        /// coordinate that its index records
        #[arg(short, long)]
        input: Option<String>,
        #[arg(long)]
        index: Option<String>,
        #[arg(short = 'T', long)]
        reference: Option<PathBuf>,
        /// padding (bp) added to both sides of every control region, as `count` uses (>= 400)
        #[arg(long, default_value_t = count::DEFAULT_PAD, value_parser = clap::value_parser!(i64).range(count::READLEN_MAX as i64..))]
        pad: i64,
        /// print the samtools steps that add the unmapped bin to the cut, for a count with --unmapped
        #[arg(long)]
        unmapped: bool,
        /// output BED; '-' or absent for stdout
        #[arg(short, long)]
        out: Option<PathBuf>,
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

fn at_least_one(s: &str) -> std::result::Result<usize, String> {
    match s.parse::<usize>() {
        Ok(0) => Err("must be at least 1".into()),
        Ok(v) => Ok(v),
        Err(e) => Err(e.to_string()),
    }
}

fn fraction(s: &str) -> std::result::Result<f64, String> {
    match s.parse::<f64>() {
        Ok(v) if (0.0..=1.0).contains(&v) => Ok(v),
        Ok(_) => Err("must be between 0 and 1".into()),
        Err(e) => Err(e.to_string()),
    }
}

/// The bundled libcurl has no compiled-in CA store; point it at the system one so that https://
/// inputs work out of the box. An explicit CURL_CA_BUNDLE always wins. (s3:// and gs:// are not
/// read: htslib is built without its S3/GCS plugins.)
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

fn main() {
    ensure_ca_bundle();
    // htslib loads its I/O plugins at the first open, leaving errno from the search; done here, so
    // that errno after a failed open is the open's own (a 404 is final, a 503 is retried)
    unsafe { rust_htslib::htslib::hfile_has_plugin(c"libcurl".as_ptr()) };
    if let Err(e) = run(Cli::parse()) {
        // htslib's messages carry URLs, signed ones too
        eprintln!("Error: {}", count::scrub_urls(&format!("{:?}", e)));
        // a remote input that kept failing: EX_TEMPFAIL, as the stall watchdog, so the caller retries
        std::process::exit(if e.downcast_ref::<count::TempFail>().is_some() { 75 } else { 1 });
    }
}

fn run(cli: Cli) -> Result<()> {
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
        Cmd::Plan { controls: controls_path, sinks, input, index, reference, pad, unmapped, out } => {
            let inp = input.map(|path| count::Input::new(path, index, reference));
            let probe = match &inp {
                Some(i) => Some(i.probe(true, false, count::DEFAULT_RETRIES)?),
                None => None,
            };
            // without an input every contig is kept, each under an id of its own: the controls loader
            // checks for overlapping regions per id, so one id for all would fault regions on different
            // chromosomes against each other
            let ids: std::cell::RefCell<std::collections::HashMap<String, i32>> = Default::default();
            let tid_of = |name: &str| match &probe {
                Some(p) => p.header.tid(name.as_bytes()).map(|t| t as i32),
                None => {
                    let mut m = ids.borrow_mut();
                    let next = m.len() as i32;
                    Some(*m.entry(name.to_string()).or_insert(next))
                }
            };
            let ctrl = controls::Controls::load(&controls_path, &tid_of)?;
            let sk = count::load_sinks(&sinks, &tid_of)?;
            let iv = sk.intervals();
            let skipped: u64 = sk.skipped.values().map(|s| s.intervals).sum();
            let absent = ctrl.regions.iter().filter(|r| r.absent).count();
            let plan = count::fetch_plan(&ctrl, &iv, pad);
            let mut text = String::new();
            for (c, s, e) in &plan {
                text.push_str(&format!("{}\t{}\t{}\n", c, s, e));
            }
            let dest = out.as_ref().filter(|p| p.as_os_str() != "-");
            match dest {
                Some(p) => std::fs::write(p, &text).with_context(|| format!("cannot write {}", p.display()))?,
                None => std::io::stdout().write_all(text.as_bytes())?,
            }
            let bp: i64 = plan.iter().map(|(_, s, e)| e - s).sum();
            eprintln!(
                "[plan] {} intervals, {:.1} Mb: {} control regions padded by {} bp, {} sink intervals{}{}",
                plan.len(),
                bp as f64 / 1e6,
                ctrl.regions.len() - absent,
                pad,
                iv.len(),
                if absent as u64 + skipped > 0 {
                    format!(
                        "; dropped {} control regions and {} sink intervals on contigs the input lacks{}",
                        absent,
                        skipped,
                        if skipped > 0 { format!(" (sink intervals by class: {})", sk.skipped_summary()) } else { String::new() }
                    )
                } else {
                    String::new()
                },
                if inp.is_none() { " (no input given: every contig kept)" } else { "" }
            );
            let (src, bed) = (
                inp.as_ref().map(|i| i.display()).unwrap_or_else(|| "in.cram".into()),
                dest.map(|p| p.display().to_string()).unwrap_or_else(|| "plan.bed".into()),
            );
            if let (Some(i), Some(p)) = (&inp, &probe) {
                match p.unplaced_unmapped {
                    Some(n) => eprintln!(
                        "[plan] {} holds {} reads without a coordinate (from its index): a cut along the plan leaves them out{}",
                        i.display(),
                        n,
                        if n > 0 && !unmapped { "; --unmapped prints the steps that add them" } else { "" }
                    ),
                    None => eprintln!(
                        "[plan] {} is a CRAM, whose index does not record the reads without a coordinate (`samtools view -c {} '*'` counts them); a cut along the plan leaves them out",
                        i.display(),
                        i.display()
                    ),
                }
            }
            if unmapped {
                eprintln!(
                    "[plan] to count the cut with --unmapped, add the reads without a coordinate, which `samtools view -L` never outputs \
                     (for a CRAM, give samtools the reference too):\n  samtools view -b -M -L {bed} -o cut.regions.bam {src}\n  \
                     samtools view -b -o cut.unmapped.bam {src} '*'\n  samtools merge -o cut.bam cut.regions.bam cut.unmapped.bam\n  samtools index cut.bam"
                );
            }
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
            pad,
            allow_missing_sinks,
            stall_timeout,
        } => {
            let t0 = std::time::Instant::now();
            let panel = panel::Panel::load_many(&panel_path)?;
            let mut inp = count::Input::new(input, index, reference);
            if stall_timeout > 0 {
                count::spawn_watchdog(stall_timeout, inp.display());
            }
            if mode == Mode::Fetch && sinks.is_none() {
                bail!("--mode fetch needs --sinks (class sink intervals for this reference build)");
            }
            // header: contig ids and the sample name
            let probe = inp.probe(mode == Mode::Fetch, mode == Mode::Fetch, retries)?;
            inp.cram = probe.cram;
            if inp.cram && inp.reference.is_none() && std::env::var_os("REF_PATH").is_none_or(|v| v.is_empty()) {
                bail!(
                    "{} is a CRAM, and neither -T nor a non-empty REF_PATH is set: htslib would download its reference from the EBI server \
                     (slowly, into ~/.cache/hts-ref, or not at all on a node without internet; REF_CACHE alone does not stop that). \
                     Pass -T with the FASTA the CRAM was made with (for DRAGEN CRAMs, DRAGEN's hg38), or set REF_PATH (and REF_CACHE) \
                     to local copies. REF_PATH set to any local directory, even an empty one, lets htslib use a populated REF_CACHE \
                     and the path in the @SQ UR tags, and is enough for a CRAM that needs no reference (made with no_ref or embed_ref)",
                    inp.display()
                );
            }
            let header = probe.header;
            let mut eof_marker = probe.eof_marker;
            // a truncated file still decodes; what is lost is whatever sorted last - for a coordinate-sorted
            // GRCh38 analysis-set file: the unmapped reads, then the HLA, decoy, chrEBV, alt and unplaced contigs
            // (where much of the rDNA and DJ lands, e.g. chrUn_GL000220v1, chr22_KI270733v1_random), then chrM, chrY, chrX
            if eof_marker == "absent" {
                truncated(&inp, allow_truncated, false)?;
            }
            let tid_of = |name: &str| header.tid(name.as_bytes()).map(|t| t as i32);
            let ctrl = controls::Controls::load(&controls_path, &tid_of)?;
            let params = count::Params { min_hits, min_frac, bin, l_grid, threads, place_bin, unmapped, pad, retries };
            let mut sink_contigs: Vec<String> = Vec::new();
            let mut sinks_missing: Vec<String> = Vec::new();
            let mut sinks_skipped = Default::default();
            let (acc, st) = match mode {
                Mode::Scan => count::scan(&inp, &panel, &ctrl, &params)?,
                Mode::Fetch => {
                    // intervals on contigs the header lacks cannot be fetched: left out, and recorded
                    let sk = count::load_sinks(sinks.as_ref().unwrap(), &tid_of)?;
                    // every loaded class must keep sinks, or its reads are counted only where they
                    // happen to fall inside other intervals, an undercount nothing would report
                    if sk.named {
                        let have: std::collections::HashSet<&str> = sk.kept.iter().map(|r| r.name.as_str()).collect();
                        sinks_missing = panel.classes.iter().map(|c| c.name.clone()).filter(|n| !have.contains(n.as_str())).collect();
                        if !sinks_missing.is_empty() {
                            let lost: Vec<&str> = sinks_missing.iter().map(|s| s.as_str()).filter(|n| sk.skipped.contains_key(*n)).collect();
                            let msg = format!(
                                "the sinks BED has no interval {}for the loaded class(es) {}{}: a fetch would count their reads only where they fall inside other intervals",
                                if lost.is_empty() { "" } else { "on this file's contigs " },
                                sinks_missing.join(", "),
                                if lost.is_empty() {
                                    String::new()
                                } else {
                                    format!(
                                        " (every interval of {} is on a contig absent from the alignment header: sinks from another reference or pipeline?)",
                                        lost.join(", ")
                                    )
                                }
                            );
                            if allow_missing_sinks {
                                eprintln!("[count] WARNING: {}; recorded in sinks_missing_classes (--allow-missing-sinks)", msg);
                            } else {
                                bail!(
                                    "{}. Learn sinks for them from whole-file scans (`ngsdose sinks scan*.json.gz --classes ...`), load a panel without them, \
                                     or pass --allow-missing-sinks to record the gap and continue",
                                    msg
                                );
                            }
                        }
                    } else {
                        eprintln!("[count] note: the sinks BED carries no class column; whether every loaded class has sinks is not checked");
                    }
                    if !sk.skipped.is_empty() {
                        eprintln!(
                            "[count] WARNING: sink intervals on contigs absent from this file's header are left out (recorded as sinks_skipped): {}. \
                             Their classes lose whatever reads those intervals hold; sinks are specific to a reference and an aligner",
                            sk.skipped_summary()
                        );
                    }
                    let iv = sk.intervals();
                    if iv.is_empty() {
                        bail!("no sink interval matches the alignment header - wrong reference build?");
                    }
                    sink_contigs = iv.iter().map(|(c, _, _)| c.clone()).collect();
                    sink_contigs.sort();
                    sink_contigs.dedup();
                    sinks_skipped = sk.skipped;
                    count::fetch(&inp, &header, &panel, &ctrl, &iv, &params)?
                }
            };
            // a stream whose marker could not be checked up front says at its end whether it had one
            if eof_marker == "unchecked" && mode == Mode::Scan {
                match st.end_marker {
                    Some(present) => {
                        eof_marker = if present { "present" } else { "absent" };
                        if !present {
                            truncated(&inp, allow_truncated, true)?;
                        }
                    }
                    None => eprintln!(
                        "[count] WARNING: the end-of-file marker of {} could not be checked, so a truncated stream would go unnoticed{}",
                        inp.display(),
                        if inp.cram { " (a CRAM stream is checked at its end when decoded with -@ 1)" } else { "" }
                    ),
                }
            }
            if st.ctrl == 0 {
                bail!(
                    "no read of {} has its 5' end in a control region, so nothing can be estimated from it: an empty or wrong input, \
                     contig names other than the controls', an index that belongs to another file, or a fetch that returned nothing",
                    inp.display()
                );
            }
            let sample = sample.or_else(|| sample_from_header(&header)).unwrap_or_else(|| {
                eprintln!("[count] WARNING: no --sample and no @RG SM in the header: the sample is called 'unknown'");
                "unknown".into()
            });
            let header_contigs: Vec<(String, u64)> = (0..header.target_count())
                .map(|t| (String::from_utf8_lossy(header.tid2name(t)).to_string(), header.target_len(t).unwrap_or(0)))
                .collect();
            let mut o = count::make_output(
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
            o.sinks_missing_classes = sinks_missing;
            o.sinks_skipped = sinks_skipped;
            o.pipeline = count::pipeline(&header);
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

/// Refuse an input without an end-of-file marker, unless told to count it. `at_end`: the marker
/// could not be checked before reading, and the stream, read to its end, had none.
fn truncated(inp: &count::Input, allow: bool, at_end: bool) -> Result<()> {
    let how = if at_end { "ended without an end-of-file marker (it could not be checked before reading)" } else { "has no end-of-file marker" };
    if allow {
        eprintln!("[count] WARNING: {} {} (truncated?); continuing as asked", inp.display(), how);
        Ok(())
    } else {
        bail!("{} {}: the file is truncated (an interrupted copy?). --allow-truncated overrides", inp.display(), how)
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn count(extra: &[&str]) -> std::result::Result<Cli, clap::Error> {
        let base = ["ngs-dose", "count", "-i", "x.bam", "-p", "p.tsv.gz", "-c", "c.fa.gz"];
        Cli::try_parse_from(base.iter().chain(extra))
    }

    #[test]
    fn count_parameters_out_of_range_are_refused_by_the_parser() {
        assert!(count(&[]).is_ok());
        for bad in [
            &["--bin", "0"][..],
            &["--place-bin", "0"],
            &["--place-bin=-1"],
            &["--place-bin", "1000000001"],
            &["--place-bin", "1000000000000000000"],
            &["--min-frac", "2"],
            &["--min-frac", "-0.1"],
            &["--min-frac", "NaN"],
            &["--min-hits", "0"],
            &["--retries", "0"],
            &["--l-grid", "100,0"],
            &["--pad=-5"],
            &["--pad", "399"],
        ] {
            assert!(count(bad).is_err(), "{:?} accepted", bad);
        }
        for good in
            [&["--bin", "1"][..], &["--place-bin", "1"], &["--place-bin", "1000000000"], &["--min-frac", "1"], &["--pad", "400"], &["--retries", "1"]]
        {
            assert!(count(good).is_ok(), "{:?} refused", good);
        }
        let plan = |pad: &str| Cli::try_parse_from(["ngs-dose", "plan", "-c", "c.fa.gz", "--sinks", "s.bed", pad]);
        assert!(plan("--pad=100").is_err() && plan("--pad=-5000").is_err() && plan("--pad=600").is_ok());
    }

    #[test]
    fn count_defaults_are_the_documented_ones() {
        let Cmd::Count { pad, retries, bin, place_bin, min_hits, min_frac, .. } = count(&[]).unwrap().cmd else { panic!() };
        assert_eq!((pad, retries, bin, place_bin, min_hits, min_frac), (600, 5, 50, 1000, 4, 0.0));
    }
}
