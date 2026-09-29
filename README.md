# NGS-DOSE

Sequence-class dosage from short-read WGS: how many 45S rDNA units, 5S units, and copies of
other multi-copy classes a person carries, measured from the GRCh38 CRAMs that biobanks already
hold. The rDNA units, the distal junction and the telomeric repeat are measured by a targeted
fetch: about a minute per genome with no download, reading only the control regions and the
*sinks* where each class's reads land. Sinks are specific to the aligner and the reference, so
they are learned from whole-file scans of a small subset of each alignment pipeline (at biobank
scale, the 0.1–1% of CRAMs that can be scanned whole): two scans were enough for 45S, 5S and the
junction; the telomeric repeat's were learned from 372 (sets learned from 30 hold 99.4–99.6% of it
at the lowest in held-out genomes). The bundle ships sinks only for the NYGC bwa-mem pipeline;
DRAGEN (UK Biobank, All of Us) is not done yet. The satellite families are measured by whole-file
scan today. In the NYGC scans their reads also land in a stable set of intervals, and an
experimental sinks file learned from 100 cohort scans (`resources/experimental/sinks.satellites.bed`)
holds at least 99.8% of nine of the ten families in each of 1,648 other scans; no satellite fetch
has yet been compared with its scan. What a fetch reads is chosen from a menu, and every option has
a price in bytes of the CRAM (`ngsdose fetchplan`; see [Choosing what a fetch
reads](#choosing-what-a-fetch-reads)). A new class starts as a k-mer panel loaded in whole-file
scans; it becomes fetchable only once its sinks have been learned from those scans and checked on
others. 83 such candidate classes are ready to be scanned, and none of them has sinks yet.

**NGS-DOSE is not a CNV caller.** It measures the part of the genome that variant and CNV callers
mask. [NGS-PCA](https://github.com/jlanej/NGS-PCA) describes a cohort's coverage from the bins it
retains; NGS-DOSE measures what lives in the bins it excludes.

## How it works

```
CRAM ──ngs-dose count──▶ counts.json ──ngsdose estimate──▶ per-sample CN ──ngsdose cohort / adjust / trios──▶ cohort table
        (Rust, htslib)    70-340 kB         (Python, numpy)
```

1. **Count fragment 5′ ends**, in single-copy control regions and in every read that carries
   class-diagnostic 31-mers. A read is assigned to a class, and placed on its unit, by k-mers —
   not by where the aligner put it. `scan` reads the whole file, including the reads the aligner
   left without a position; `fetch` retrieves only the control regions and the few *sinks* where
   a class's reads are known to land (learned from scans), plus the unmapped bin with
   `--unmapped`. It returns the same counts to within what the sinks miss: with the bundle's
   sinks, over the cohort's first 1,748 scans, 0.05% of 45S reads on average, 0.02% of 5S, 0.24%
   of the distal junction and 0.14% of the telomeric repeat (at most 0.11%, 0.07%, 0.35% and
   0.61% in one genome).
2. **Model the library**: a Poisson spline of end density on fragment-scale GC, fitted per sample
   on the controls. Copy number of a window of the unit is `2 × observed / expected`.
3. **Calibrate the unit**: parts of the rDNA unit drop out far beyond what any genome-wide GC
   curve predicts (28S reads ~0.75× of 18S on NovaSeq), and *which* parts depends on the
   sequencing chemistry. Per-window efficiencies are learned by median polish over the samples
   given to `ngsdose cohort` (one fit per class). To get per-library-type efficiencies, run it
   separately per library type or chemistry, or apply a saved table with `--efficiencies`. The
   scale is pinned on *anchor* windows where different chemistries were shown to agree.
4. **Check it against known truth in every sample**: held-out autosomal sequence (2 copies),
   chrX (1 or 2), chrY (1 or 0), and the acrocentric distal junction (10 copies). All four are
   measured under the same control-fitted fragment-GC model and denominator as the positional
   classes (rDNA). The autosomal, chrX and chrY truths are counted by the aligned 5′ end in their
   regions; the distal junction goes through the same k-mer path as the classes. Mitochondrial
   genomes and EBV episomes per cell come along as covariates of the tissue or culture the DNA
   was taken from.
5. **Cohort layer**: adjustment on coverage PCs (NGS-PCA's, or the control regions' own) - as
   many as clear the noise edge of their spectrum, with a sweep against the known truths to
   confirm or overrule that number - and transmission reliability in trios to decide which
   estimator carries the most real variance.

The reasoning, and the measurements on real data behind each step, are in
[docs/DESIGN.md](docs/DESIGN.md).

## Validation on the 1000 Genomes cohort

The method was validated on the expanded 1000 Genomes cohort (3,202 genomes, 602 trios, NYGC
30× NovaSeq CRAMs): a twelve-genome pilot in which every sample also has an older library of the
same cell line on another instrument, and the cohort run that followed it, complete on 2026-09-28
(every genome scanned and fetched; 602 of the pedigree's 603 trios, the last naming a parent the
release never sequenced). All of that work — the pipeline that drives the cohort
through the method, the pilot, the counts files, the page built from them, the evidence
write-up, and the comparisons with ddPCR, with HPRC assemblies and with published estimates —
lives in its own repository, [NGS-DOSE-1000G](https://github.com/jlanej/NGS-DOSE-1000G), with the
page live at [jlanej.github.io/NGS-DOSE-1000G](https://jlanej.github.io/NGS-DOSE-1000G/). In
short, on the complete run: sequence of known copy number reads at its known copy number in every
genome (held-out autosomal sequence 1.997 ± 0.008 copies in 3,202); the targeted fetch reproduces
the whole-file scan (0.9996 of its 45S estimate, range 0.9988–0.9999); in 602 trios the 45S copy
number is inherited with a reliability of 0.95 (0.86–1.04) while the culture's and the library's
properties are not (−0.31 to 0.15); the calibrated estimate reproduces across sequencing
technologies where a read-depth ratio does not (intraclass correlation 0.98 against 0.19); and
against ddPCR it reads 0.96× the assay (r = 0.94); and the distal junction, read as whole numbers
of copies along its 400 kb, finds every partial copy that the HPRC assemblies of 28 of the genomes
resolve and is Mendelian at every position in 553 of 556 trios, while what the whole numbers leave
is kept, so that a change in part of the cells can show. The numbers, and what each finding
rules out, are in that repository's `docs/EVIDENCE.md` and `docs/DJ.md`. This repository holds the method alone, so that it can be applied
to any cohort.

## Quick start

From source (the engine needs Rust 1.88 or newer, a C compiler, libclang and CMake. htslib and
its zlib-ng, bzip2, xz and curl, and on Linux OpenSSL, are built from bundled sources; OpenSSL
also needs perl and make. The package needs numpy):

```bash
cargo build --release                     # the engine: target/release/ngs-dose
pip install -e .                          # the modelling layer: ngsdose
```

Without compiling: the container has the engine, the package, the GRCh38 bundle and `aria2c` -
a cluster needs nothing else but Apptainer, SLURM and the cohort's own scripts (the 1000 Genomes
ones are in NGS-DOSE-1000G),

```bash
apptainer pull ngs-dose.sif docker://ghcr.io/jlanej/ngs-dose:latest
apptainer exec ngs-dose.sif ngs-dose count --help
```

(`:latest` moves with every push to `main`; for a cohort run, pull one image by
`:sha-<commit>` or `@sha256:<digest>`), and a version tag (vX.Y.Z; none pushed yet, so the
container is the only no-compile route today) makes a
[release](https://github.com/jlanej/NGS-DOSE/releases) that carries the engine for Linux (glibc
2.28 and newer, curl and TLS built in) and for macOS on Apple silicon, the Python wheel, and the
resource bundle (`export NGSDOSE_RESOURCES=/path/to/resources/GRCh38`). The experimental
resources (`resources/experimental`: the satellite sinks, the sub-option BEDs, the candidate units)
are looked for beside the bundle's directory, or where `NGSDOSE_EXPERIMENTAL` points; a bare copy
of the bundle's directory alone works, but `ngsdose estimate` then says so and reports a fetch made
with a sinks BED other than the bundle's own as unverified rather than measured, since the
sub-options are unknown to it.

```bash
B=resources/GRCh38
# one genome, straight from the AWS Open Data mirror (~1 min; CRAM decoding needs the reference:
# a CRAM is refused without -T or a non-empty REF_PATH)
target/release/ngs-dose count -m fetch -@ 16 \
    -i https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram \
    -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -p $B/panel.k31.tsv.gz -c $B/controls.fa.gz --sinks $B/sinks.bed -o NA12878.json.gz
```

A remote open or fetched interval that still fails after `--retries` attempts (and a stalled
connection) ends with exit status 75, for the scheduler to retry the sample later (a lost index
request is retried like any other); a read error in the middle of a remote scan, which cannot
resume, still exits 1 at once, as does a file with no index beside it. Other errors exit 1, and
bad arguments 2. The engine's messages and the counts' `input` carry a URL with its query string
redacted; for an input or `--index` URL that has one (a signed URL), htslib's own messages, which
would print the signature in full on every failed open, are turned off before the first open, and
one line says so.

```bash
ngsdose estimate NA12878.json.gz -o estimates/ -t single_sample.tsv   # writes estimates/NA12878.estimate.json.gz
```

`estimate` checks each counts file against the bundle. A class it cannot measure (no sinks in
the fetch, sinks on contigs the input lacked, a different panel) is NA in the table, with the
reason in `<class>.status`, and a warning; a file that cannot be read or fitted is listed, the
table holds the others, and the exit status is 1.

```bash
# cohort: window calibration, coverage-PC adjustment, transmission reliability
ngsdose cohort estimates/*.estimate.json.gz -t cohort.tsv --save-efficiencies efficiencies.json \
    --segments segments.tsv             # and the whole numbers of copies called along the distal junction, one row per segment,
                                        # with the stretches that read a fraction of a copy off them
ngsdose adjust cohort.tsv --pcs ngspca/svd.pcs.txt -c rDNA45S.cn -o cohort.adjusted.tsv   # PCs above the Marchenko-Pastur edge; --n-pc N overrides
ngsdose pcsweep cohort.tsv --pcs ngspca/svd.pcs.txt -c rDNA45S.cn -p pedigree.txt -o sweep.tsv   # known-truth error and transmission for every number of PCs
ngsdose trios cohort.adjusted.tsv -p pedigree.txt -c rDNA45S.cn rDNA45S.cn_single rDNA45S.18S.flat \
    --compare-to rDNA45S.18S.flat       # reliabilities with bootstrap CIs, and paired differences
```

A cohort takes one estimate per sample (`cohort` refuses a sample given twice). `trios` and
`pcsweep` centre values within the populations of the pedigree's population column, or of a
`--population` file, and warn when the samples are not labelled with at least two populations of
five or more (a cohort whose samples all share one label is not warned about); the bootstrap resamples
whole families. A column the table lacks gets an NA row with the reason in `trios` and is left
out of a `pcsweep`; the run stops only when no column is left (`adjust` still stops, since its
output would lack what was asked for).

```bash
ngsdose selftest        # simulation checks of the statistics; needs no data
```

A cohort's validation page - known truths in every sample, fetch against scan, transmission in
trios, the coverage PCs - is built from these outputs by NGS-DOSE-1000G's `report` package; it is
written for that cohort, and shows what any cohort's page needs.

```bash
# whole-file scan, the full-accuracy mode: placement-independent, the reference for the experimental
# satellite families (no satellite fetch has been compared with its scan yet), the only mode for a
# class without sinks, and a record of where class reads were aligned and what else is in those bins
E=resources/experimental
target/release/ngs-dose count -m scan -@ 10 -i sample.cram -T ref.fa -c $B/controls.fa.gz -p $B/panel.k31.tsv.gz \
    -p $E/satellites.CHM13v2.k31.panel.tsv.gz -p $E/telomere.k31.panel.tsv.gz -o sample.scan.json.gz
```

A different aligner or decoy set places multi-copy reads differently. Before trusting `fetch`
on a new pipeline, scan a subset of whole CRAMs and check the sinks, or learn that pipeline's own:

```bash
ngsdose sinks scan*.json.gz --evaluate $B/sinks.bed      # fraction of each class the sinks capture, and the class reads left in the unmapped bin, per sample
ngsdose sinks scan*.json.gz --classes TEL -o sinks.bed   # or re-learn them (--classes TEL keeps the telomeric repeat fetchable)
```

`ngsdose sinks` takes whole-file scans only, and refuses scans aligned to references whose
contig lengths differ. Naming satellite families in `--classes` learns sinks for them as well.
In NYGC bwa-mem scans, sinks learned from 30 scans held at least 99.8% of HSat1A, HSat2, HSat3,
the α-satellite HORs, β-satellite, ACRO, SST1, CER and SATR in each of 200 other genomes (at least
99.85% in two of three random draws of the 30 and the 200; the HORs in 59.6 Mb of intervals, the
others in 0.2–3.5 Mb), and 99.4–99.6% of the telomeric repeat at the lowest, in about 0.7 Mb.
HSat1B reached only 96.6–96.8%: in the batch of 698 related samples part of it is left unmapped,
and more of it is scattered. DRAGEN re-alignments of one genome (HG00096) were checked under two
4.x versions with alt-masked references; as shares of the class's reads in its NYGC scan, the
unmapped bin, which a fetch reads only with `--unmapped`, held under DRAGEN 4.2.7 90% of 45S, 56%
of the distal junction and 64–90% of HSat1A, HSat1B, β-satellite, ACRO and the telomeric repeat,
but almost none of 5S, HSat2 or the α-satellite HORs; under DRAGEN 4.4.7 it held 72% of 45S, 53%
of the junction, 66–90% of the same five classes and 38% of 5S, with again almost none of HSat2
or the HORs (DESIGN.md §12).
Check capture with `--evaluate` on scans the sinks were not learned from. With `--stats FILE`,
`--evaluate` also writes each interval's share of its class (median, 10th percentile and largest
in any one scan) and the capture curve that `ngsdose fetchplan --capture` trims by; `--held-out`
marks the file as held out (with `-o`, the same statistics of the training scans are written and
marked in-sample). The bundle's `sinks.bed` has such statistics from 1,375 cohort scans none of
its intervals was learned from (`resources/GRCh38/sinks.stats.tsv`).
`resources/experimental/sinks.satellites.bed` holds the ten satellite families' NYGC sinks learned
this way from 100 cohort scans (25 per release batch and inferred sex), with the statistics of
1,648 held-out scans beside it (`sinks.satellites.stats.tsv`; `resources/experimental/README.md`).
Both are for NYGC bwa-mem alignments to the GRCh38 analysis set only.

An engine since 645ae55 refuses a fetch with a loaded panel class that has no interval in the
sinks, or whose intervals all lie on contigs the input's header lacks, and one from 0.1.1 on also
a class that loses any interval to such a contig: its reads would be counted only where they fall
inside other intervals, and `ngsdose estimate` would not measure the class. `--allow-missing-sinks` (also since 645ae55) instead
records the gap in the counts (`sinks_missing_classes`), and `--classes A,B` (engines from the
fetch-menu change of 2026-09-26 on) counts only the named classes, so only their sinks are read and checked. fae1124, the
cohort's engine, has none of the three: it does not refuse, it undercounts such a class without
saying so, which is why `fetchplan` refuses to plan such a fetch for it (below). Sink intervals on contigs the header
lacks cannot be read: a class that loses any of them is refused like a class without sinks,
unless `--allow-missing-sinks`, which records them per class (`sinks_skipped`) and fetches the
rest; `ngsdose estimate` reports such a class as NA rather than as an undercount, so the flag
spends a fetch on a class that will not be measured (engines from 0.1.1 on; the build of
2026-09-26 warned and fetched).

Where the engine cannot read the CRAMs itself, `ngs-dose plan` writes the intervals a fetch reads
as BED, to cut them out with samtools and count the cut in fetch mode with the same bundle, sinks
and padding. That reproduces the fetch of the whole file exactly, except for the unmapped bin,
which the plan leaves out and `samtools view -L` never outputs. To reproduce a
`count --unmapped`, `ngs-dose plan --unmapped` prints the samtools steps that add the reads
without a coordinate to the cut, and `plan -i` reports how many there are when the input is a
BAM (a CRAM index does not record them). A cut counted in scan mode would pass for a whole-file
scan, which it is not; `ngsdose sinks` refuses such a file, when learning and when evaluating.

```bash
target/release/ngs-dose plan -c $B/controls.fa.gz --sinks $B/sinks.bed -i sample.cram -T ref.fa -o plan.bed
samtools view -T ref.fa -b -M -L plan.bed -o sample.cut.bam sample.cram && samtools index -c sample.cut.bam
target/release/ngs-dose count -m fetch -i sample.cut.bam -p $B/panel.k31.tsv.gz -c $B/controls.fa.gz --sinks $B/sinks.bed -o sample.json.gz
```

Each counts file also records the input's pipeline from its header (`pipeline`: the @PG
programs, and a hash of the @SQ names, lengths and M5 checksums), so that counts from different
alignments can be told apart. A fetch's record is compared with the pipeline the bundle's sinks
were learned under (`bundle.json`, `sinks_learned_from.pipeline`: the aligner named in @PG and the
hash of the analysis set's @SQ lines): `ngsdose estimate` warns when either differs and records
the verdict in `sinks_pipeline`, since sinks learned under one aligner and reference do not hold
another's reads (DESIGN.md §12). Counts from engines before the record, and headers stripped of
their @PG lines or of contigs, are not compared.

The 1000 Genomes cohort run - its pipeline for SLURM or a plain loop, its pilot, its counts files
and the page built from them - is [NGS-DOSE-1000G](https://github.com/jlanej/NGS-DOSE-1000G),
published as the run proceeds. `resources/build/` rebuilds the GRCh38 bundle from public inputs,
pinned by sha256, into `work/bundle/GRCh38`, and records in `build_manifest.tsv` whether each
output is identical to the shipped file.

## Choosing what a fetch reads

A fetch decodes whole CRAM slices: every slice whose reads overlap a control region or a sink
interval is read in full. What it costs is therefore set by where the intervals fall in the file,
not by their length: the 0.5 Mb of 45S sinks cost almost as much as the 13 Mb of control regions,
because they share their slices with piles of other reads. `ngsdose fetchplan` prices each option from a CRAM's own index
(`.crai`) as `ngs-dose count -m fetch` reads it: one indexed fetch per run of touching or
overlapping intervals, each decoding every slice that overlaps it with its container's
compression header, so a slice under several runs is decoded once per run. Beside that figure the
plan reports the floor with every slice decoded once (`cum_mb_floor`), which a reader that sorted
the plan's slices would reach: a few percent below in men and 11-13% in women on `core_tel`,
whose few sparse chrY slices are decoded once per chrY truth region. It writes the files that
`ngs-dose count -m fetch` takes. The options are
the rows of `resources/fetch_menu.tsv`; its header documents the format, and its tiers and presets
are meant to be edited. The menu is found beside the bundle's directory (the repository, the image
and the release tarball keep it there) or in a source checkout; with `NGSDOSE_RESOURCES` set and no
menu beside it, `fetchplan` stops rather than take the checkout's, whose rows would name another
bundle's files (`--menu` names one outright).

The options, and the capture of their sinks in held-out cohort scans (min / median):

| option | status | tier | capture, min / median |
| --- | --- | --- | --- |
| controls: 800 control regions, 180 truth regions, chrM, chrEBV (always read) | shipped | A | |
| `rDNA5S` | shipped | A | 99.93 / 99.98% |
| `DJ` | shipped | A | 99.65 / 99.76% |
| `rDNA45S` | shipped | A | 99.89 / 99.95% |
| `TEL` | shipped (no fetch compared yet) | B | 99.39 / 99.86% |
| unmapped bin (`--unmapped`) | shipped | C | |
| `SST1` | experimental | C | 99.90 / 99.97% |
| `CER` | experimental | D | 99.95 / 99.98% |
| `SATR` | experimental | D | 99.83 / 99.93% |
| `ACRO` | experimental | D | 99.92 / 99.96% |
| `HSat1A` | experimental | D | 99.98 / 99.99% |
| `HSat2` | experimental | C | 99.99 / 100.00% |
| `bSat` | experimental | D | 99.93 / 99.96% |
| `HSat1B` | experimental | D | 96.59 / 99.02% (97.92 / 99.42% with the unmapped bin) |
| `aSatHOR` | experimental | C | 99.85 / 99.89% |
| `HSat3` | experimental | D | 99.91 / 99.97% |
| `DYZ3`: the chrY centromeric HOR array, 4 of `aSatHOR`'s intervals | experimental | C | sub-option |
| `DXZ1`: the chrX centromeric HOR array, 4 of `aSatHOR`'s intervals | experimental | C | sub-option |
| `DYZ2`: an index of male-specific `HSat1B` (78 male-only intervals, about a quarter of it) | experimental | C | sub-option |
| `DYZ1`: Yq12 `HSat3`, 1 of its intervals | experimental | C | sub-option |
| 83 candidate classes | candidate | | no sinks: scan only |

What each option costs, alone and in plans, is in [docs/fetch_examples.md](docs/fetch_examples.md):
`resources/build/fetch_examples.sh` runs `fetchplan` for a set of worked examples on the indexes
of 13 NYGC 1000 Genomes CRAMs (median 16.5 GB) and writes each plan's table, count flags and total.
The costs, like the sinks, are those of NYGC bwa-mem alignments to the GRCh38 analysis set. On
them, the controls and all of the bundle's sinks (`--preset core_tel`) read 546.9 MB per genome,
3.33% of the CRAM, and every option with sinks 3,416.7 MB (21.59%), medians over the 13 (522.7 and
3,099.5 MB with every slice decoded once).

`shipped` means that the sinks are in the bundle. For rDNA45S, rDNA5S and DJ a fetch returns
what the scan placed in them (fetch / scan reads inside the sinks 1.00000–1.00019 over the cohort's
1,748 genomes); no fetch with the TEL sinks has been compared with its scan, so TEL, like the
satellites, rests on scan placements. The captures of the bundle's four classes are over 1,375
cohort scans that none of their sinks was learned from (the 1,748 counted by 2026-09-25, less the
372 the TEL sinks were learned from and HG02258); those of the satellites over the 1,648 cohort
scans not among the 100 their sinks were learned from. A placement bin counts as captured only if
all of it lies inside a sink, and a fully unmapped read counts as missed. Both sets of statistics
are held out (their headers say `# held-out: yes`), carry each interval's largest share in any one
scan (`share_max`), ship, and are named by the menu (`resources/GRCh38/sinks.stats.tsv`,
`resources/experimental/sinks.satellites.stats.tsv`), so every class of the menu has statistics a
capture target can trim by, per byte of its CRAM slices. Two kinds of option are not trimmed: the
four sub-options, which are fetched whole, and `HSat1B` at targets above 0.9717 at the default
10th-percentile statistic (its capture with all its sinks there), where all its intervals are kept;
with `--capture-stat median`, whose full capture for it is 0.990, a target up to 0.990 trims it. The sub-options
are named parts of a satellite family's sinks: they are counted as their family and reported by
`ngsdose estimate` as `DXZ1.mass_Mb` and so on (`resources/experimental/README.md` says what each
measures and how well).

```bash
ngsdose fetchplan --list                           # options, tiers and presets
ngsdose fetchplan --preset biobank_lite --capture 0.995 --engine target/release/ngs-dose \
    --crai HG00096.crai HG00097.crai ... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai -o plan
target/release/ngs-dose count -m fetch -i sample.cram -T ref.fa -c "$(cat plan.controls.txt)" --sinks plan.sinks.bed \
    $(sed 's/^/-p /' plan.panels.txt) $(cat plan.count_flags.txt) -o sample.json.gz
ngsdose estimate sample.json.gz --fetch-sinks plan.sinks.bed -o estimates/
```

- **Selecting.** `--classes` and `--preset` name options. `--budget-mb` adds options in tier order
  (A to D), cheapest first within a tier, and stops at the first that does not fit, so that a
  lower tier never displaces a higher one; `--fill` skips it and goes on. The controls count
  toward the budget. Candidates are never fetched: they have no sinks yet, and a plan only lists
  their panels in `PREFIX.scan_panels.txt` for the whole-file scans.
- **Capture targets.** `--capture 0.995` (or `--capture-class TEL=0.99`) keeps, per class, the
  intervals of highest yield until the capture reaches the target at the 10th percentile of the
  statistics' scans (`--capture-stat median` for the median), and drops the rest. Without an
  index, yield is a class's share per read of any class in the interval. With `--crai` it is the
  share per byte of the interval's own CRAM slices, so an interval of low share on costly slices
  (a decoy slice shared with other contigs, a pile-up bin) is among the first dropped; the expected
  capture is then a lower bound (the capture of all the class's intervals less the largest share
  each dropped interval held in any one scan, or the per-read curve up to the first interval
  dropped, whichever is larger), and where the per-read order reaches the target with fewer bytes,
  that set is kept instead (`order` in the plan). Checked against each held-out scan's own capture
  at targets of 0.95–0.999, the bound was never above the truth for the bundle's classes. Classes
  share intervals: one that another selected class keeps is read anyway, so every class that has
  it keeps it, and with `--crai` so is any interval whose CRAM slices the plan reads for another
  option (each class is then trimmed again with those slices free). The pile-up bin
  chr2:32,909,000–32,921,000 is among `rDNA45S`'s first drops, but it holds 1.3% of `TEL`'s reads
  in the median held-out scan (up to 3.4%), so `TEL` keeps it at targets above about 0.979, and
  with `TEL` in the plan `rDNA45S` keeps it too, at no extra cost. Per option the plan reports
  `mb_saved`, what its trimming saved of the plan with the rest of the plan as it is (slices
  another option reads are no saving, so these need not add up to the total saving, which a note
  gives), and `capture_lost`. The statistics come from `ngsdose sinks HELD_OUT --evaluate BED --held-out
  --stats FILE`; `--held-out` records that the scans were not used to learn the sinks, and the
  plan says for each class whether its expected capture is held out. Examples 4, 6, 11 and 13 of
  [docs/fetch_examples.md](docs/fetch_examples.md) are trimmed plans.
- **What it writes.** `PREFIX.sinks.bed`, `PREFIX.panels.txt` (one `-p` each),
  `PREFIX.controls.txt` (the `-c` FASTA whose regions were costed: `controls.fa.gz`, or
  `controls.lite200.fa.gz` with `--controls controls.lite200.bed`) and `PREFIX.count_flags.txt`
  for the fetch; `PREFIX.scan_panels.txt`, every selected panel, candidates included, for the
  whole-file scans; and `PREFIX.plan.tsv`, the cost table, headed by every option that shaped
  the plan. The flags hold `--unmapped` when the unmapped bin is selected. When a loaded panel also
  defines classes that were not selected (fetching `HSat2` alone loads the satellite panel),
  `fetchplan` runs `ENGINE count --help` on the engine given by `--engine` (default: `ngs-dose` on
  `PATH`) and writes what it takes. An engine with `--classes` (from 2026-09-26 on) gets
  `--classes=` and the classes to count: the others are neither counted nor listed. One with only
  `--allow-missing-sinks` (645ae55 up to that change, e.g. 7772e32) gets that flag: the counts file lists
  the others in `sinks_missing_classes` and `estimate` leaves them NA (`no_sinks_in_fetch`). The
  cohort's engine, fae1124, takes neither, and `fetchplan` refuses the plan; `--unmarked-companions`
  writes it anyway, and the fetch then counts those classes only inside the plan's intervals
  without saying so, so `estimate` must get `--fetch-sinks PREFIX.sinks.bed` to mark them
  `no_sinks_in_fetch`. A plan of whole panel files, such as `core`, `core_tel` or `satellites`,
  loads no unselected class and needs neither flag, so fae1124 can run it. Without an engine to
  probe, `fetchplan` writes `--allow-missing-sinks`, which fae1124 rejects as an unknown option, so that fetch fails at once rather than undercounting.
- **Counting some classes only.** `ngs-dose count --classes A,B` loads and merges every panel as
  without it, so the k-mers shared between panels are dropped as in the full load, and counts each
  named class exactly as the full load does in scan mode; in a fetch, the named class's reads that
  lie inside other classes' sink intervals are not read (1 of 54,767 rDNA45S reads in the test
  fixture), which `classes_selected` records; the counts file records the selection
  (`classes_selected`) and lists only those classes. In a fetch it reads only the named classes'
  sink intervals (and rows without a class) when the BED names its classes, so the bundle's own
  `sinks.bed` can be fetched for fewer classes and fewer bytes; `ngs-dose plan --classes` writes
  the same intervals. It does not lower peak memory, since every panel is still loaded.
- **Costs and sinks belong to one pipeline.** These were measured on NYGC bwa-mem CRAMs of the
  GRCh38 analysis set, and the experimental sinks were learned from scans of that pipeline.
  Another pipeline costs its options with its own indexes (`--crai`, the median over several; each
  from the same place as its CRAM, since an index of another file of the same sample gives wrong
  numbers) and fetches with sinks learned from its own scans (`--sinks`).
- **The controls are the floor.** They cost more than any shipped class. The bundle's
  `controls.lite200.bed` / `controls.lite200.fa.gz` keep 200 of the 800 control regions and all 182
  truth and dosage regions: 120.9 MB of CRAM slices instead of 260.0 (medians over the 13). Over
  the cohort's first 1,748 genomes, 45S copy number with them differs from that with all 800 by a median +0.21% (SD 0.30%, largest 1.1%), against a 3.8–3.9%
  SD between independent libraries of the same cell line; 45S transmission reliability is unchanged
  (0.9993 and 0.9992 in 385 trios). These figures come from GC tables rebuilt from each scan's
  per-region counts, not yet from lite fetches of whole CRAMs. The price: the aneuploidy test
  flags 20 of the 32 samples the full set flags, the all-window 45S estimate rises by 1.5%, and
  HSat1B moves with an SD of 2.3% (up to 23% in one genome). A cohort should not mix the two
  (`ngsdose estimate` warns when its inputs do). `ngsdose estimate` accepts counts made with the
  lite set, as it accepts any subset the bundle names (`controls.<name>.bed`), and refuses any
  other set of control regions; the table records `controls_used` and `controls_subset`. Plan with
  `--controls controls.lite200.bed`, and `PREFIX.controls.txt` names the FASTA to fetch with.

### New classes: scan first

Sinks are learned, never assumed from a class's reference coordinates, because the aligner puts
reads where it can place them, and that depends on the pipeline and on the person. A new class
goes through four steps before any fetch reads it:

1. **Load its panel in whole-file scans.** `resources/experimental/candidates/` holds 83 candidate
   classes in seven panel files (5,076,728 k-mers; their sha256 are listed in its README):
   tandem macrosatellites, multi-copy and deleted genes, chrY and chrX arrays, RNA-gene arrays,
   viruses and Mycoplasma, and long coding VNTRs, each with the truth available for 1000 Genomes
   samples (`candidates.tsv`). No k-mer is shared with a shipped panel or between two of the
   files, and loading all seven changes no count of a shipped class (the test fixture with three
   engine builds, and a whole 30× NA12878 scan). With the final files, fae1124 scanned that CRAM
   in 119.5 s and 1.65 GB of peak memory, against 106.1 s and 1.43 GB without them (one run on a
   shared machine; repeated runs over the review rounds gave 109.6–119.5 s and 1.49–1.81 GB with
   them), most of the memory for the Mycoplasma panel, which is a file of its own. A plan must
   fetch at least one class with sinks, so pair the candidates with a fetchable preset:
   `ngsdose fetchplan --preset core candidates -o PREFIX` (or `core_tel`, and `candidates_A` or a
   group in place of `candidates`) writes every selected panel, candidates included, to
   `PREFIX.scan_panels.txt`; `--preset candidates` alone is refused and writes nothing.
2. **Learn the sinks** from a training set of those scans:
   `ngsdose sinks TRAIN/*.json.gz --classes CLASS ... -o new.bed` (positional classes are learned
   whether named or not).
3. **Check them on held-out scans**, and write the statistics a capture target needs:
   `ngsdose sinks HELD_OUT/*.json.gz --evaluate new.bed --held-out --stats new.stats.tsv`.
4. **Point its menu row at that BED and those statistics** (status `experimental`; every
   candidate already has a row, of status `candidate`). `fetchplan` then prices it and can fetch it.

One whole-file scan of NA12878 with every candidate loaded shows why step 1 comes first: 58 of the
72 classes with at least 25 reads had at least 98% of them at their GRCh38 reference copies, and
14 did not. `GSTT1` had 65% (GRCh38 carries it only on an alt contig), `D4Z4_4qA` 81%, `CCL3L` 85%
and eight others 90.5–97.6%: reads of those classes were placed away from the reference copies. The
other three had few or none of their reads there, and these are off-target floors rather than
placements: `UGT2B17` (0%; NA12878 carries its deletion, and its 174 reads are off-target hits),
`DAZ` (0.9% of 109 reads) and `RBMY` (0% of 47), a deleted gene and two chrY classes in a female.
The first build also gave `HHV7`, which has no GRCh38 copy, 283 reads at chromosome ends from
simple-repeat k-mers; the periodic-k-mer filter now applied to every candidate panel removed them,
and the final files give it none. `resources/experimental/candidates/reference_copies.GRCh38.bed`
lets `ngsdose sinks --evaluate` measure, on held-out scans, what fetching the reference copies
would lose. From the placement bins of that one female genome, fetching all 83 would add about
109 MB (0.66%) to `core_tel` on the 13 NYGC CRAMs with every slice decoded once (the floor; what the
engine reads is more, as it decodes a slice once per run of intervals): a planning figure, not sinks. The chrY
candidates (`TSPY`, `RBMY`, `DAZ`, `BPY2`, `CDY1`, `CDY2`, `DYZ19`) have no bins in a female apart
from off-target reads (`DAZ` 0.27 MB), so their real sinks are not priced. NGS-DOSE-1000G's
pipeline carries the same route (`CANDIDATE_PANELS` for the scans, `pipeline/06_learn_sinks.sh`
for steps 2–4, `FETCH_PRESET`, `FETCH_CLASSES`, `FETCH_BUDGET_MB` and `FETCH_CAPTURE` for the
fetch).

A candidate's copy number comes out of `ngsdose estimate` as for the satellites when the class is
compositional (the VNTRs, `DYZ19`, `MYCO`, `EBV2`: a sequence mass). A positional candidate needs
its unit sequence (`resources/experimental/candidates/units/`, or a directory in
`NGSDOSE_EXTRA_UNITS`), and its estimate is marked `experimental`: every usable window, one sample
at a time, no anchor windows and no cohort calibration. Such an estimate can carry a level offset;
for comparison, the all-window 45S estimate is 0.935–0.978 of the anchored one (median 0.958) in 60
cohort genomes, and DJ's 0.996–1.025. Without a unit the class is reported as skipped, and the
reason says where a unit was looked for.

## Output columns (per sample)

| column | meaning |
| --- | --- |
| `rDNA45S.cn` | cohort-calibrated diploid copy number (`ngsdose cohort`); NA for a sample whose estimate lacks the class or has no usable window |
| `rDNA45S.cn_single` | single-sample headline: anchor windows under the fragment-GC model, or every usable window when the anchors hold fewer than 1,000 fragment ends (`.cn_basis` says which, `.n_anchor` how many ends) |
| `rDNA45S.18S`, `.28S`, … / `.flat` | per-feature estimates with / without the GC model; `18S.flat` is the estimator used in the UK Biobank literature |
| `rDNA5S.cn` | 5S units |
| `DJ.cn`, `DJ.cn_unit` | distal junction (expected 10): the level on the core of the unit (the intervals where junction copies differ left out), on a scale pinned to the cohort's mode when the cohort has fifty genomes or more; and the level over the whole unit on the same scale, which is what `DJ.cn` was before the class had rules (`calibration.json`; `--no-class-rules` restores it) |
| `DJ.copies`, `DJ.partial`, `DJ.variants` | whole numbers of copies called along the unit (`segments.py`): the copies the genome is described against (ten where it holds ten over 40 kb or more of the core); the copies that hold (+) or lack (−) an end of the unit, as `+1:0-316kb`; every event, the local ones and the polymorphic intervals' included. `none` where there is none |
| `DJ.call`, `DJ.call_gap`, `DJ.scale_f`, `DJ.tilt` | `settled` (the genome sits on its whole numbers), `fractional` (the nearest whole numbers are clear and the genome is not on them: its level, or a stretch of its unit) or `uncertain` (the level lies between two whole numbers throughout, so that another whole number of copies explains the profile nearly as well: `call_gap` is how many log units that reading is behind; below 3 is uncertain; empty where no other reading is within reach); the genome's scale, and its lean (the log change of its profile across the unit), both fitted with the chain |
| `DJ.off`, `DJ.off_z`, `DJ.fractional`, `DJ.fractional_z` | what the whole numbers leave, for a change in part of the cells to show in: the level less the whole numbers called, in copies, and the same in robust SDs of the cohort's scales (3 or more makes the call fractional); the stretches of the unit that read a fraction of a copy off their whole number, as `+0.53:106-400kb` (height against the copies the genome is described against), and each stretch's z against the same contrast in the cohort's other genomes (4 or more, and a quarter of a copy, for a stretch to be kept). With fewer than fifty genomes no stretch is looked for, and the level is judged against the spread a saved efficiency table or the bundle records. What a fraction is, the measurement cannot say: a change in part of the cells, or a library unlike the cohort's |
| `HSat3.mass_Mb`, `ACRO.mass_Mb`, `TEL.mass_Mb`, … | the experimental panels: diploid sequence mass of a satellite family (scan mode, or a fetch through satellite sinks learned for the pipeline: `resources/experimental/sinks.satellites.bed` for NYGC bwa-mem, not yet compared with scans) or of the telomeric repeat (either mode; its sinks are in the bundle). The same column for a compositional candidate class (`VNTR_ACAN.mass_Mb`, `MYCO.mass_Mb`, …) when its panel was loaded. How far each is validated: `resources/experimental/README.md` |
| `rDNA45S.status`, … | `ok`, or why the class is NA: `no_sinks_in_fetch` (a class the fetch's sinks gave no interval: listed in `sinks_missing_classes` by an engine with `--allow-missing-sinks`, or, for a fae1124 fetch of a `--unmarked-companions` plan, found from the plan's BED given to `--fetch-sinks`), `sinks_skipped`, `panel_mismatch`, `subset_only` (a fetch that read the family only at its sub-options), `unverified` (`aSatHOR`, `HSat1B` or `HSat3`, the families with sub-options, counted by a fetch whose sinks BED `estimate` does not know, such as a fetchplan `PREFIX.sinks.bed`: it cannot tell whether the fetch read the whole family or only its sub-options, so the family is not measured until `--fetch-sinks PREFIX.sinks.bed` names the BED), or `skipped: <reason>` for a positional class that neither the bundle nor an experimental unit covers (the reason says where a unit was looked for). `experimental` for a positional class estimated from an experimental unit: all usable windows, no anchor, no cohort calibration; its `.cn_se_rel` is the single-sample relative error |
| `DXZ1.mass_Mb`, `.reads`, `.status` (and `DYZ3`, `DYZ1`, `DYZ2`) | the sub-options of `resources/experimental/subsets/`: the family's reads placed inside the named intervals, and their mass (reads × the family's mass per read). Status `ok`; `unverified` for a fetch whose sinks BED `estimate` does not know (reads given, no mass, and the family itself `unverified` too; pass `--fetch-sinks PREFIX.sinks.bed`); `not_fetched`, `not_counted`, `class_not_measured`, `not_compositional`, `name_clash` or `no_placements` (counts from an engine that wrote no placements) when there is nothing to measure |
| `controls_used`, `controls_subset` | how many control regions the counts held (800 with the bundle's controls, 200 with `controls.lite200`), and the name of the bundle's lighter set when they were made with one |
| `*.adj` | an estimate with coverage PCs regressed out (`ngsdose adjust`); NA for a sample without those PCs, or for a column with too few usable values |
| `truth.auto`, `truth.chrX`, `truth.chrY` | held-out known-copy-number sequence (expected 2; 1 or 2; 1 or 0) |
| `chrM.copies`, `chrEBV.copies` | mitochondrial genomes and EBV episomes per cell: covariates of the state of the tissue or cell line, measured like the truths |
| `engine` | engine version and the commit it was built from, as recorded in the counts file |
| `eof_marker` | `present`; `absent` if the input lacked its end-of-file block (the engine refuses such files unless `--allow-truncated`); `unchecked` where it could not be checked (a CRAM stream decoded with several threads, or a remote check that kept failing); NA for counts files without the field. `estimate` warns on `absent` and `unchecked` |
| `unmapped_fetched` | fetch mode: whether the unmapped bin was read (`--unmapped`) |
| `pipeline_sq_sha256` | the hash of the input's @SQ lines, for counts made by engines that record the pipeline |
| `depth`, `ctrl_dup_frac`, `gc_rel_35`, `gc_rel_65`, `ctrl_region_sd`, `flagged_chromosomes` | library and sample QC; an aneuploid chromosome, or a large arm-level gain, is reported and excluded from the denominator |
| `untestable_chromosomes` | chromosomes with too few control regions (fewer than 3 after outlier trimming) to be tested for aneuploidy: chr22 in the GRCh38 bundle |
| `*.profile_sd`, `*.profilePC*` | how far the sample's window profile departs from the cohort's |
| `ctrlPC1…`, `ctrlPC_mp` | components of the control regions' residual depth across the cohort: internal technical covariates, usable by `ngsdose adjust` when NGS-PCA has not been run; `ctrlPC_mp` is how many of them stand above the noise edge (what `adjust` uses by default). None below 10 samples with control residuals; NA for a sample without them, or whose control regions are not the same set in the same order as the others' (the `region_order_sha256` an estimate records: another controls file or bundle revision), which are named |
| `gc_curve_max_se` | how well the sample's GC curve is determined (large at very low depth) |

For a class estimated from an experimental unit, the per-sample estimate files
(`*.estimate.json.gz`) also record `estimator`, `cn_basis` (`all`), `cn_se_rel` (the window spread
over the square root of the number of windows, plus the Poisson term), `unit_source` (the unit's
path, relative to the install root when under it) and `unit_sha256`.

## Status

Engine, estimator, cohort layer (with the distal junction's rules: core, pinned scale, polymorphic
intervals, whole numbers of copies along the unit) and the GRCh38 bundle (45S, 5S, DJ) are implemented and tested:
Rust unit tests, a simulated genome with known truth run end to end, a 2% subsample of real
NA12878 reads, a mock trio cohort (that subsample sixty times over) through the whole cohort
layer and bundle-integrity checks, all in CI (Linux with Python 3.10, 3.12 and 3.13; macOS with
Python 3.12). Every push to `main` publishes the container image once its smoke test on the
NA12878 subsample passes, and a version tag makes a release with prebuilt engines, the Python
package and the resource bundle (`.github/workflows/`). Every assumption the counts files rest
on was audited against data before the cohort run, and reviews extended the audit during it;
DESIGN.md §15 lists what held and what was wrong, each with its fix. Validated on
the 1000 Genomes cohort in NGS-DOSE-1000G (a twelve-genome pilot with independent library
replicates, then the cohort run: trios, two counting modes, ddPCR, assemblies, published
estimates); its results page covers all 3,202 genomes and 602 complete trios (the run finished on
2026-09-28).
Experimental panels - ten satellite families, measured by scan, with experimental NYGC sinks
that no fetch has yet been compared against, and the telomeric repeat, which the bundle's sinks
make fetchable (no TEL fetch has yet been compared with its scan) - ship under
`resources/experimental/`, and 83 candidate classes, which have no sinks yet, under
`resources/experimental/candidates/`. Against HPRC release-2 assemblies of 200 cohort members (the results
page), the median estimate/assembly ratio is 0.94 for HSat3 and 1.02 for the α-satellite HORs,
0.91 and 0.86 for HSat1A and HSat1B, and 0.43–0.79 for CER, β-satellite and ACRO, which read low
by their k-mer recall. Across people HSat1B tracks the assemblies at r = 1.00, ACRO at 0.95,
β-satellite at 0.93 and CER at 0.91, HSat1A and HSat3 at 0.79 and 0.76; the HORs reach 0.70
because people differ by only 5% while the two measurements agree to 3% per genome. HSat2 does not
track the assemblies with at most 2% of its arrays in marked gaps (r = 0.39 in 79) although it is
inherited. SST1 and SATR are annotated so differently in HPRC and in
CHM13 that their absolute ratios mean nothing. `resources/experimental/README.md` says how far
each is validated, and NGS-DOSE-1000G holds the current numbers.
Not yet done: sinks for DRAGEN-aligned data (UK Biobank, All of
Us), which must be learned from whole-file scans of a subset of those CRAMs; real satellite
and TEL fetches compared with their scans, without which the satellite families are measured
only where a whole file is scanned - at biobank scale, a subset; sinks for the candidate classes, from the
scans still to run; a wider orthogonal rDNA calibration than the twelve ddPCR lines. See
DESIGN.md §12–13. No licence has been chosen yet.

## Provenance and credit

Method design, statistics and implementation were developed by **Claude (Anthropic)** in
September 2026 working sessions with **@jlanej** — first as an offshoot of a review of the
acrocentric short arms (Claude Opus 5), then as the ground-up rewrite in this repository
(Claude Fable 5.1), in which the original specification was tested against real 1000 Genomes
data and largely replaced. The naming, the scoping and the decision to build it are shared; the
errors are worth attributing to the machine that made them until a human has checked each one.

It builds on prior work that should be cited ahead of this repository: NGS-PCA for the exclusion
set behind the controls and for the coverage components; htslib; the T2T-CHM13 assembly and
CenSat annotation; the 1000 Genomes 30× resource (Byrska-Bishop et al., *Cell* 2022); Benjamini &
Speed (*NAR* 2012) for the fragment-GC model; Hall, Turner & Queitsch (*Sci Rep* 2021), whose
per-sample table for the same CRAMs is the external check used here; and Rodriguez-Algarra,
Evans & Rakyan (*Cell Genomics* 2024) for the demonstration that rDNA copy number carries
phenotypic signal. Unit sequences are GenBank KY962518.1 and X12811.1.

Two honesty notes. The transmission-reliability framing is standard midparent regression for an
inherited trait; no claim of novelty is made for it. The claims of what is specific to NGS-DOSE in
DESIGN.md §14 were checked against a literature search on 2026-09-26, which narrowed several of
them (its references were verified against PubMed or Europe PMC); it was a broad search, not a
systematic review. And an
AI system is not an author under prevailing journal policy — if this becomes a paper, the
appropriate form is a contributions or acknowledgements statement describing what was
machine-generated, with human authors taking responsibility for verification.
