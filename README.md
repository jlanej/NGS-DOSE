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
scan today. In the NYGC scans their reads also land in a stable set of intervals that learned
sinks capture (see below), but no satellite sinks ship, so where none have been learned they are
measured only in the scanned subset.

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

The method is being validated on the expanded 1000 Genomes cohort (3,202 genomes, 602 trios, NYGC
30× NovaSeq CRAMs): a twelve-genome pilot in which every sample also has an older library of the
same cell line on another instrument, and the cohort run that follows it, which is under way
(1,748 genomes counted as of 2026-09-25). All of that work — the pipeline that drives the cohort
through the method, the pilot, the counts files, the page built from them, the evidence
write-up, and the comparisons with ddPCR, with HPRC assemblies and with published estimates —
lives in its own repository, [NGS-DOSE-1000G](https://github.com/jlanej/NGS-DOSE-1000G), with the
page live at [jlanej.github.io/NGS-DOSE-1000G](https://jlanej.github.io/NGS-DOSE-1000G/). In
short, on the run's first 735 genomes and 149 trios (the page carries the current count):
sequence of known copy number reads at its known copy number in every genome; the targeted fetch
reproduces the whole-file scan; the 45S copy number is inherited with a reliability near 1 in the
trios while the culture's and the library's properties are not; the calibrated estimate
reproduces across sequencing technologies where a read-depth ratio does not; and against ddPCR it
reads about 0.96× the assay. The numbers, and what each finding rules out, are in that
repository's `docs/EVIDENCE.md`. This repository holds the method alone, so that it can be applied
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
resource bundle (`export NGSDOSE_RESOURCES=/path/to/resources/GRCh38`).

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
bad arguments 2.

```bash
ngsdose estimate NA12878.json.gz -o estimates/ -t single_sample.tsv   # writes estimates/NA12878.estimate.json.gz
```

`estimate` checks each counts file against the bundle. A class it cannot measure (no sinks in
the fetch, sinks on contigs the input lacked, a different panel) is NA in the table, with the
reason in `<class>.status`, and a warning; a file that cannot be read or fitted is listed, the
table holds the others, and the exit status is 1.

```bash
# cohort: window calibration, coverage-PC adjustment, transmission reliability
ngsdose cohort estimates/*.estimate.json.gz -t cohort.tsv --save-efficiencies efficiencies.json
ngsdose adjust cohort.tsv --pcs ngspca/svd.pcs.txt -c rDNA45S.cn -o cohort.adjusted.tsv   # PCs above the Marchenko-Pastur edge; --n-pc N overrides
ngsdose pcsweep cohort.tsv --pcs ngspca/svd.pcs.txt -c rDNA45S.cn -p pedigree.txt -o sweep.tsv   # known-truth error and transmission for every number of PCs
ngsdose trios cohort.adjusted.tsv -p pedigree.txt -c rDNA45S.cn rDNA45S.cn_single rDNA45S.18S.flat \
    --compare-to rDNA45S.18S.flat       # reliabilities with bootstrap CIs, and paired differences
```

A cohort takes one estimate per sample (`cohort` refuses a sample given twice). `trios` and
`pcsweep` centre values within the populations of the pedigree's population column, or of a
`--population` file, and warn when the samples are not labelled with at least two populations of
five or more (a cohort whose samples all share one label is not warned about); the bootstrap resamples
whole families.

```bash
ngsdose selftest        # simulation checks of the statistics; needs no data
```

A cohort's validation page - known truths in every sample, fetch against scan, transmission in
trios, the coverage PCs - is built from these outputs by NGS-DOSE-1000G's `report` package; it is
written for that cohort, and shows what any cohort's page needs.

```bash
# whole-file scan, the full-accuracy mode: placement-independent, the mode for the experimental
# satellite families wherever no satellite sinks have been learned (none ship), and a record of
# where class reads were aligned and what else is in those 1-kb bins
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
Check capture with `--evaluate` on scans the sinks were not learned from.

A fetch refuses a loaded panel class that has no interval in the sinks, or whose intervals all
lie on contigs the input's header lacks: its reads would be counted only where they fall inside
other intervals. `--allow-missing-sinks` records the gap in the counts (`sinks_missing_classes`)
instead. Sink intervals on contigs the header lacks are always left out, with a warning, and
recorded per class (`sinks_skipped`); `ngsdose estimate` reports such a class as NA rather than
as an undercount.

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
alignments can be told apart.

The 1000 Genomes cohort run - its pipeline for SLURM or a plain loop, its pilot, its counts files
and the page built from them - is [NGS-DOSE-1000G](https://github.com/jlanej/NGS-DOSE-1000G),
published as the run proceeds. `resources/build/` rebuilds the GRCh38 bundle from public inputs,
pinned by sha256, into `work/bundle/GRCh38`, and records in `build_manifest.tsv` whether each
output is identical to the shipped file.

## Output columns (per sample)

| column | meaning |
| --- | --- |
| `rDNA45S.cn` | cohort-calibrated diploid copy number (`ngsdose cohort`); NA for a sample whose estimate lacks the class or has no usable window |
| `rDNA45S.cn_single` | single-sample headline: anchor windows under the fragment-GC model, or every usable window when the anchors hold fewer than 1,000 fragment ends (`.cn_basis` says which, `.n_anchor` how many ends) |
| `rDNA45S.18S`, `.28S`, … / `.flat` | per-feature estimates with / without the GC model; `18S.flat` is the estimator used in the UK Biobank literature |
| `rDNA5S.cn`, `DJ.cn` | 5S units; distal junction (expected 10) |
| `HSat3.mass_Mb`, `ACRO.mass_Mb`, `TEL.mass_Mb`, … | the experimental panels: diploid sequence mass of a satellite family (scan mode, or a fetch with satellite sinks learned for the pipeline; none ship) or of the telomeric repeat (either mode; its sinks are in the bundle). How far each is validated: `resources/experimental/README.md` |
| `rDNA45S.status`, … | `ok`, or why the class is NA: `no_sinks_in_fetch`, `sinks_skipped`, `panel_mismatch`, or `skipped: <reason>` for a positional class the bundle has no panel entry or unit for |
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
| `ctrlPC1…`, `ctrlPC_mp` | components of the control regions' residual depth across the cohort: internal technical covariates, usable by `ngsdose adjust` when NGS-PCA has not been run; `ctrlPC_mp` is how many of them stand above the noise edge (what `adjust` uses by default). None below 10 samples with control residuals; NA for a sample without them |
| `gc_curve_max_se` | how well the sample's GC curve is determined (large at very low depth) |

## Status

Engine, estimator, cohort layer and the GRCh38 bundle (45S, 5S, DJ) are implemented and tested:
Rust unit tests, a simulated genome with known truth run end to end, a 2% subsample of real
NA12878 reads, a mock trio cohort (that subsample sixty times over) through the whole cohort
layer and bundle-integrity checks, all in CI (Linux with Python 3.10, 3.12 and 3.13; macOS with
Python 3.12). Every push to `main` publishes the container image once its smoke test on the
NA12878 subsample passes, and a version tag makes a release with prebuilt engines, the Python
package and the resource bundle (`.github/workflows/`). Every assumption the counts files rest
on was audited against data before the cohort run, and reviews extended the audit during it;
DESIGN.md §15 lists what held and what was wrong, each with its fix. Being validated on
the 1000 Genomes cohort in NGS-DOSE-1000G (a twelve-genome pilot with independent library
replicates, then the cohort run: trios, two counting modes, ddPCR, assemblies, published
estimates); its results page covered 1,259 genomes and 252 complete trios as of 2026-09-24.
Experimental panels - ten satellite families, measured by scan (no satellite sinks ship), and
the telomeric repeat, which the bundle's sinks make fetchable - ship under
`resources/experimental/`. Against HPRC release-2 assemblies of 96 cohort members (the results
page as of 2026-09-24), the median estimate/assembly ratio is 0.94–1.00 for HSat3 and the
α-satellite HORs, 0.85–0.87 for HSat1A and HSat1B, and 0.43–0.79 for CER, β-satellite and ACRO,
which read low by their k-mer recall. Across people HSat1B tracks the assemblies at r = 0.99, and
ACRO, β-satellite, CER and HSat1A at 0.87–0.95; the HORs reach only 0.65 because people differ by
about 4%. HSat2 does not track the assemblies with at most 2% of its arrays in marked gaps (r = 0.09 in 47;
0.18 in the 43 with none) although it is inherited.
(Recomputed on the same table with the revised gap accounting, which also counts standalone gap
records next to an array, the HORs give r = 0.63 in 92 genomes and HSat2 0.02 in 38; the page
shows these once it is regenerated.) SST1 and SATR are annotated so differently in HPRC and in
CHM13 that their absolute ratios mean nothing. `resources/experimental/README.md` says how far
each is validated, and NGS-DOSE-1000G holds the current numbers.
Not yet done: the rest of the cohort run (1,748 of 3,202 genomes counted, and 385 of the 602
trios complete among them, as of 2026-09-25); sinks for DRAGEN-aligned data (UK Biobank, All of
Us), which must be learned from whole-file scans of a subset of those CRAMs; shipped satellite
sinks, without which the satellite families are measured only where a whole file is scanned - at
biobank scale, a subset; a wider orthogonal rDNA calibration than the twelve ddPCR lines. See
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
