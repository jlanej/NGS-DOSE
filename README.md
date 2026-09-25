# NGS-DOSE

Sequence-class dosage from short-read WGS: how many 45S rDNA units, 5S units, and copies of
other multi-copy classes a person carries, measured from the GRCh38 CRAMs that biobanks already
hold — in about a minute per genome, without downloading it.

**NGS-DOSE is not a CNV caller.** It measures the part of the genome that variant and CNV callers
mask. [NGS-PCA](https://github.com/jlanej/NGS-PCA) describes a cohort's coverage from the bins it
retains; NGS-DOSE measures what lives in the bins it excludes.

## How it works

```
CRAM ──ngs-dose count──▶ counts.json ──ngsdose estimate──▶ per-sample CN ──ngsdose cohort / adjust / trios──▶ cohort table
        (Rust, htslib)    70-250 kB         (Python, numpy)
```

1. **Count fragment 5′ ends**, in single-copy control regions and in every read that carries
   class-diagnostic 31-mers. A read is assigned to a class, and placed on its unit, by k-mers —
   not by where the aligner put it. `scan` reads the whole file; `fetch` retrieves only the
   control regions and the few *sinks* where a class's reads are known to land (learned from
   scans), and returns the same counts to within what the sinks miss: 0.05% of 45S reads, 0.02%
   of 5S, 0.2% of the distal junction in the samples scanned so far.
2. **Model the library**: a Poisson spline of end density on fragment-scale GC, fitted per sample
   on the controls. Copy number of a window of the unit is `2 × observed / expected`.
3. **Calibrate the unit**: parts of the rDNA unit drop out far beyond what any genome-wide GC
   curve predicts (28S reads ~0.75× of 18S on NovaSeq), and *which* parts depends on the
   sequencing chemistry. Per-window efficiencies are learned per library type, with the scale
   pinned on *anchor* windows where different chemistries were shown to agree.
4. **Check it against known truth in every sample**: held-out autosomal sequence (2 copies),
   chrX (1 or 2), chrY (1 or 0), and the acrocentric distal junction (10 copies), measured by
   the same code paths as the classes. Mitochondrial genomes and EBV episomes per cell come
   along as covariates of the tissue or culture the DNA was taken from.
5. **Cohort layer**: adjustment on coverage PCs (NGS-PCA's, or the control regions' own) - as
   many as clear the noise edge of their spectrum, with a sweep against the known truths to
   confirm or overrule that number - and transmission reliability in trios to decide which
   estimator carries the most real variance.

The reasoning, and the measurements on real data behind each step, are in
[docs/DESIGN.md](docs/DESIGN.md).

## Validation on the 1000 Genomes cohort

The method is validated on the expanded 1000 Genomes cohort (3,202 genomes, 602 trios, NYGC 30×
NovaSeq CRAMs): a twelve-genome pilot in which every sample also has an older library of the same
cell line on another instrument, and the cohort run that follows it. All of that work — the
pipeline that drives the cohort through the method, the pilot, the counts files, the page built
from them, the evidence write-up, and the comparisons with ddPCR, with HPRC assemblies and with
published estimates — lives in its own repository, [NGS-DOSE-1000G](https://github.com/jlanej/NGS-DOSE-1000G),
with the page live at [jlanej.github.io/NGS-DOSE-1000G](https://jlanej.github.io/NGS-DOSE-1000G/).
In short: sequence of known copy number reads at its known copy number in every genome; the
targeted fetch reproduces the whole-file scan; the 45S copy number is inherited with a reliability
near 1 in the trios while the culture's and the library's properties are not; the calibrated
estimate reproduces across sequencing technologies where a read-depth ratio does not; and against
ddPCR it reads about 0.96× the assay. The numbers, and what each finding rules out, are in that
repository's `docs/EVIDENCE.md`. This repository holds the method alone, so that it can be applied
to any cohort.

## Quick start

From source (the engine needs a Rust toolchain and libclang; the package needs numpy):

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

and a version tag makes a [release](https://github.com/jlanej/NGS-DOSE/releases) that carries the
engine for Linux (glibc 2.28 and newer, curl and TLS built in) and for macOS on Apple silicon,
the Python wheel, and the resource bundle (`export NGSDOSE_RESOURCES=/path/to/resources/GRCh38`).

```bash
B=resources/GRCh38
# one genome, straight from the AWS Open Data mirror (~1 min; CRAM decoding needs the reference)
target/release/ngs-dose count -m fetch -@ 16 \
    -i https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram \
    -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -p $B/panel.k31.tsv.gz -c $B/controls.fa.gz --sinks $B/sinks.bed -o NA12878.json.gz
```

```bash
ngsdose estimate NA12878.json.gz -o estimates/ -t single_sample.tsv
```

```bash
# cohort: window calibration, coverage-PC adjustment, transmission reliability
ngsdose cohort estimates/*.estimate.json.gz -t cohort.tsv --save-efficiencies efficiencies.json
ngsdose adjust cohort.tsv --pcs ngspca/svd.pcs.txt -c rDNA45S.cn -o cohort.adjusted.tsv   # PCs above the Marchenko-Pastur edge; --n-pc N overrides
ngsdose pcsweep cohort.tsv --pcs ngspca/svd.pcs.txt -c rDNA45S.cn -p pedigree.txt -o sweep.tsv   # known-truth error and transmission for every number of PCs
ngsdose trios cohort.adjusted.tsv -p pedigree.txt -c rDNA45S.cn rDNA45S.cn_single rDNA45S.18S.flat \
    --compare-to rDNA45S.18S.flat       # reliabilities with bootstrap CIs, and paired differences
```

```bash
ngsdose selftest        # simulation checks of the statistics; needs no data
```

A cohort's validation page - known truths in every sample, fetch against scan, transmission in
trios, the coverage PCs - is built from these outputs by NGS-DOSE-1000G's `report` package; it is
written for that cohort, and shows what any cohort's page needs.

```bash
# whole-file scan, the full-accuracy mode: placement-independent, the only mode for dispersed sequence (the
# experimental satellite families; the telomeric repeat is fetchable), and a record of where class reads were aligned and
# what else is in those 1-kb bins
E=resources/experimental
target/release/ngs-dose count -m scan -@ 10 -i sample.cram -T ref.fa -c $B/controls.fa.gz -p $B/panel.k31.tsv.gz \
    -p $E/satellites.CHM13v2.k31.panel.tsv.gz -p $E/telomere.k31.panel.tsv.gz -o sample.scan.json.gz
```

A different aligner or decoy set places multi-copy reads differently. Before trusting `fetch`
on a new pipeline, scan a few whole CRAMs and check the sinks:

```bash
ngsdose sinks scan*.json.gz --evaluate $B/sinks.bed      # fraction of each class the sinks capture, per sample
ngsdose sinks scan*.json.gz -o sinks.bed                 # or re-learn them
```

A fetch refuses a loaded panel class that has no interval in the sinks (its reads would be counted only
where they fall inside other intervals); `--allow-missing-sinks` records the gap in the counts instead.
Where the engine cannot read the CRAMs itself, `ngs-dose plan` writes the intervals a fetch reads as BED,
to cut them out with samtools and count the result in scan mode:

```bash
target/release/ngs-dose plan -c $B/controls.fa.gz --sinks $B/sinks.bed -i sample.cram -T ref.fa -o plan.bed
samtools view -b -M -L plan.bed -o sample.fetch.bam sample.cram && samtools index -c sample.fetch.bam
```

The 1000 Genomes cohort run - its pipeline for SLURM or a plain loop, its pilot, its counts files
and the page built from them - is [NGS-DOSE-1000G](https://github.com/jlanej/NGS-DOSE-1000G),
published as the run proceeds. `resources/build/` rebuilds the GRCh38 bundle from public inputs.

## Output columns (per sample)

| column | meaning |
| --- | --- |
| `rDNA45S.cn` | cohort-calibrated diploid copy number (`ngsdose cohort`) |
| `rDNA45S.cn_single` | single-sample headline: anchor windows under the fragment-GC model |
| `rDNA45S.18S`, `.28S`, … / `.flat` | per-feature estimates with / without the GC model; `18S.flat` is the estimator used in the UK Biobank literature |
| `rDNA5S.cn`, `DJ.cn` | 5S units; distal junction (expected 10) |
| `HSat3.mass_Mb`, `ACRO.mass_Mb`, `TEL.mass_Mb`, … | the experimental panels: diploid sequence mass of a satellite family (scan mode) or of the telomeric repeat (either mode; its sinks are in the bundle). How far each is validated: `resources/experimental/README.md` |
| `*.adj` | an estimate with coverage PCs regressed out (`ngsdose adjust`) |
| `truth.auto`, `truth.chrX`, `truth.chrY` | held-out known-copy-number sequence (expected 2; 1 or 2; 1 or 0) |
| `chrM.copies`, `chrEBV.copies` | mitochondrial genomes and EBV episomes per cell: covariates of the state of the tissue or cell line, measured like the truths |
| `engine` | engine version and the commit it was built from, as recorded in the counts file |
| `eof_marker` | `present` unless the input lacked its end-of-file block (the engine refuses such files unless told otherwise) |
| `depth`, `ctrl_dup_frac`, `gc_rel_35`, `gc_rel_65`, `ctrl_region_sd`, `flagged_chromosomes` | library and sample QC; an aneuploid chromosome is reported and excluded from the denominator |
| `*.profile_sd`, `*.profilePC*` | how far the sample's window profile departs from the cohort's |
| `ctrlPC1…`, `ctrlPC_mp` | components of the control regions' residual depth across the cohort: internal technical covariates, usable by `ngsdose adjust` when NGS-PCA has not been run; `ctrlPC_mp` is how many of them stand above the noise edge (what `adjust` uses by default) |
| `gc_curve_max_se` | how well the sample's GC curve is determined (large at very low depth) |

## Status

Engine, estimator, cohort layer and the GRCh38 bundle (45S, 5S, DJ) are implemented and tested:
Rust unit tests, a simulated genome with known truth run end to end, a 2% subsample of real
NA12878 reads, a mock trio cohort (that subsample sixty times over) through the whole cohort
layer and bundle-integrity checks, all in CI
(Linux and macOS, Python 3.10 to 3.13). Every push to `main` publishes the container image, and
a version tag makes a release with prebuilt engines, the Python package and the resource bundle
(`.github/workflows/`). Before the cohort run every assumption the counts files
rest on was audited against data; nine were wrong and are fixed (DESIGN.md §15). Validated on the
1000 Genomes cohort in NGS-DOSE-1000G (a twelve-genome pilot with independent library replicates,
then the cohort run: trios, two counting modes, ddPCR, assemblies, published estimates).
Experimental panels - ten satellite families, which are dispersed and need scan mode, and the
telomeric repeat, which the aligner concentrates and either mode measures - ship under
`resources/experimental/`: in a first comparison with HPRC assemblies
of the same people (two samples) HSat3, HSat1A and the α-satellite HORs come out within 7% of
the assembly; the rest are relative measures or undecided, and that README says which and why.
Not yet done: the rest of the cohort run (1,259 of 3,202 genomes counted as of 2026-09-24),
sinks for DRAGEN-aligned data, a wider orthogonal rDNA calibration than the twelve ddPCR lines. See DESIGN.md §12–13. No licence has been chosen yet.

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

Two honesty notes. The transmission-reliability framing is standard midparent regression applied
to a trait whose heritability is one by construction; no claim of novelty is made for it, and no
systematic literature search backs the claims of novelty that *are* made in DESIGN.md §14. And an
AI system is not an author under prevailing journal policy — if this becomes a paper, the
appropriate form is a contributions or acknowledgements statement describing what was
machine-generated, with human authors taking responsibility for verification.
