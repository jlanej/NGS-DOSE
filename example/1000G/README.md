# 1000 Genomes 30× example

The expanded 1000 Genomes cohort (3,202 samples, 602 trios; NYGC, TruSeq PCR-free, NovaSeq
2×150, GRCh38) is the validation cohort for NGS-DOSE: it is public, it has trios, many samples
have HPRC assemblies, and NGS-PCA has already been run on it (`ngspca/`).

| path | what |
| --- | --- |
| `pilot/` | four trios, each sample also as an independent older library; run locally; counts files, the evaluation and plotting scripts, the report and the figure are here |
| `00_setup.sh`, `01_count.sh`, `02_cohort.sh`, `03_compare_modes.sh`, `04_hprc_satellites.sh`, `config.sh` | the full-cohort run, for SLURM or a plain loop |
| `ngs-dose.def` | fallback Apptainer definition of the container image |
| `hprc_r2_censat.keys.txt`, `hprc_satellites.py` | the HPRC release-2 CenSat annotations (S3 keys; 205 samples, 200 of them in this cohort) and the comparison of satellite estimates against them |
| `ngspca/` | NGS-PCA output for this cohort (200 coverage PCs, `AUTO_HQ_median`), produced by [NGS-PCA's 1000G example](https://github.com/jlanej/NGS-PCA/tree/master/example/1000G_highcov) |

## The pilot

```bash
REF=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa bash pilot/run_pilot.sh
```

counts every sample in fetch mode straight from the public CRAMs (AWS Open Data mirror for the
NYGC data, EBI for the HGSVC replicates) and writes [`pilot/pilot_report.md`](pilot/pilot_report.md).
`python pilot/evaluate_pilot.py` alone regenerates the report from the committed counts files.

The replicates are what make the pilot informative. HG00512/3/4, HG00731/2/3 and NA19238/39/40
were sequenced years earlier for the HGSVC (HiSeq 2500, 2×126, ~78×, bwakit + postalt), and the
CEU trio for Illumina's Platinum pedigree (HiSeq 2000, 2×101, ~54×): different library
preparation, chemistry, read length, insert size, depth and alignment pipeline, and a GC
response that is the reverse of NYGC's. Agreement between the two is reproducibility of the
*measurement*, not of the file — and an upper bound on its error, because the two DNA batches
come from different cultures of the cell line.

`python pilot/evaluate_pilot.py --write-anchors` additionally rewrites
`resources/GRCh38/anchors.json` from the replicate pairs; `python pilot/plot_pilot.py` draws
`pilot_figure.png` (needs matplotlib).

## The cohort

```bash
export WORK_DIR=/scratch/$USER/ngs_dose_1000G
cd example/1000G                                        # submit from here: config.sh and logs/ are relative to it
bash 00_setup.sh                                        # reference, pedigree, manifests (all 3,202; the 200 with HPRC assemblies)
N=$(wc -l < $WORK_DIR/manifest.tsv)
sbatch --array=0-$(( (N - 1) / 10 ))%25 01_count.sh     # fetch mode, the whole cohort: ~0.5 GB and ~1 min per sample
sbatch 02_cohort.sh                                     # estimate, calibrate, adjust, transmission

N=$(wc -l < $WORK_DIR/manifest.hprc.tsv)                # whole-file scans of the 200 samples with HPRC assemblies
MODE=scan MANIFEST=$WORK_DIR/manifest.hprc.tsv sbatch --array=0-$(( (N - 1) / 10 ))%10 01_count.sh
bash 03_compare_modes.sh                                # what the sinks capture in every scanned sample; sinks re-learned
sbatch 04_hprc_satellites.sh                            # satellite array mass against the assemblies of the same people
```

**Fetch the cohort, scan a subset.** Fetch mode reads the control regions and the class sinks
through the index, straight from the public bucket: about 1.5 TB of transfer for all 3,202
samples, nothing staged on disk, and everything the transmission analysis needs (45S, 5S, DJ,
the known-truth regions, chrM and chrEBV). A scan reads every record — 15 GB per sample, 48 TB
for the cohort, the download campaign that NGS-PCA's example had to engineer around — and what
it adds is validation and the experimental classes: it is placement-independent, it records
where every class read was aligned (so `03_compare_modes.sh` can say how much of each class the
shipped sinks capture in each person, across populations and both sexes), and it is the only
mode that measures the satellite families, for which the 200 HPRC assemblies are a truth
(`04_hprc_satellites.sh`). `MODE=scan` on the full manifest works too, if the bandwidth is there.

Operational notes:

- `01_count.sh` is idempotent (finished samples are skipped), tries each sample three times,
  and exits non-zero if any sample still failed; re-submit to retry those. The engine retries
  failed intervals itself, gives up with exit status 75 when a connection has gone silent for
  five minutes (a dead HTTPS connection waits for ever rather than failing; `timeout` is a
  backstop where it exists, and `--stall-timeout 0` turns the watchdog off, e.g. for files that
  have to be recalled from tape), writes its output under a temporary name and renames it, and
  refuses a BAM/CRAM without an end-of-file marker (a truncated copy decodes without error and
  loses exactly the contigs the rDNA is on).
  `%25` caps concurrent array tasks: the bottleneck is the network, and hundreds of parallel
  readers against one S3 prefix earn throttling rather than speed.
- The manifest's second column may be a local path instead of a URL (e.g. CRAMs staged by
  NGS-PCA's download stage); fetch mode then looks for `<cram>.crai` next to it.
- With `export SIF=/path/to/ngs-dose.sif` both the engine and the Python steps run through
  Apptainer. `00_setup.sh` pulls the image if the file is not there yet
  (`docker://ghcr.io/jlanej/ngs-dose:latest`, published by `.github/workflows/container.yml` on
  version tags or on demand; a private package needs `apptainer remote login` first). Where
  nothing has been published, `ngs-dose.def` builds the same image from a checkout
  (`apptainer build --fakeroot`, from the repository root). Bind the repository and `$WORK_DIR` if
  your site does not do so by default (`export APPTAINER_BIND=$WORK_DIR,$PWD/../..`). Without
  it: `cargo build --release` (needs libclang, e.g. `module load llvm`) and `pip install .`.
- Every counts file records the SHA-256 of the panel, controls and sinks it was made with and
  the lengths of the contigs it used; `ngsdose estimate` warns if a cohort mixes resource sets
  and refuses a file aligned to another reference build.
- QC columns to look at first in `cohort.tsv`: `truth.auto` (2), `truth.chrX` and `truth.chrY`
  (sex; 1.6 and the like is mosaic loss, common in LCLs), `DJ.cn` (10), `flagged_chromosomes`
  (aneuploidy), `eof_marker`, `gc_curve_max_se`.

### What the cohort run is for

`02_cohort.sh` ends with the table the design is waiting on — transmission reliability of every
candidate estimator (`rDNA45S.18S.flat`, the estimator in the literature; the single-sample
anchor estimate; the calibrated `cn`; 5S; DJ) and of three traits whose answer is known in
advance (`truth.auto`: no variance but error; `chrEBV.copies` and `chrM.copies`: large variance,
none of it transmitted through the nuclear genome - if these come out "reliable", families
share batches and every other reliability in the table is inflated by as much), unadjusted, adjusted on NGS-PCA's coverage PCs and
adjusted on the internal control PCs, each with a family-bootstrap confidence interval, the
spousal correlation, a permuted-family null, and the *paired* bootstrap of every estimator
against the 18S depth ratio (separate intervals are about ±0.09 at 602 trios; the paired
difference is much sharper). Further checks that need the cohort and are not scripted yet:

- `DJ.cn` across 3,202 samples: a tight distribution at 10 is the accuracy claim; integer
  outliers (8, 9, 11) are candidate acrocentric rearrangements to look at;
- `truth.auto` and `truth.chrX` by sex as cohort-wide QC, and `flagged_chromosomes` as an LCL
  aneuploidy screen;
- sample-by-sample comparison with Hall et al. 2021 (their Supplementary Data 1 covers 2,419 of
  these samples), including the prediction that the ratio between the two pipelines tracks each
  sample's duplicate-flag rates;
- whether the leading coverage PCs explain estimate variance, and whether reliability rises or
  falls when they are removed;
- inherited one-copy steps along the DJ (seen in three of the four pilot trios) as a Mendelian
  test on integer states: a step in a child and in neither parent is a de novo event or an
  error, and a parental step is transmitted half the time;
- whether DJ, the female X and the leading control PC move together (the S-phase hypothesis);
- whether a child's departure from the midparent tracks the state of the culture it was
  sequenced from (`chrEBV.copies`, `chrM.copies`, the leading control PCs): the part of the
  non-transmitted variance that is biology of the cell line rather than measurement error;
- the two published claims that need exactly this cohort and these controls: a 5S-45S
  correlation (Gibbons et al. 2015; not seen by Hall et al. 2021) and an inverse relation
  between rDNA and mitochondrial DNA abundance (Gibbons et al. 2014) - both now testable within
  one library type, against known-truth regions, with transmission as the arbiter of what is
  signal.

## Data sources and attribution

- NYGC 30× CRAMs: Byrska-Bishop et al., *Cell* 185:3426 (2022); AWS Open Data mirror `s3://1000genomes/1000G_2504_high_coverage/`.
- HGSVC and Illumina Platinum-pedigree CRAMs: IGSR data collections `hgsv_sv_discovery` and `illumina_platinum_pedigree` (EBI).
- `pilot/hall2021_MOESM1.txt`: Supplementary Data 1 of Hall, Turner & Queitsch, *Sci Rep* 11:449 (2021),
  doi:10.1038/s41598-020-80049-y, distributed under CC BY 4.0; unmodified.
