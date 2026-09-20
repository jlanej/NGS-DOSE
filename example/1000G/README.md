# 1000 Genomes 30× example

The expanded 1000 Genomes cohort (3,202 samples, 602 trios; NYGC, TruSeq PCR-free, NovaSeq
2×150, GRCh38) is the validation cohort for NGS-DOSE: it is public, it has trios, many samples
have HPRC assemblies, and NGS-PCA has already been run on it (`ngspca/`).

| path | what |
| --- | --- |
| `pilot/` | four trios, each sample also as an independent older library; run locally; counts files, the evaluation and plotting scripts, the report and the figure are here |
| `00_setup.sh`, `01_stage_and_dose.sh`, `01a_dispatch_staged.sh`, `01b_dose_sample.sh`, `01_count.sh`, `02_cohort.sh`, `03_compare_modes.sh` (+ `compare_modes.py`), `04_hprc_satellites.sh`, `config.sh` | the full-cohort run, for SLURM or a plain loop |
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

Full accuracy for everyone: every CRAM is staged, scanned whole, fetched as well (three seconds
more on a local file), verified, and removed.

Nothing is installed on the cluster: the image holds the engine, the `ngsdose` package, the
resource bundle, `aria2c` and these scripts, and the host needs Apptainer, SLURM, `curl` and the
coreutils. (From a checkout instead: leave `SIF` unset, `cargo build --release`, `pip install .`.)

```bash
apptainer pull ngs-dose.sif docker://ghcr.io/jlanej/ngs-dose:latest
apptainer exec ngs-dose.sif ngs-dose-example $PWD/dose1000G    # the scripts submit jobs, so they live on the host
cd dose1000G                                            # submit from here: config.sh and logs/ are relative to it
export SIF=$PWD/../ngs-dose.sif WORK_DIR=/scratch/$USER/ngs_dose_1000G

bash 00_setup.sh                                        # reference, pedigree, manifest (sample, URL, MD5)
LIMIT=5 bash 01_stage_and_dose.sh                       # smoke test: five samples staged, counted, removed
sbatch 01_stage_and_dose.sh                             # the cohort: a few transfers at a time, one 01b_dose_sample.sh job per CRAM

MODE=scan sbatch 02_cohort.sh                           # estimate, calibrate, adjust, transmission - on the scans
MODE=fetch sbatch 02_cohort.sh                          # the same on the targeted fetches, for the comparison
bash 03_compare_modes.sh                                # sink capture per sample, sinks re-learned, fetch / scan per sample
sbatch 04_hprc_satellites.sh                            # satellite array mass against the HPRC assemblies of the same people
```

`01_stage_and_dose.sh` keeps a few multi-connection `aria2c` transfers going (the image's
`aria2c`; it checks each file against the MD5 of the sequence index), holds at most
`MAX_LOCAL_CRAMS` files on disk, submits one job per landed CRAM and can be re-run at any time:
finished samples, queued jobs and files already on disk are recognised. Compute nodes need no
internet unless the manager itself runs as a job; on a cluster where they have none, run it on a
login or transfer node (`bash 01_stage_and_dose.sh` in tmux).

`01b_dose_sample.sh SAMPLE LINE [MD5]` has the calling convention of NGS-PCA's
`01b_mosdepth_sample.sh` (verify the MD5 if given, process, check the outputs, delete the CRAM),
so the aria2 download manager that staged this cohort for NGS-PCA can drive it: it needs its
per-sample job script and its "already done" test to be settable (`$WORK_DIR/counts_scan/<sample>.json.gz`
instead of the mosdepth output), nothing else. `01a_dispatch_staged.sh` is for files staged by
other means (Globus, rsync): it neither downloads nor deletes, skips files aria2 is still
writing, and submits each landed CRAM once. About 14 CPU-minutes and 1.3 GB of memory per sample; the counts files
are ~230 kB (scan) and ~70 kB (fetch), so the whole cohort is under 1 GB.

Without staging, `01_count.sh` counts straight from the public bucket, a block of the manifest
per array task: `MODE=fetch` moves ~0.5 GB per sample (1.5 TB for the cohort, about a minute
each) and is how a biobank would run; `MODE=scan` streams the whole file (15 GB per sample over
one connection):

```bash
N=$(wc -l < $WORK_DIR/manifest.tsv)
MODE=fetch sbatch --array=0-$(( (N - 1) / 10 ))%25 01_count.sh
```

**What the full scan buys, and why this cohort should have it.** A scan is placement-independent;
it is the only mode that measures the satellite families (200 samples of the cohort have HPRC
assemblies to hold them against); and it records, on the 1-kb grid of `mosdepth --by 1000`,
where every class read was aligned *and every other read in those bins* - which is the data
that decides how far the measurement can be simplified for a biobank:

1. **Targeted fetch** (exists): with both modes run on every sample, `03_compare_modes.sh` gives
   fetch / scan per sample across 26 populations and both sexes, and re-learns the sinks at
   n = 3,202 on the finer grid (tighter sinks, fewer requests).
2. **A depth proxy from bins a cohort already has.** NGS-PCA's mosdepth run over this cohort
   (1-kb bins, all contigs, no MAPQ filter, duplicates excluded) is kept. In NA12878, 99.8% of
   the aligned 45S reads sit in 273 of those bins. Whether a fixed set of such bins, normalised
   to the autosomal median, tracks the full estimate - and at what loss of transmission
   reliability, which `ngsdose trios --compare-to` measures directly - needs exactly three
   things: the scans, the kept mosdepth outputs, and the per-bin census that says how much of
   each bin's depth is the class and how its duplicate-flag rate differs from the genome's
   (5.5% against 10.8% in NA12878: a depth tool that drops flagged reads over-reads rDNA by the
   difference, by a different amount in every sample). All three exist after this run.
3. Anything else - fewer or smaller sinks, fewer controls, lower depth - can be tried on the
   counts files alone, because they hold positions, not summaries.

Operational notes:

- `01b_dose_sample.sh` is idempotent under a requeue, verifies that both counts files exist and
  are whole before it removes a CRAM (`KEEP_CRAMS=1` keeps them), and leaves the CRAM in place
  if anything failed. `01_count.sh` is idempotent too (finished samples are skipped), tries each sample three times,
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
  every push to `main`, and as `:X.Y.Z` on version tags; a package that has not been made public
  needs `apptainer remote login` first). Where
  nothing has been published, `ngs-dose.def` builds the same image from a checkout
  (`apptainer build --fakeroot`, from the repository root). The scripts bind `$WORK_DIR`, the
  repository, the reference directory and `$CRAM_DIR` (`APPTAINER_BINDS` in `config.sh`).
  Without a container: `cargo build --release` (needs libclang, e.g. `module load llvm`) and
  `pip install .`.
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
difference is much sharper). The number of PCs is the Marchenko–Pastur default, separately for the two PC
sets (`N_PC=20 N_CTRL_PC=5 sbatch 02_cohort.sh` overrides them), and `pcsweep.*.tsv` with its `.summary.txt` is the evidence for or
against that default: for every number of PCs, the cross-validated error of the known truths
and the reliability of the classes. Further checks that need the cohort and are not scripted yet:

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
