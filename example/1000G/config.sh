#!/usr/bin/env bash
# Shared configuration for the 1000 Genomes 30x example. Override any variable by exporting it.
EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$EX_DIR/../.." && pwd)"

WORK_DIR="${WORK_DIR:-/scratch/${USER}/ngs_dose_1000G}"
LOG_DIR="${LOG_DIR:-$WORK_DIR/logs}"
MANIFEST="${MANIFEST:-$WORK_DIR/manifest.tsv}"      # SAMPLE <TAB> CRAM (https:// URL or local path)

# Two ways to run. With SIF set, everything - engine, python package, resource bundle, aria2c -
# comes from the container image and nothing is installed on the host (these scripts themselves
# can be copied out of the image: `apptainer exec $SIF ngs-dose-example DIR`). Without it: a
# checkout with `cargo build --release` and `pip install .`.
SIF="${SIF:-}"
NGSDOSE_BIN="${NGSDOSE_BIN:-$REPO/target/release/ngs-dose}"
if [ -n "$SIF" ]; then IMAGE_ROOT=/opt/ngs-dose; else IMAGE_ROOT="$REPO"; fi
BUNDLE="${BUNDLE:-$IMAGE_ROOT/resources/GRCh38}"
# panels to load; the satellite panel is experimental and only meaningful in scan mode
SATELLITES="${SATELLITES:-$IMAGE_ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz}"
# staged CRAMs for the full-accuracy path (01b_dose_sample.sh): $CRAM_DIR/<sample>.cram(.crai)
CRAM_DIR="${CRAM_DIR:-$WORK_DIR/crams}"

# CRAM decoding needs the reference: a local FASTA (recommended) or htslib's REF_PATH/REF_CACHE
REF_FASTA="${REF_FASTA:-$WORK_DIR/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}"

# what the container has to see (bind sources must exist: 00_setup.sh creates them)
APPTAINER_BINDS="${APPTAINER_BINDS:-$WORK_DIR,$EX_DIR,$(dirname "$REF_FASTA"),$CRAM_DIR}"
ngsdose_engine() { if [ -n "$SIF" ]; then apptainer exec --bind "$APPTAINER_BINDS" "$SIF" ngs-dose "$@"; else "$NGSDOSE_BIN" "$@"; fi; }
# any other command of the image (ngsdose, python3, aria2c, test): run where the bundle paths are valid
ngsdose_py() { if [ -n "$SIF" ]; then apptainer exec --bind "$APPTAINER_BINDS" "$SIF" "$@"; else "$@"; fi; }
have_file() { ngsdose_py test -s "$1"; }

# scan  = read every record: the full-accuracy mode. Placement-independent, measures the
#         (experimental) satellite classes too, records where every class read was aligned and
#         what else sits in those 1-kb bins, and is what sinks are learned from and checked
#         against. ~14 CPU-min per 30x genome; 15 GB per sample, so stage the CRAMs
#         (01b_dose_sample.sh) rather than stream them if the whole cohort is to be scanned.
# fetch = controls + learned sinks through the index: ~0.5 GB and ~1 min per sample straight from
#         the public bucket (3 s from a local file), no CRAM on disk: 45S, 5S, DJ, the known-truth
#         regions, chrM and chrEBV. The biobank-scale mode; 01b_dose_sample.sh runs it on every
#         staged CRAM as well, so the cohort says what it costs in accuracy, sample by sample.
# Counts go to $WORK_DIR/counts_$MODE; MODE selects which of them 01_count.sh makes and
# 02_cohort.sh analyses.
MODE="${MODE:-scan}"
COUNTS_DIR="${COUNTS_DIR:-$WORK_DIR/counts_$MODE}"
EST_DIR="${EST_DIR:-$WORK_DIR/estimates_$MODE}"
THREADS="${THREADS:-8}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-10}"
SAMPLE_TIMEOUT="${SAMPLE_TIMEOUT:-5400}"            # seconds; a stalled HTTPS connection can hang
S3_HTTPS_BASE="${S3_HTTPS_BASE:-https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage}"

# NGS-PCA outputs for this cohort (coverage PCs); the copy in this repository is the default
NGSPCA_DIR="${NGSPCA_DIR:-$EX_DIR/ngspca}"
PEDIGREE="${PEDIGREE:-$WORK_DIR/20130606_g1k_3202_samples_ped_population.txt}"
PEDIGREE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
