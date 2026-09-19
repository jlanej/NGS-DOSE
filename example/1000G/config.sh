#!/usr/bin/env bash
# Shared configuration for the 1000 Genomes 30x example. Override any variable by exporting it.
EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$EX_DIR/../.." && pwd)"

WORK_DIR="${WORK_DIR:-/scratch/${USER}/ngs_dose_1000G}"
LOG_DIR="${LOG_DIR:-$WORK_DIR/logs}"
MANIFEST="${MANIFEST:-$WORK_DIR/manifest.tsv}"      # SAMPLE <TAB> CRAM (https:// URL or local path)

# the engine and the python package: local installs (cargo build --release; pip install .), or
# the container for both (`apptainer exec $SIF ...`) when SIF is set
NGSDOSE_BIN="${NGSDOSE_BIN:-$REPO/target/release/ngs-dose}"
SIF="${SIF:-}"
ngsdose_py() { if [ -n "$SIF" ]; then apptainer exec "$SIF" "$@"; else "$@"; fi; }
BUNDLE="${BUNDLE:-$REPO/resources/GRCh38}"
# panels to load; the satellite panel is experimental and only meaningful in scan mode
SATELLITES="${SATELLITES:-$REPO/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz}"

# CRAM decoding needs the reference: a local FASTA (recommended) or htslib's REF_PATH/REF_CACHE
REF_FASTA="${REF_FASTA:-$WORK_DIR/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}"

# fetch = controls + learned sinks through the index: ~0.5 GB and ~1 min per sample straight from
#         the public bucket, no CRAM on disk. 45S, 5S, DJ, the known-truth regions, chrM and chrEBV.
#         The whole cohort is ~1.5 TB of transfer - against ~48 TB for whole files - which is why
#         it is the default here: on most clusters the WAN, not the CPU, is the scarce resource.
# scan  = read every record: placement-independent, measures the (experimental) satellite classes
#         too, and is what sinks are learned from and checked against. ~14 CPU-min per 30x genome
#         and 15 GB per sample over HTTPS. Run it on a subset: 00_setup.sh writes manifest.hprc.tsv,
#         the 200 samples with HPRC assemblies (03_compare_modes.sh, 04_hprc_satellites.sh).
# Counts go to $WORK_DIR/counts_$MODE, so both can be run side by side.
MODE="${MODE:-fetch}"
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
