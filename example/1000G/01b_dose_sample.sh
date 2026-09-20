#!/usr/bin/env bash
# =============================================================================
# 01b_dose_sample.sh - NGS-DOSE for one sample whose CRAM is already local
# =============================================================================
# The full-accuracy path: a whole-file scan of a staged CRAM (every class, the
# satellite panel, where every class read was aligned and what else sits in those
# 1-kb bins), then - three seconds more - a targeted fetch of the same file, so
# that every sample of the cohort says exactly what the fast mode would have
# returned for it. Outputs are verified, and only then is the CRAM removed.
#
# Same calling convention as NGS-PCA's 01b_mosdepth_sample.sh, so whatever stages
# CRAMs for that pipeline (its aria2 download manager, its stage watcher, a Globus
# batch) can dispatch this job instead:
#   sbatch 01b_dose_sample.sh SAMPLE_ID MANIFEST_LINE_NUM [CRAM_MD5]
# Expects $CRAM_DIR/SAMPLE_ID.cram and .cram.crai. KEEP_CRAMS=1 leaves them in place.
# =============================================================================
#SBATCH --job-name=ngsdose_sample
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dose_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
SAMPLE="${1:?usage: 01b_dose_sample.sh SAMPLE_ID MANIFEST_LINE_NUM [CRAM_MD5]}"
LINE_NUM="${2:-0}"; CRAM_MD5="${3:-}"
CRAM="$CRAM_DIR/$SAMPLE.cram"; CRAI="$CRAM.crai"
SCAN_OUT="$WORK_DIR/counts_scan/$SAMPLE.json.gz"; FETCH_OUT="$WORK_DIR/counts_fetch/$SAMPLE.json.gz"
mkdir -p "$WORK_DIR/counts_scan" "$WORK_DIR/counts_fetch"
echo "== NGS-DOSE: $SAMPLE (manifest line $LINE_NUM), started $(date)"

tidy() { [ "${KEEP_CRAMS:-0}" = 1 ] || { rm -f "$CRAM" "$CRAI"; echo "   removed $CRAM"; }; }
valid() { [ -s "$1" ] && gzip -t "$1" 2>/dev/null; }          # a counts file that exists and is whole

if valid "$SCAN_OUT" && valid "$FETCH_OUT"; then echo "   already done"; tidy; exit 0; fi   # idempotent under a requeue
[ -s "$CRAM" ] && [ -s "$CRAI" ] || { echo "ERROR: $CRAM or its index is missing" >&2; exit 1; }
if [ -n "$CRAM_MD5" ]; then
  got="$( (md5sum "$CRAM" 2>/dev/null || md5 -r "$CRAM") | awk '{print $1}')"
  [ "$got" = "$CRAM_MD5" ] || { echo "ERROR: MD5 mismatch for $CRAM ($got, expected $CRAM_MD5)" >&2; exit 1; }
  echo "   MD5 verified"
fi

common=(-i "$CRAM" --index "$CRAI" -T "$REF_FASTA" -c "$BUNDLE/controls.fa.gz" -p "$BUNDLE/panel.k31.tsv.gz" -s "$SAMPLE" -@ "${SLURM_CPUS_PER_TASK:-$THREADS}")
if ! valid "$SCAN_OUT"; then
  scan=("${common[@]}" -m scan)
  # a scan without a panel that was asked for is a scan to be done again: stop rather than skip it
  for x in $EXTRA_PANELS; do have_file "$x" || { echo "ERROR: extra panel $x not found" >&2; exit 1; }; scan+=(-p "$x"); done
  ngsdose_engine count "${scan[@]}" -o "$SCAN_OUT"
fi
valid "$FETCH_OUT" || ngsdose_engine count "${common[@]}" -m fetch --sinks "$BUNDLE/sinks.bed" -o "$FETCH_OUT"
valid "$SCAN_OUT" && valid "$FETCH_OUT" || { echo "ERROR: counts missing or damaged for $SAMPLE; the CRAM is kept" >&2; exit 1; }
tidy
echo "== $SAMPLE complete, $(date)"
