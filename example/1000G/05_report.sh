#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_report
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/report_%j.out
# The cohort page, from whatever counts exist so far: known truths in every sample, fetch against
# scan, transmission through the trios, the cell-line covariates, satellites against assemblies.
# Re-run at any time; a sample whose counts have not changed is not estimated again.
#   sbatch 05_report.sh                      -> $WORK_DIR/report/index.html (+ report.json, data/)
#   REPORT_OUT=/path/to/NGS-DOSE-1000G/docs bash 05_report.sh   -> straight into the results repository's Pages folder
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
OUT="${REPORT_OUT:-$WORK_DIR/report}"; mkdir -p "$OUT"
args=(-p "$PEDIGREE" --hall "$EX_DIR/pilot/hall2021_MOESM1.txt" --pilot "$EX_DIR/pilot" --pcs "$NGSPCA_DIR/svd.pcs.txt" -o "$OUT" -j "${SLURM_CPUS_PER_TASK:-8}" --cache "${REPORT_CACHE:-$WORK_DIR/report_cache}")
[ -d "$WORK_DIR/counts_scan" ] && [ -n "$(ls "$WORK_DIR"/counts_scan/*.json.gz 2>/dev/null)" ] && args+=(--scan "$WORK_DIR/counts_scan")
[ -d "$WORK_DIR/counts_fetch" ] && [ -n "$(ls "$WORK_DIR"/counts_fetch/*.json.gz 2>/dev/null)" ] && args+=(--fetch "$WORK_DIR/counts_fetch")
[ -d "$WORK_DIR/hprc_censat" ] && args+=(--censat "$WORK_DIR/hprc_censat")     # from 04_hprc_satellites.sh, if it has run
[ -s "${NGSPCA_QC:-}" ] && args+=(--qc "$NGSPCA_QC")                          # NGS-PCA's sample_qc.tsv for the cohort, if at hand
ngsdose_py ngsdose report "${args[@]}"
echo "report: $OUT/index.html"
