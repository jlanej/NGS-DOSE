#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_hprc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/hprc_%j.out
# Scan-mode satellite estimates against the HPRC release-2 assemblies of the same people.
# Needs MODE=scan counts for the HPRC samples (manifest.hprc.tsv, written by 00_setup.sh):
#   N=$(wc -l < $WORK_DIR/manifest.hprc.tsv)
#   MODE=scan MANIFEST=$WORK_DIR/manifest.hprc.tsv sbatch --array=0-$(( (N - 1) / 10 ))%10 01_count.sh
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
MODE=scan; COUNTS_DIR="$WORK_DIR/counts_scan"
CENSAT="$WORK_DIR/hprc_censat"; mkdir -p "$CENSAT"

# CenSat annotations: two small BED files per sample (open access, HPRC S3 bucket)
while read -r key; do
  f="$CENSAT/$(basename "$key")"
  [ -s "$f" ] || curl -sSL --retry 5 -o "$f" "https://s3-us-west-2.amazonaws.com/human-pangenomics/$key"
done < "$EX_DIR/hprc_r2_censat.keys.txt"

counts=()
while IFS=$'\t' read -r sample _; do
  [ -s "$COUNTS_DIR/$sample.json.gz" ] && counts+=("$COUNTS_DIR/$sample.json.gz")
done < "$WORK_DIR/manifest.hprc.tsv"
[ "${#counts[@]}" -gt 0 ] || { echo "no scan-mode counts for the HPRC samples in $COUNTS_DIR" >&2; exit 1; }
echo "${#counts[@]} HPRC samples with scan-mode counts"
ngsdose_py ngsdose estimate "${counts[@]}" -r "$BUNDLE" -j "${SLURM_CPUS_PER_TASK:-4}" -t "$WORK_DIR/hprc_satellite_estimates.tsv"
ngsdose_py python3 "$EX_DIR/hprc_satellites.py" --estimates "$WORK_DIR/hprc_satellite_estimates.tsv" --censat "$CENSAT" --out "$WORK_DIR/hprc_satellites.tsv"
