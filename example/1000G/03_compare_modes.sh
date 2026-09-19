#!/usr/bin/env bash
# After MODE=scan counts exist (the HPRC subset, or everything): a scan records where every class
# read was aligned, so it says how much of each class the shipped sinks capture in each sample -
# that is, what fetch mode retrieves - and whether new sinks appear with more people.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
OUT="$WORK_DIR/sinks_check"; mkdir -p "$OUT"
ngsdose_py ngsdose sinks "$WORK_DIR"/counts_scan/*.json.gz --evaluate "$BUNDLE/sinks.bed" > "$OUT/capture_per_sample.tsv"
awk -F'\t' 'NR>1 && ($2=="rDNA45S" || $2=="rDNA5S" || $2=="DJ") {n[$2]++; s[$2]+=$5; if (!($2 in m) || $5<m[$2]) m[$2]=$5} END{for (c in n) printf "%s\tsamples=%d\tmean_capture=%.5f\tmin_capture=%.5f\n", c, n[c], s[c]/n[c], m[c]}' "$OUT/capture_per_sample.tsv"
# sinks re-learned from the whole cohort, for comparison with the shipped table
ngsdose_py ngsdose sinks "$WORK_DIR"/counts_scan/*.json.gz -o "$OUT/sinks.relearned.bed"
