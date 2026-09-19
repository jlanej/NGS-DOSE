#!/usr/bin/env bash
# Dispatch 01b_dose_sample.sh for every manifest sample whose CRAM and index have landed in
# $CRAM_DIR and whose counts are not there yet. It never downloads and never deletes, so it can
# run beside whatever is staging the files (NGS-PCA's aria2 download manager, Globus, rsync);
# a file that aria2 is still writing (it keeps a .aria2 control file beside it) is left alone.
# Run it again, or from cron, as files arrive.   RUN_LOCAL=1 runs the jobs here instead of sbatch.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
cd "$EX_DIR"; mkdir -p logs "$WORK_DIR/dispatched"
n=0; line=0
while IFS=$'\t' read -r sample _ md5; do
  line=$((line + 1))
  cram="$CRAM_DIR/$sample.cram"
  [ -s "$cram" ] && [ -s "$cram.crai" ] && [ ! -e "$cram.aria2" ] && [ ! -e "$cram.crai.aria2" ] || continue
  [ -s "$WORK_DIR/counts_scan/$sample.json.gz" ] && [ -s "$WORK_DIR/counts_fetch/$sample.json.gz" ] && continue
  if [ "${RUN_LOCAL:-0}" = 1 ]; then
    bash "$EX_DIR/01b_dose_sample.sh" "$sample" "$line" "${md5:-}"
  else
    mark="$WORK_DIR/dispatched/$sample"          # one job per sample; remove the mark to dispatch again
    [ -e "$mark" ] && continue
    sbatch --parsable "$EX_DIR/01b_dose_sample.sh" "$sample" "$line" "${md5:-}" > "$mark"
  fi
  n=$((n + 1))
done < "$MANIFEST"
echo "dispatched $n sample(s)"
