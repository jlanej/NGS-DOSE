#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_stage
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=96:00:00
#SBATCH --output=logs/stage_%j.out
# Stage the cohort a few CRAMs at a time and hand each one to 01b_dose_sample.sh (whole-file scan +
# targeted fetch, then the CRAM is removed) the moment it has landed and its MD5 has been checked.
# The WAN is the scarce, failure-prone resource, so it gets a few gentle multi-connection
# transfers; counting is ordinary compute and gets whatever the scheduler grants. At most
# MAX_LOCAL_CRAMS files are on disk at any time (~16 GB each).
#
#   sbatch 01_stage_and_dose.sh            # or, on a login / transfer node, in tmux:  bash 01_stage_and_dose.sh
#   LIMIT=5 bash 01_stage_and_dose.sh      # the first five samples only: a smoke test
#   RUN_LOCAL=1 ...                        # count here instead of submitting jobs
#
# Idempotent: samples with both counts files are skipped; a sample whose job is still queued or
# running is left alone; a CRAM that is already on disk is verified, not downloaded again; aria2c
# resumes a partial file. aria2c comes from the container image when SIF is set (curl otherwise,
# or when the host has no aria2c).
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
DOWNLOAD_SLOTS="${DOWNLOAD_SLOTS:-4}"                # concurrent transfers
DOWNLOAD_CONNECTIONS="${DOWNLOAD_CONNECTIONS:-8}"    # connections per transfer
MAX_LOCAL_CRAMS="${MAX_LOCAL_CRAMS:-40}"
LIMIT="${LIMIT:-0}"
cd "$EX_DIR"; mkdir -p logs "$CRAM_DIR" "$WORK_DIR/dispatched" "$WORK_DIR/counts_scan" "$WORK_DIR/counts_fetch"

have_aria2() { ngsdose_py aria2c --version >/dev/null 2>&1; }
md5_of() { (md5sum "$1" 2>/dev/null || md5 -r "$1") | awk '{print $1}'; }
done_sample() { [ -s "$WORK_DIR/counts_scan/$1.json.gz" ] && [ -s "$WORK_DIR/counts_fetch/$1.json.gz" ]; }
job_alive() { [ -s "$WORK_DIR/dispatched/$1" ] && command -v squeue >/dev/null 2>&1 && squeue -h -j "$(cat "$WORK_DIR/dispatched/$1")" 2>/dev/null | grep -q .; }

download() {  # url destination [md5]
  local url="$1" dest="$2" md5="${3:-}" try
  for try in 1 2 3 4; do
    if have_aria2; then
      ngsdose_py aria2c --quiet=true --console-log-level=warn --summary-interval=0 -x "$DOWNLOAD_CONNECTIONS" -s "$DOWNLOAD_CONNECTIONS" -k 16M \
        --file-allocation=none --auto-file-renaming=false --allow-overwrite=true --continue=true --max-tries=5 --retry-wait=20 \
        ${md5:+--checksum=md5=$md5} -d "$(dirname "$dest")" -o "$(basename "$dest")" "$url" && return 0
    else
      curl -sSL --fail --retry 5 -C - -o "$dest" "$url" && { [ -z "$md5" ] || [ "$(md5_of "$dest")" = "$md5" ]; } && return 0
      [ -z "$md5" ] || rm -f "$dest"                  # a complete file with the wrong MD5 is not worth resuming
    fi
    echo "  download attempt $try failed: $(basename "$dest")" >&2; sleep $(( try * 30 ))
  done
  return 1
}

stage_one() {  # sample url md5 line
  local sample="$1" url="$2" md5="$3" line="$4" cram="$CRAM_DIR/$1.cram"
  if [ -s "$cram" ] && [ ! -e "$cram.aria2" ] && [ "$(md5_of "$cram")" = "$md5" ]; then
    echo "[$sample] already on disk, MD5 verified"
  else
    download "$url" "$cram" "$md5" || { echo "[$sample] FAILED to download" >&2; return 1; }
  fi
  [ -s "$cram.crai" ] || download "$url.crai" "$cram.crai" || { echo "[$sample] FAILED to download the index" >&2; return 1; }
  if [ "${RUN_LOCAL:-0}" = 1 ]; then
    bash "$EX_DIR/01b_dose_sample.sh" "$sample" "$line"            # MD5 checked above
  else
    sbatch --parsable "$EX_DIR/01b_dose_sample.sh" "$sample" "$line" > "$WORK_DIR/dispatched/$sample" && echo "[$sample] staged, job $(cat "$WORK_DIR/dispatched/$sample")"
  fi
}

line=0; started=0
while IFS=$'\t' read -r sample url md5; do
  line=$((line + 1))
  done_sample "$sample" && continue
  job_alive "$sample" && continue
  case "$url" in *://*) ;; *) echo "[$sample] the manifest gives a local path, nothing to stage" >&2; continue ;; esac
  [ "$LIMIT" -gt 0 ] && [ "$started" -ge "$LIMIT" ] && break
  # back-pressure: the per-sample jobs delete their CRAM when they are done
  while [ "$(find "$CRAM_DIR" -maxdepth 1 -name '*.cram' | wc -l)" -ge "$MAX_LOCAL_CRAMS" ]; do sleep 60; done
  while [ "$(jobs -rp | wc -l)" -ge "$DOWNLOAD_SLOTS" ]; do sleep 5; done
  stage_one "$sample" "$url" "$md5" "$line" &
  started=$((started + 1))
done < "$MANIFEST"
wait
left=0; while IFS=$'\t' read -r sample _; do done_sample "$sample" || left=$((left + 1)); done < "$MANIFEST"
echo "staging pass finished: $started sample(s) started; $left of $(( $(wc -l < "$MANIFEST") )) without counts yet (queued jobs included). Run again to continue."
