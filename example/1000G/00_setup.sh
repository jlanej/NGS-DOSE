#!/usr/bin/env bash
# One-time setup: directories, reference, pedigree, and the sample -> CRAM URL manifest.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
mkdir -p "$WORK_DIR/reference" "$COUNTS_DIR" "$EST_DIR" "$LOG_DIR" "$WORK_DIR/crai"
mkdir -p "$EX_DIR/logs"      # SLURM opens `--output=logs/...` relative to the submit directory before the job starts

# the container, when one is asked for (SIF=/path/to/ngs-dose.sif): pull the published image
if [ -n "$SIF" ] && [ ! -s "$SIF" ]; then
  apptainer pull "$SIF" "${NGSDOSE_IMAGE:-docker://ghcr.io/jlanej/ngs-dose:latest}" || {
    echo "ERROR: could not pull the image. Build it instead, from the repository root:" >&2
    echo "  apptainer build --fakeroot $SIF example/1000G/ngs-dose.def" >&2; exit 1; }
fi

if [ ! -s "$REF_FASTA" ]; then
  curl -sSL -o "$REF_FASTA" https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  curl -sSL -o "$REF_FASTA.fai" https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
fi
[ -s "$PEDIGREE" ] || curl -sSL -o "$PEDIGREE" "$PEDIGREE_URL"

# manifest: SAMPLE <TAB> CRAM_URL, from the two NYGC sequence indexes (2,504 + 698 samples)
: > "$MANIFEST"
for pair in "1000G_2504_high_coverage.sequence.index:data" \
            "additional_698_related/1000G_698_related_high_coverage.sequence.index:additional_698_related/data"; do
  idx="${pair%%:*}"; sub="${pair##*:}"
  # the file name comes from the index's own path column, the sample name is trimmed: three
  # sample names in the 2504 index carry a trailing space
  curl -sSL "$S3_HTTPS_BASE/$idx" | awk -F'\t' -v b="$S3_HTTPS_BASE/$sub" '!/^#/ {n=split($1, p, "/"); s=$10; gsub(/[ \t\r]+/, "", s); print s "\t" b "/" $3 "/" p[n]}' >> "$MANIFEST"
done
sort -u -o "$MANIFEST" "$MANIFEST"
echo "manifest: $(wc -l < "$MANIFEST") samples -> $MANIFEST"
[ "$(wc -l < "$MANIFEST")" -ge 3000 ] || { echo "ERROR: expected 3,202 samples in the manifest; is the network reachable from this node?" >&2; exit 1; }
if awk -F'\t' 'NF != 2 || $0 ~ / / || $2 !~ /\.cram$/ || seen[$1]++ {bad = 1} END {exit !bad}' "$MANIFEST"; then
  echo "ERROR: malformed manifest line (whitespace, duplicate sample, or not a .cram)" >&2; exit 1
fi

# the samples that also have an HPRC release-2 assembly: the whole-file scan subset
awk -F/ '{print $3}' "$EX_DIR/hprc_r2_censat.keys.txt" | sort -u | awk 'NR==FNR {h[$1]=1; next} $1 in h' - "$MANIFEST" > "$WORK_DIR/manifest.hprc.tsv"
echo "manifest: $(wc -l < "$WORK_DIR/manifest.hprc.tsv") samples with HPRC assemblies -> $WORK_DIR/manifest.hprc.tsv"
