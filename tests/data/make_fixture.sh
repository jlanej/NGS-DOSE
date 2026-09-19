#!/usr/bin/env bash
# Rebuild tests/data/NA12878.subsample.bam from the full NA12878 CRAM of the 1000 Genomes 30x
# resource (ERR3239334). Needed only when the bundle's control regions or sinks change.
#
#   tests/data/make_fixture.sh NA12878.final.cram GRCh38_full_analysis_set_plus_decoy_hla.fa
#
# A 2% template subsample (samtools -s 1.02: seed 1, so the same templates every time) of the
# reads placed in the bundle's control regions (with their flanks and the engine's fetch
# padding) and class sinks; read names replaced by serial numbers, qualities and tags dropped,
# contigs without reads dropped from the header.
set -euo pipefail
cram="$1"; ref="$2"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bundle="${BUNDLE:-$here/../../resources/GRCh38}"
out="${OUT:-$here/NA12878.subsample.bam}"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

{ awk 'BEGIN{OFS="\t"} {s=$2-1600; if (s<0) s=0; print $1, s, $3+1600}' "$bundle/controls.bed"; cut -f1-3 "$bundle/sinks.bed"; } \
  | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} $1==c && $2<=e {if ($3>e) e=$3; next} {if (c!="") print c, s, e; c=$1; s=$2; e=$3} END{print c, s, e}' > "$tmp/regions.bed"
samtools view -T "$ref" -s 1.02 -M -L "$tmp/regions.bed" "$cram" \
  | awk 'BEGIN{OFS="\t"} {if (!($1 in id)) id[$1]="r" n++; print id[$1], $2, $3, $4, $5, $6, $7, $8, $9, $10, "*", "RG:Z:N"}' > "$tmp/body.sam"
{ printf '@HD\tVN:1.6\tSO:coordinate\n'
  samtools view -H "$cram" | awk -F'\t' 'NR==FNR {keep[$3]=1; if ($7!="=") keep[$7]=1; next} $1=="@SQ" && (substr($2, 4) in keep) {print $1 "\t" $2 "\t" $3}' "$tmp/body.sam" -
  printf '@RG\tID:N\tSM:NA12878\n'
  cat "$tmp/body.sam"; } | samtools view -b -o "$out" -
samtools index -c "$out"
echo "$(samtools view -c "$out") reads -> $out"
