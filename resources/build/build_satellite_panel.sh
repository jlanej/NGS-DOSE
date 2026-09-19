#!/usr/bin/env bash
# EXPERIMENTAL compositional panel for satellite families, from the T2T-CHM13v2.0 CenSat annotation.
# Classes: HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR. A k-mer is kept if it occurs >= 10 times
# in the class's CHM13 arrays, in no other class, and nowhere in CHM13 outside CenSat-annotated
# satellite. One genome's arrays are the only source, so recall on other people's arrays is
# unknown until it has been measured against assemblies (HPRC). Scan mode only: satellite reads
# are spread over centromere models, decoys and much else, and no sinks have been learned.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"
CHM13="$WORK/chm13v2.0.fa"; CENSAT="$WORK/chm13v2.0_censat_v2.1.bed"
OUT="${OUT:-$ROOT/resources/experimental}"; T="$WORK/satellite_panel"; mkdir -p "$OUT" "$T"
[ -s "$CENSAT" ] || curl -sSL -o "$CENSAT" https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed

: > "$T/classes.tsv"
for pair in hsat1A:HSat1A hsat1B:HSat1B hsat2:HSat2 hsat3:HSat3 bsat:bSat hor:aSatHOR; do
  key="${pair%%:*}"; name="${pair##*:}"
  awk -v c="$key" 'BEGIN{OFS=""} {n=$4; sub(/[_(].*/,"",n); if (n==c && $3-$2>=2000) print $1,":",$2+1,"-",$3}' "$CENSAT" > "$T/$name.regions.txt"
  samtools faidx "$CHM13" -r "$T/$name.regions.txt" > "$T/$name.fa"
  printf '%s\tcompositional\t%s\t0\t10\n' "$name" "$T/$name.fa" >> "$T/classes.tsv"
done
# background = CHM13 minus everything CenSat calls satellite ("ct" transition regions stay in the background)
awk 'BEGIN{OFS="\t"} {n=$4; sub(/[_(].*/,"",n); if (n!="ct") print $1,$2,$3,n}' "$CENSAT" > "$T/censat_mask.bed"
"$BIN" panel -m "$T/classes.tsv" -k 31 -b "$CHM13:$T/censat_mask.bed" -o "$OUT/satellites.CHM13v2.k31.panel.tsv.gz"
