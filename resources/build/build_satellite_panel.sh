#!/usr/bin/env bash
# EXPERIMENTAL compositional panels for dispersed sequence - what only a whole-file scan can measure,
# so what a cohort that is scanned once should carry. From the T2T-CHM13v2.0 CenSat annotation:
#   HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR (active and inactive higher-order repeats), and
#   four smaller families of the acrocentric short arms and pericentromeres: ACRO (ACRO1
#   composites), SST1, CER, SATR.
# A k-mer is kept if it occurs >= 10 times in the class's CHM13 arrays, in no other class, and
# nowhere in CHM13 outside CenSat-annotated satellite. One genome's arrays are the only source, so
# recall on other people's arrays has to be measured against assemblies (HPRC:
# example/1000G/04_hprc_satellites.sh). Scan mode only: these reads are spread over centromere
# models, decoys and much else, and no sinks have been learned.
# Not included, and why (share of 150-bp reads from the family's own CHM13 arrays that carry the four
# k-mers a read needs): gamma satellite 13%, divergent alpha HORs 10%, HSat4 0% - too few recurring
# k-mers to measure; monomeric alpha 45%, and as a class of its own it takes 17% of aSatHOR's k-mers
# with it, because a k-mer shared between classes is dropped from both.
# Second panel: TEL, the six canonical 31-mers of (TTAGGG)n, unfiltered (every genome has them at
# its chromosome ends, which is the point). A read is assigned with as few as four of them, i.e.
# 34 bp of perfect repeat, which interstitial telomeric sequence also has; the engine records the
# fraction of each read's k-mers that hit (hit_frac), so that a TelSeq-style threshold can be
# applied afterwards.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"
CHM13="$WORK/chm13v2.0.fa"; CENSAT="$WORK/chm13v2.0_censat_v2.1.bed"
OUT="${OUT:-$ROOT/resources/experimental}"; T="$WORK/satellite_panel"; mkdir -p "$OUT" "$T"
[ -s "$CENSAT" ] || curl -sSL -o "$CENSAT" https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed

: > "$T/classes.tsv"
# CenSat name -> class: by the category before "_" / "(", or, inside the catch-all "censat"
# category, by the family named first in the parentheses (chrY composites that merely contain an
# SST1 subunit are named after their RBMY/TSPY/DAZ subunits and are left out)
family() {  # class-name  awk-condition on (cat, fam)
  awk -v OFS="" '{n=$4; cat=n; sub(/[_(].*/,"",cat); fam=""; if (n ~ /\(/) {fam=n; sub(/^[^(]*\(/,"",fam); sub(/[,)].*/,"",fam)}
                  if (('"$2"') && $3-$2>=2000) print $1,":",$2+1,"-",$3}' "$CENSAT" > "$T/$1.regions.txt"
  [ -s "$T/$1.regions.txt" ] || { echo "no CenSat regions for $1" >&2; exit 1; }
  samtools faidx "$CHM13" -r "$T/$1.regions.txt" > "$T/$1.fa"
  printf '%s\tcompositional\t%s\t0\t10\n' "$1" "$T/$1.fa" >> "$T/classes.tsv"
  printf '%-9s %4d arrays %8.2f Mb\n' "$1" "$(wc -l < "$T/$1.regions.txt")" "$(awk -F'[:-]' '{s+=$3-$2+1} END{print s/1e6}' "$T/$1.regions.txt")"
}
family HSat1A   'cat=="hsat1A"'
family HSat1B   'cat=="hsat1B"'
family HSat2    'cat=="hsat2"'
family HSat3    'cat=="hsat3"'
family bSat     'cat=="bsat"'
family aSatHOR  'cat=="hor"'
family SATR     'cat=="censat" && fam ~ /^SATR/'
family ACRO     'cat=="censat" && fam ~ /^ACRO/'
family SST1     'cat=="censat" && fam ~ /^SST1/'
family CER      'cat=="censat" && fam == "CER"'
# background = CHM13 minus everything CenSat calls satellite ("ct" transition regions stay in the background)
awk 'BEGIN{OFS="\t"} {n=$4; sub(/[_(].*/,"",n); if (n!="ct") print $1,$2,$3,n}' "$CENSAT" > "$T/censat_mask.bed"
"$BIN" panel -m "$T/classes.tsv" -k 31 -b "$CHM13:$T/censat_mask.bed" -o "$OUT/satellites.CHM13v2.k31.panel.tsv.gz"

# telomeric repeat: no background to filter against
python3 -c "print('>TTAGGG_x200'); print('TTAGGG' * 200)" > "$T/TEL.fa"
printf 'TEL\tcompositional\t%s\t0\t10\n' "$T/TEL.fa" > "$T/tel.classes.tsv"
"$BIN" panel -m "$T/tel.classes.tsv" -k 31 -o "$OUT/telomere.k31.panel.tsv.gz"
