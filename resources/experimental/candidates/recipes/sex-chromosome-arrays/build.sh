#!/usr/bin/env bash
# Candidate panel: sex-chromosome ampliconic and tandem arrays (13 classes).
#
#   Y: TSPY (Yp11.2 array, 20,308-bp unit), RBMY (RBMY1 gene), DAZ (DAZ1 gene), BPY2, CDY1, CDY2
#      (positional); DYZ19 (Yq11.221 125-bp satellite, compositional, from the CHM13 array)
#   X: OPN1 (opsin array unit, exon 5 +-300 bp set to N), OPN1LW and OPN1MW (exon 5 +-300 bp of each
#      gene type: after the shared k-mers drop, only the LW/MW-diagnostic ones are left; each keeps
#      only k-mers found in every copy of its type in both assemblies, tools/core_keep.py), GAGE
#      (9,556-bp array unit), CT45 (17,252-bp array unit), SPANXB (SPANXB1 +-8.2 kb; only the k-mers
#      absent from SPANXA/C/D survive)
#
# The circular units (TSPY, OPN1, GAGE, CT45) are exact periods: each ends at the base before the next
# array unit starts (where the unit's first 31-mer recurs), so every 31-mer across the circular join
# occurs in the array.
#
# Every unit is cut from the local GRCh38 analysis set (or, for DYZ19, CHM13v2.0) with samtools faidx;
# the coordinates are below. A k-mer is kept only if it occurs nowhere in GRCh38 (primary, alts, decoys,
# HLA, EBV) or CHM13v2.0 outside the class loci in masks/, and in no shipped panel (panel.k31,
# satellites, telomere): the shipped k-mers are passed as a third background, so a k-mer shared with a
# shipped class cannot survive (the running cohort's classes would otherwise lose it at load time).
#
# The masks were made with tools/loci.py (1-kb chunks of each unit, minimap2 asm20, against GRCh38
# chrX + chrY + all 3,341 non-primary contigs, and CHM13 chrX + chrY; clusters holding >= 25% of the
# unit's chunks at >= 90% identity, padded 1 kb). DERIVE_LOCI=1 re-derives them into $T/loci and
# diffs them against masks/ (needs minimap2; ~5 min, ~4 GB). By hand, on top of the tool's output:
#   SPANXB  only its own locus (SPANXA1/A2/C/D are paralogs at 96-98%, not copies: their k-mers must go)
#   DYZ19   GRCh38 chrY:20053678-20352049 (both sides of the 50-kb gap), chrY_KI270740v1_random and the
#           eleven hs38d1 decoys that are DYZ19 (89-100% of their 31-mers are in the CHM13 array)
#           (a twelfth, chrUn_JTFH01001423v1_decoy, is DYZ19-like at 89% identity, but only 513 of its
#           1,406 distinct 31-mers (36%) are in the array: it stays background, and the class k-mers it
#           holds go)
#   CDY1/2  both families' loci (they are 98.8% identical; the shared k-mers drop as shared)
#
# One mask for all classes lets through a k-mer that one class shares with a paralog copy inside ANOTHER
# class's masked locus: the site at OPN1LW unit 778-808 also occurs in OPN1MW2/MW3 (inside the OPN1 array
# mask, but not in the OPN1MW unit), and OPN1LW then counted MW2/MW3 reads; CDY1 unit 2720-2750 also
# occurs at both CDY2 loci, just outside the CDY2 unit. So CDY1, CDY2, OPN1LW and OPN1MW are each taken
# from a build whose masks hold only that class's own copies (masks/*.own_copies.bed), where every other
# copy is background (tools/merge_own.py), as the macrosatellite recipe does for the D4Z4 classes.
#
# Needs: ngs-dose, samtools, python3 (minimap2 for DERIVE_LOCI=1). Inputs: the GRCh38 analysis set and
# CHM13v2.0 in WORK, the shipped panels. Writes OUT/panel.tsv.gz, OUT/units/, OUT/rep.tsv (the union-mask
# build; T/own.<class>.rep.tsv.gz for the four own-copy builds) and the shared-k-mer checks; RECALL=1 also
# OUT/recall.tsv and OUT/usable.tsv. ~6 min, ~2.2 GB (two own-copy builds at a time).
# Called by resources/build/build_candidate_panels.sh.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/sex-chromosome-arrays}"; T="${T:-$OUT/work}"
G="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"; C="$WORK/chm13v2.0.fa"
SHIPPED=("$ROOT/resources/GRCh38/panel.k31.tsv.gz" "$ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz" "$ROOT/resources/experimental/telomere.k31.panel.tsv.gz")
mkdir -p "$OUT/units" "$T"

# ---- 1. units ---------------------------------------------------------------------------------
# unit NAME FASTA REGION LENGTH DESCRIPTION
unit() {
  samtools faidx "$2" "$3" | sed "1s|.*|>$1 $5|" > "$OUT/units/$1.fa"
  local n; n=$(grep -v '>' "$OUT/units/$1.fa" | tr -d '\n' | wc -c | tr -d ' ')
  [ "$n" -eq "$4" ] || { echo "error: $1 is $n bp, expected $4" >&2; exit 1; }
}
unit TSPY   "$G" chrY:9462691-9482998     20308 "GRCh38 chrY:9462691-9482998, one 20,308-bp TSPY array unit (holds TSPY1): the next unit starts at chrY:9,482,999"
unit RBMY   "$G" chrY:21534879-21549326   14448 "GRCh38 chrY:21534879-21549326, RBMY1A1 gene (RefSeq txStart-txEnd)"
unit DAZ    "$G" chrY:23129355-23199094   69740 "GRCh38 chrY:23129355-23199094, DAZ1 gene (RefSeq txStart-txEnd)"
unit BPY2   "$G" chrY:22984263-23005465   21203 "GRCh38 chrY:22984263-23005465, BPY2 gene (RefSeq txStart-txEnd)"
unit CDY1   "$G" chrY:25622117-25624902    2786 "GRCh38 chrY:25622117-25624902, CDY1 gene (RefSeq txStart-txEnd)"
unit CDY2   "$G" chrY:18025789-18028561    2773 "GRCh38 chrY:18025789-18028561, CDY2A locus: the span homologous to the CDY1 unit"
unit GAGE   "$G" chrX:49532177-49541732    9556 "GRCh38 chrX:49532177-49541732, one GAGE array unit (GAGE12C start to GAGE12D start)"
unit CT45   "$G" chrX:135829229-135846480 17252 "GRCh38 chrX:135829229-135846480, one CT45 array unit (CT45A7 start to the next unit, which starts at chrX:135,846,481)"
unit SPANXB "$G" chrX:140994932-141011369 16438 "GRCh38 chrX:140994932-141011369, SPANXB1 +-8.2 kb (16,438 bp, the SPANXB1 CNV unit length of Lucotte 2018)"
unit OPN1LW "$G" chrX:154155994-154156833   840 "GRCh38 chrX:154155994-154156833, OPN1LW exon 5 +-300 bp"
unit OPN1MW "$G" chrX:154193108-154193947   840 "GRCh38 chrX:154193108-154193947, OPN1MW exon 5 +-300 bp"
unit DYZ19  "$C" chrY:20961204-21226263  265060 "CHM13v2.0 chrY:20961204-21226263, CenSat v2.1 censat_Y_59(DYZ19_125bp_LTR12Bfrag)"
samtools faidx "$G" chrX:154182596-154219733 > "$T/OPN1MWunit.fa"
# the opsin unit: OPN1MW gene start to the next unit (37,138 bp; the unit's first 31-mer recurs at
# chrX:154,219,734, 21 bp before the OPN1MW2 gene start), with OPN1MW exon 5 +-300 bp
# (unit 10513-11352, the OPN1MW class) set to N so that the exon-5 reads go to OPN1LW / OPN1MW
"$PY" - "$T/OPN1MWunit.fa" "$OUT/units/OPN1.fa" <<'PYEOF'
import sys
s = ''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>'))
assert len(s) == 37138, len(s)
a, b = 154193107 - 154182595, 154193947 - 154182595
t = s[:a] + 'N' * (b - a) + s[b:]
with open(sys.argv[2], 'w') as fh:
    fh.write('>OPN1 GRCh38 chrX:154182596-154219733 (OPN1MW gene start to the next unit), unit 10513-11352 (OPN1MW exon 5 +-300 bp) set to N\n')
    fh.write('\n'.join(t[i:i + 60] for i in range(0, len(t), 60)) + '\n')
PYEOF

# ---- 2. class loci (the background masks) -------------------------------------------------------
if [ "${DERIVE_LOCI:-0}" = 1 ]; then
  mkdir -p "$T/loci"; ( cd "$T/loci"
  awk '$1!~/^chr([0-9]+|M)$/{print $1}' "$G.fai" > g38.list
  [ -s g38.xy_np.fa ] || samtools faidx "$G" -r g38.list > g38.xy_np.fa
  [ -s t2t.xy.fa ] || samtools faidx "$C" chrX chrY > t2t.xy.fa
  cp "$T/OPN1MWunit.fa" OPN1.src.fa
  for c in TSPY RBMY DAZ BPY2 GAGE CT45 SPANXB OPN1 CDY1 CDY2 DYZ19; do
    src="$OUT/units/$c.fa"; [ "$c" = OPN1 ] && src=OPN1.src.fa
    opt=(); case $c in CDY1|CDY2) opt=(--size 500 --step 250);; DYZ19) opt=(--size 2000 --step 2000 --min-frac 0.02 --gap 20000);; esac
    "$PY" "$HERE/tools/loci.py" "$src" g38.xy_np.fa "$c" "${opt[@]}" > "g38.$c.bed" 2> "g38.$c.log"
    "$PY" "$HERE/tools/loci.py" "$src" t2t.xy.fa "$c" "${opt[@]}" > "t2t.$c.bed" 2> "t2t.$c.log"
  done )
  echo "loci re-derived in $T/loci; compare with masks/ (SPANXB: first locus only; DYZ19: see the header)" >&2
fi
MASK_G="$HERE/masks/GRCh38.sexarr_loci.bed"; MASK_C="$HERE/masks/CHM13v2.sexarr_loci.bed"
OWN_G="$HERE/masks/GRCh38.own_copies.bed"; OWN_C="$HERE/masks/CHM13v2.own_copies.bed"
own() {  # own CLASS BED: that class's rows of an own-copies BED
  awk -v c="$1" 'BEGIN {OFS = "\t"} !/^#/ && $4 == c' "$2"
}

# ---- 2b. OPN1LW / OPN1MW keep BEDs: unit positions whose k-mer occurs once in every copy of the type
# (1 LW + 3 MW windows in GRCh38, 1 + 2 in CHM13), so each copy of a type gives the class the same k-mers
mkdir -p "$OUT/keep"
for c in OPN1LW OPN1MW; do
  { own $c "$OWN_G" | awk '{print $1 ":" $2 + 1 "-" $3}' | xargs samtools faidx "$G"
    own $c "$OWN_C" | awk '{print $1 ":" $2 + 1 "-" $3}' | xargs samtools faidx "$C" | sed 's/^>/>CHM13_/'; } > "$T/$c.copies.fa"
  "$PY" "$HERE/tools/core_keep.py" "$OUT/units/$c.fa" "$T/$c.copies.fa" > "$OUT/keep/$c.bed"
done

# ---- 3. the shipped panels as a background --------------------------------------------------------
"$PY" "$HERE/tools/shared.py" --fasta "$T/shipped_kmers.fa" "${SHIPPED[@]}"

# ---- 4. the panel -----------------------------------------------------------------------------------
cat > "$OUT/manifest.tsv" <<'EOF'
TSPY	positional	units/TSPY.fa	1
RBMY	positional	units/RBMY.fa	0
DAZ	positional	units/DAZ.fa	0
BPY2	positional	units/BPY2.fa	0
CDY1	positional	units/CDY1.fa	0
CDY2	positional	units/CDY2.fa	0
DYZ19	compositional	units/DYZ19.fa	0	10
OPN1	positional	units/OPN1.fa	1
OPN1LW	positional	units/OPN1LW.fa	0	.	keep/OPN1LW.bed
OPN1MW	positional	units/OPN1MW.fa	0	.	keep/OPN1MW.bed
GAGE	positional	units/GAGE.fa	1
CT45	positional	units/CT45.fa	1
SPANXB	positional	units/SPANXB.fa	0
EOF
# without the shipped background first, to count what the shipped panels would remove
"$BIN" panel -m "$OUT/manifest.tsv" -k 31 -b "$G:$MASK_G" -b "$C:$MASK_C" --max-bg 0 -o "$T/noshipped.panel.tsv.gz" 2> "$T/noshipped.log"
"$PY" "$HERE/tools/shared.py" "$T/noshipped.panel.tsv.gz" "${SHIPPED[@]}" > "$OUT/shared_with_shipped.tsv"
"$BIN" panel -m "$OUT/manifest.tsv" -k 31 -b "$G:$MASK_G" -b "$C:$MASK_C" -b "$T/shipped_kmers.fa" --max-bg 0 \
  -o "$T/union.panel.tsv.gz" --report "$T/rep.tsv.gz" 2> "$OUT/panel.log"
gzip -dc "$T/rep.tsv.gz" > "$OUT/rep.tsv"
grep -v scanning "$OUT/panel.log"
# 4b. the paralog-pair classes, each from a build that masks only its own copies (same manifest, so
# the class ids agree and the within-panel shared k-mers are the same)
ownbuild() {
  own "$1" "$OWN_G" > "$T/own.$1.GRCh38.bed"; own "$1" "$OWN_C" > "$T/own.$1.CHM13.bed"
  "$BIN" panel -m "$OUT/manifest.tsv" -k 31 -b "$G:$T/own.$1.GRCh38.bed" -b "$C:$T/own.$1.CHM13.bed" -b "$T/shipped_kmers.fa" \
    --max-bg 0 -o "$T/own.$1.panel.tsv.gz" --report "$T/own.$1.rep.tsv.gz" 2> "$T/own.$1.log"
  awk -v c="$1" '$1 == c' "$T/own.$1.log"
}
ownbuild CDY1 & p1=$!; ownbuild CDY2 & p2=$!; wait $p1; wait $p2
ownbuild OPN1LW & p1=$!; ownbuild OPN1MW & p2=$!; wait $p1; wait $p2
"$PY" "$HERE/tools/merge_own.py" "$T/union.panel.tsv.gz" "$OUT/panel.tsv.gz" \
  CDY1="$T/own.CDY1.panel.tsv.gz" CDY2="$T/own.CDY2.panel.tsv.gz" OPN1LW="$T/own.OPN1LW.panel.tsv.gz" OPN1MW="$T/own.OPN1MW.panel.tsv.gz"

# ---- 5. checks: shared k-mers (must be 0), read recall per source copy, usable fraction ---------------
"$PY" "$HERE/tools/shared.py" "$OUT/panel.tsv.gz" "${SHIPPED[@]}" > "$OUT/shared_final.tsv"
awk -F'\t' '$3 != 0 {bad = 1} END {exit bad}' "$OUT/shared_final.tsv" || { echo "error: the panel shares k-mers with a shipped panel" >&2; exit 1; }
[ "${RECALL:-0}" = 1 ] || { echo "panel: $OUT/panel.tsv.gz"; exit 0; }
sed "s|@G@|$G|; s|@C@|$C|" "$HERE/copies.tsv.in" > "$T/copies.tsv"
"$PY" "$HERE/tools/recall_tiled.py" "$OUT/panel.tsv.gz" "$T/copies.tsv" --other "${SHIPPED[@]}" > "$OUT/recall.tsv" 2> "$OUT/recall_summary.txt"
PYTHONPATH="$ROOT" "$PY" "$HERE/tools/usable.py" "$OUT/panel.tsv.gz" > "$OUT/usable.tsv"
echo "panel: $OUT/panel.tsv.gz"
