#!/usr/bin/env bash
# Candidate panel "rna-arrays": five positional classes of RNA-gene arrays, for scans first
# (where their reads land is learned from whole-file scans, as for every class; nothing here
# assumes the reads sit at the reference coordinates).
#   TDNA1Q23  the 1q23.3 tRNA-gene macrosatellite: 7,380-bp tandem unit with five tRNA genes
#   RNU1      canonical U1 snRNA gene units at 1p36.13 (not the variant vU1 genes at 1q12-q21)
#   SNORD3    the U3 snoRNA module (SNORD3A/B/C/D) of the chr17:19.02-19.24 Mb REPA/REPB duplicons
#   SNORD116  the SNORD116 cluster at 15q11.2 (30 genes; the whole diverged array is the unit)
#   SNORD115  the SNORD115 cluster at 15q11.2 (48 genes; likewise)
# Units come from the local GRCh38 analysis set and CHM13v2.0 (samtools faidx), coordinates in
# copies/copies.tsv (checked against a fresh derivation, scripts/locate.py). The panel is built
# with the engine's panel builder against both assemblies with every copy of every class masked,
# then filtered (scripts/post_filter.py) so that it shares no k-mer with the shipped panels and no
# class keeps a k-mer found at another class's loci, and finally (scripts/neighbors.py) so that no
# k-mer is one substitution away from any sequence outside the class loci in either assembly.
#
# Needs: samtools, minimap2, python3, the ngs-dose binary; the two assemblies (WORK); the shipped
# panels. Writes OUT/panel.tsv.gz, OUT/units/, OUT/masks/, OUT/rep.tsv; RECALL=1 also OUT/recall.tsv.
# ~7 min, ~3.2 GB memory (the builder streams both genomes twice). Called by
# resources/build/build_candidate_panels.sh.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/rna-arrays}"
G38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"
C13="$WORK/chm13v2.0.fa"
SHIPPED=("$ROOT/resources/GRCh38/panel.k31.tsv.gz"
         "$ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz"
         "$ROOT/resources/experimental/telomere.k31.panel.tsv.gz")
TMP="$OUT/tmp"; mkdir -p "$TMP" "$OUT/units" "$OUT/masks" "$OUT/copies"
for f in "$G38" "$G38.fai" "$C13" "$C13.fai" "${SHIPPED[@]}"; do [ -s "$f" ] || { echo "error: $f missing" >&2; exit 1; }; done
CLASSES="TDNA1Q23 RNU1 SNORD3 SNORD116 SNORD115"

# ---- 1. units ------------------------------------------------------------------------------
# faidx_to OUT NAME [-i] REGION: one FASTA record named NAME
faidx_to() { local out=$1 name=$2; shift 2; samtools faidx "$@" | sed "1s/.*/>$name/" > "$out"; }
faidx_to "$OUT/units/TDNA1Q23.fa" TDNA1Q23_GRCh38_chr1_161447001_161454380 "$G38" chr1:161447001-161454380   # starts 227 bp before tRNA-Glu-CTC-1-2
faidx_to "$OUT/units/RNU1.fa"     RNU1_CHM13_chr1_16118190_16130189_rc   "$C13" -i chr1:16118190-16130189  # gene -3 kb .. +9 kb, gene on +
faidx_to "$OUT/units/SNORD3.fa"   SNORD3_GRCh38_chr17_19111334_19115875  "$G38" chr17:19111334-19115875   # U3 module with SNORD3D
faidx_to "$OUT/units/SNORD116.fa" SNORD116_GRCh38_chr15_25051001_25109000 "$G38" chr15:25051001-25109000  # SNORD116-1 .. -30
faidx_to "$OUT/units/SNORD115.fa" SNORD115_GRCh38_chr15_25170001_25270500 "$G38" chr15:25170001-25270500  # SNORD115-1 .. -48
for c in $CLASSES; do grep -q N <(tail -n +2 "$OUT/units/$c.fa") && { echo "error: unit $c contains N" >&2; exit 1; }; done
printf '%s\n' 'unit lengths:'; for c in $CLASSES; do printf '  %-9s %7d\n' "$c" "$(awk '!/^>/{n+=length($0)} END{print n}' "$OUT/units/$c.fa")"; done

# ---- 2. copies: check the table against a fresh derivation, then extract copies and masks ----
"$PY" "$HERE/scripts/locate.py" "$G38" "$C13" "$OUT/units" "$TMP/locate" > "$TMP/derived.tsv"
grep -v '^#' "$HERE/copies/copies.tsv" | cut -f1-6 > "$TMP/table.tsv"
if ! diff <(sort "$TMP/table.tsv") <(sort "$TMP/derived.tsv") > "$TMP/copies.diff"; then
  echo "error: copies/copies.tsv differs from the derivation (scripts/locate.py):" >&2; cat "$TMP/copies.diff" >&2; exit 1
fi
echo "copies/copies.tsv agrees with the derivation ($(wc -l < "$TMP/table.tsv") copies)"
: > "$TMP/grch38.bed"; : > "$TMP/chm13.bed"
for c in $CLASSES; do : > "$OUT/copies/$c.fa"; : > "$TMP/$c.loci.fa"; done
grep -v '^#' "$HERE/copies/copies.tsv" | while IFS=$'\t' read -r cls asm chrom s e strand pad what; do
  fa=$G38; [ "$asm" = CHM13 ] && fa=$C13
  inv=(); [ "$strand" = - ] && inv=(-i)
  samtools faidx ${inv[@]+"${inv[@]}"} "$fa" "$chrom:$((s + 1))-$e" | sed "1s/.*/>${cls}_${asm}_${chrom}_${s}_${e}_${strand/-/minus}/" | sed 's/_+$/_plus/' >> "$OUT/copies/$cls.fa"
  ms=$((s - pad)); [ $ms -lt 0 ] && ms=0
  me=$((e + pad)); len=$(awk -v c="$chrom" '$1==c{print $2}' "$fa.fai"); [ $me -gt "$len" ] && me=$len
  # the unplaced U1 contig is a U1 unit end to end: mask all of it
  [ "$chrom" = chr1_KI270713v1_random ] && { ms=0; me=$len; }
  bed=$TMP/grch38.bed; [ "$asm" = CHM13 ] && bed=$TMP/chm13.bed
  printf '%s\t%d\t%d\t%s\n' "$chrom" "$ms" "$me" "$cls" >> "$bed"
  printf '%s\t%d\t%d\t%s\t%s\n' "$chrom" "$ms" "$me" "$cls" "$asm" >> "$TMP/all_masks.tsv.part"
  samtools faidx "$fa" "$chrom:$((ms + 1))-$me" >> "$TMP/$cls.loci.fa"
done
sort -k1,1 -k2,2n "$TMP/grch38.bed" > "$OUT/masks/grch38.rna_arrays.bed"
sort -k1,1 -k2,2n "$TMP/chm13.bed"  > "$OUT/masks/chm13.rna_arrays.bed"
for c in $CLASSES; do
  awk -v c="$c" '$4==c' "$OUT/masks/grch38.rna_arrays.bed" > "$OUT/masks/$c.grch38.bed"
  awk -v c="$c" '$4==c' "$OUT/masks/chm13.rna_arrays.bed"  > "$OUT/masks/$c.chm13.bed"
done
rm -f "$TMP/all_masks.tsv.part"

# ---- 3. RNU1: core k-mers, present once in every one of the 12 reference U1 units ------------
"$PY" "$HERE/scripts/core_keep.py" "$OUT/units/RNU1.fa" "$OUT/copies/RNU1.fa" > "$OUT/masks/RNU1.core_keep.bed"

# ---- 4. panel builder: both assemblies, every class copy masked --------------------------------
cat > "$OUT/manifest.tsv" <<EOF
TDNA1Q23	positional	units/TDNA1Q23.fa	1
RNU1	positional	units/RNU1.fa	0	.	masks/RNU1.core_keep.bed
SNORD3	positional	units/SNORD3.fa	0
SNORD116	positional	units/SNORD116.fa	0
SNORD115	positional	units/SNORD115.fa	0
EOF
"$BIN" panel -m "$OUT/manifest.tsv" -k 31 --max-bg 0 \
  -b "$G38:$OUT/masks/grch38.rna_arrays.bed" \
  -b "$C13:$OUT/masks/chm13.rna_arrays.bed" \
  --report "$OUT/rep.tsv.gz" -o "$TMP/panel.builder.tsv.gz" 2> "$TMP/panel.log"
cat "$TMP/panel.log"

# ---- 5. no k-mer shared with a shipped panel, none found at another class's loci ---------------
loci=(); for c in $CLASSES; do loci+=("$c=$TMP/$c.loci.fa"); done
"$PY" "$HERE/scripts/post_filter.py" "$TMP/panel.builder.tsv.gz" "$TMP/panel.exact.tsv.gz" \
  --shipped "${SHIPPED[@]}" --loci "${loci[@]}" --summary "$TMP/post_filter.tsv"

# ---- 6. one-mismatch background: drop k-mers one substitution from sequence outside the class loci --
# (in the NA12878 scan made with the exact-filtered panel, 0.1-0.4% of each class's reads came from
# single loci elsewhere whose sequence, with the sample's allele, matched 4-7 consecutive class k-mers)
rm -rf "$TMP/nb"
"$PY" "$HERE/scripts/neighbors.py" write "$TMP/panel.exact.tsv.gz" "$TMP/nb"
"$BIN" panel -m "$TMP/nb/manifest.tsv" -k 31 --max-bg 0 \
  -b "$G38:$OUT/masks/grch38.rna_arrays.bed" \
  -b "$C13:$OUT/masks/chm13.rna_arrays.bed" \
  --report "$TMP/nb/report.tsv.gz" -o "$TMP/nb/panel.tsv.gz" 2> "$TMP/nb/build.log"
"$PY" "$HERE/scripts/neighbors.py" filter "$TMP/panel.exact.tsv.gz" "$TMP/nb" "$TMP/nb/report.tsv.gz" \
  "$OUT/panel.tsv.gz" "$TMP/neighbors.tsv"
rm -f "$TMP/nb"/nb*.fa "$TMP/nb/report.tsv.gz" "$TMP/nb/panel.tsv.gz" "$TMP/nb/neighbours.txt"   # ~1.5 GB

"$PY" - "$OUT/rep.tsv.gz" "$TMP/post_filter.tsv" "$TMP/neighbors.tsv" > "$OUT/rep.tsv" <<'EOF'
import gzip, sys
st = {}
order = []
for line in gzip.open(sys.argv[1], "rt"):
    if line.startswith("#"):
        continue
    c, pos, cnt, shared, bg = line.rstrip("\n").split("\t")
    if c not in st:
        st[c] = dict(input=0, multicopy_in_unit=0, shared_within_group=0, background=0); order.append(c)
    s = st[c]; s["input"] += 1
    if shared == "1": s["shared_within_group"] += 1
    elif int(cnt) > 1: s["multicopy_in_unit"] += 1
    elif int(bg) > 0: s["background"] += 1
pf = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(sys.argv[2]).readlines()[1:]}
nb = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(sys.argv[3]).readlines()[1:]}
print("class\tkmers_input\tdropped_multicopy_in_unit\tdropped_shared_within_group\tdropped_background_exact\tkept_by_builder\tremoved_shared_with_shipped\tremoved_at_other_class_loci\tremoved_one_mismatch_background\tfinal")
for c in order:
    s = st[c]; p = pf[c]; n = nb[c]
    print(f"{c}\t{s['input']}\t{s['multicopy_in_unit']}\t{s['shared_within_group']}\t{s['background']}\t{p[1]}\t{p[2]}\t{p[3]}\t{n[2]}\t{n[3]}")
EOF
cat "$OUT/rep.tsv"

[ "${RECALL:-0}" = 1 ] || exit 0
# ---- 7. recall per source copy (150-bp reads every 10 bp) ---------------------------------------
specs=(); for c in $CLASSES; do specs+=("$c=$OUT/copies/$c.fa"); done
"$PY" "$HERE/scripts/recall.py" "$OUT/panel.tsv.gz" "${specs[@]}" > "$OUT/recall.tsv"
cat "$OUT/recall.tsv"
