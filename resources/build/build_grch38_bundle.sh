#!/usr/bin/env bash
# Rebuild the reference-derived part of the GRCh38 bundle from public inputs: units/, build_inputs/,
# panel.k31.tsv.gz, controls.bed and controls.fa.gz. sinks.bed and anchors.json need sequencing data
# (see the end of this file). bundle.json, features.tsv and README.md are written by hand.
#
# Needs: the ngs-dose binary, python3 + numpy, curl, samtools, minimap2. ~4.2 GB of downloads
# (7.4 GB on disk), ~21 GB of memory (minimap2 indexes GRCh38), 5-15 min plus the downloads.
#
# The bundle is written to OUT (default work/bundle/GRCh38), not over resources/GRCh38: every counts
# file records the sha256 of the panel, controls and sinks it was made with, so a rebuild is compared
# with the shipped files instead of replacing them. OUT=resources/GRCh38 needs FORCE=1. From the
# inputs recorded below the rebuild reproduces the shipped units, build inputs, panel and controls
# byte for byte; $OUT/build_manifest.tsv lists every input and output with its sha256 and says, for
# each output, whether it is identical to the shipped file.
#
# Downloads go to a .part file that is moved into place only when complete, and every input is
# checked against the sha256 the shipped bundle was built from. VERIFY=0 accepts other copies (for
# example a local GRCh38 with other line widths); the manifest then records what was used.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SHIPPED="$ROOT/resources/GRCh38"
OUT="${OUT:-$ROOT/work/bundle/GRCh38}"
WORK="${WORK:-$ROOT/work/ref}"
BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"
mkdir -p "$OUT"
if [ "$(cd "$OUT" && pwd -P)" = "$(cd "$SHIPPED" && pwd -P)" ] && [ "${FORCE:-0}" != 1 ]; then
  echo "error: OUT is the shipped bundle ($SHIPPED). Build elsewhere and compare, or set FORCE=1 to replace it." >&2
  exit 1
fi
mkdir -p "$OUT/units" "$OUT/build_inputs" "$WORK"
MANIFEST="$OUT/build_manifest.tsv"
printf 'role\tfile\tsha256\tsource\n' > "$MANIFEST.part"

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }
# fetch_to FILE URL: download into FILE.part, and move it into place only once the transfer is complete
fetch_to() {
  echo "fetching $2" >&2
  curl -fsSL --retry 3 --retry-delay 5 -o "$1.part" "$2" || { rm -f "$1.part"; echo "error: download failed: $2" >&2; exit 1; }
  mv "$1.part" "$1"
}
# input FILE SHA256 SOURCE: check an input against the sha256 the shipped bundle was built from, and record it
input() {
  local got; got=$(sha256 "$1")
  if [ "$got" != "$2" ]; then
    if [ "${VERIFY:-1}" = 0 ]; then
      echo "note: $1 is not the recorded input (sha256 $got); used anyway (VERIFY=0)" >&2
    else
      echo "error: $1 has sha256 $got, expected $2. A partial or changed download? Delete it to fetch it again, or set VERIFY=0 to build from it anyway." >&2
      exit 1
    fi
  fi
  printf 'input\t%s\t%s\t%s\n' "$(basename "$1")" "$got" "$3" >> "$MANIFEST.part"
}
# contigs FASTA N NAME: the index holds N contigs, NAME among them (a truncated FASTA loses its last contigs)
contigs() {
  local n; n=$(wc -l < "$1.fai" | tr -d ' ')
  [ "$n" -eq "$2" ] && awk -v c="$3" '$1 == c {f = 1} END {exit !f}' "$1.fai" || { echo "error: $1.fai lists $n contigs, expected $2 including $3 (truncated FASTA, or a stale .fai?)" >&2; exit 1; }
}
# unit FILE LENGTH: a FASTA record of LENGTH bases, not an error page
unit() {
  [ "$(head -c1 "$1")" = ">" ] && [ "$(awk '!/^>/ {n += length($0)} END {print n + 0}' "$1")" -eq "$2" ] \
    || { echo "error: $1 is not a $2-bp FASTA record" >&2; exit 1; }
}

# ---- 1. references -----------------------------------------------------------------------
GRCH38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"
GRCH38_URL=https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
[ -s "$GRCH38" ] || fetch_to "$GRCH38" "$GRCH38_URL"
input "$GRCH38" 3b103f4742abfd54938fb0333e19ad067635c8eb86f1dbf0ce44b165c4292b50 "$GRCH38_URL"
[ -s "$GRCH38.fai" ] || samtools faidx "$GRCH38"
contigs "$GRCH38" 3366 chrEBV
CHM13="$WORK/chm13v2.0.fa"
CHM13_URL=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
if [ ! -s "$CHM13" ]; then
  [ -s "$CHM13.gz" ] || fetch_to "$CHM13.gz" "$CHM13_URL"
  input "$CHM13.gz" f274b0ee8bc1bca18d231d0eb80d764650d7f02d1fe07c43819e21168224304d "$CHM13_URL"
  gzip -t "$CHM13.gz"
  gzip -dc "$CHM13.gz" > "$CHM13.part" && mv "$CHM13.part" "$CHM13"
fi
input "$CHM13" 15a4ba1246f6021a89699bf5083da7f2bad3f79c86acd7bc1eb0ca3a13164e85 "$CHM13_URL, unpacked"
[ -s "$CHM13.fai" ] || samtools faidx "$CHM13"
contigs "$CHM13" 25 chrM
# NGS-PCA's exclusion set (10x SV blacklist + GEM 100-mer mappability < 1 + DGV + segmental duplications)
EXCLUDE_URL=https://raw.githubusercontent.com/jlanej/NGS-PCA/6543762f7df2147ade61e41da267721891d7c617/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.100.1.0.dgv.gsd.sorted.merge.bed.gz
if [ -n "${EXCLUDE:-}" ]; then
  [ -s "$EXCLUDE" ] || { echo "error: EXCLUDE=$EXCLUDE does not exist" >&2; exit 1; }
else
  EXCLUDE="$WORK/$(basename "$EXCLUDE_URL")"
  [ -s "$EXCLUDE" ] || fetch_to "$EXCLUDE" "$EXCLUDE_URL"
fi
input "$EXCLUDE" aec47097e622b0a49498f9f108150d564bdec614d96fce87f8edb35737b7d7de "$EXCLUDE_URL"
gzip -t "$EXCLUDE"

# ---- 2. unit sequences of the positional classes (GenBank) ---------------------------------
efetch_to() { fetch_to "$1" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$2&rettype=fasta&retmode=text"; }
U45="$OUT/units/rDNA45S.KY962518.1.fa"; U5="$OUT/units/rDNA5S.X12811.1.fa"
[ -s "$U45" ] || efetch_to "$U45" KY962518.1
unit "$U45" 44838                                                             # 45S unit
input "$U45" 4220cf78b947667135387098a58ac5396ecb9817e04b9a7739e7dee624ba248f "GenBank KY962518.1 (efetch)"
[ -s "$U5" ] || efetch_to "$U5" X12811.1
unit "$U5" 2231                                                               # 5S unit
input "$U5" a5c4e6314f848083d02a4a79d6a489036f0a344ae97ab39b3b02065017bf2afa "GenBank X12811.1 (efetch)"

# ---- 3. the distal junction (DJ): unit, core k-mers, GRCh38 loci ------------------------------
# CHM13v2.0: the 400 kb immediately distal to each rDNA array. chr21's copy is the unit.
DJ5="$OUT/build_inputs/CHM13v2.DJ.bed"
cat > "$DJ5" <<'BED'
chr13	5370548	5770548	DJ
chr14	1699537	2099537	DJ
chr15	2106442	2506442	DJ
chr21	2708298	3108298	DJ
chr22	4393794	4793794	DJ
BED
DJFA="$OUT/units/DJ.CHM13v2_chr21_2708299_3108298.fa"
[ -s "$DJFA" ] || { samtools faidx "$CHM13" chr21:2708299-3108298 | sed 's/^>.*/>DJ_chr21/' > "$DJFA.part" && mv "$DJFA.part" "$DJFA"; }
printf 'DJ\tpositional\t%s\t0\n' "$DJFA" > "$WORK/dj.classes.tsv"
grep -P '^chr21\t' "$DJ5" > "$WORK/dj.chr21.bed" 2>/dev/null || awk '$1=="chr21"' "$DJ5" > "$WORK/dj.chr21.bed"
"$BIN" panel -m "$WORK/dj.classes.tsv" -k 31 -b "$CHM13:$WORK/dj.chr21.bed" --report "$WORK/dj.report_a.tsv.gz" -o "$WORK/dj.a.panel.gz"
"$BIN" panel -m "$WORK/dj.classes.tsv" -k 31 -b "$CHM13:$DJ5"              --report "$WORK/dj.report_b.tsv.gz" -o "$WORK/dj.b.panel.gz"
python3 "$HERE/dj_core.py" "$WORK/dj.report_a.tsv.gz" "$WORK/dj.report_b.tsv.gz" DJ "$OUT/build_inputs/DJ.core.bed"
# GRCh38 holds ~1.9 Mb of DJ-like sequence in 23 pieces (chr21p, GL000220, KI270733, GL000195, GL000205, decoys)
python3 "$HERE/find_class_loci.py" chunk "$DJFA" 2000 1000 > "$WORK/dj.chunks.fa"
minimap2 -t 8 -c -x asm20 -N 50 -p 0.3 --secondary=yes "$GRCH38" "$WORK/dj.chunks.fa" > "$WORK/dj.chunks.GRCh38.paf"
python3 "$HERE/find_class_loci.py" loci "$WORK/dj.chunks.GRCh38.paf" DJ > "$WORK/dj.GRCh38.loci.bed"

# ---- 4. where each class legitimately lives, and the k-mer panel ------------------------------
# A panel k-mer must be absent from both genomes OUTSIDE these intervals. The rDNA loci were found
# by mapping 300-bp unit chunks with minimap2 (-x sr) and keeping array-scale loci; the dispersed
# orphan fragments (200-700 bp, 91-99.8% identity, e.g. chr1:91.39 Mb, chrX:109.05 Mb) are
# deliberately NOT exempted, so k-mers shared with them are dropped and the affected stretches
# of the unit (chiefly the 28S 3' end) become blind spots that the estimator masks exactly.
{ cat <<'BED'
chr21	8200000	8262000	rDNA45S
chr21	8384000	8474000	rDNA45S
chrUn_GL000220v1	100000	161802	rDNA45S
chr22_KI270733v1_random	117000	179772	rDNA45S
chr21	8986000	8989500	rDNA45S_fragment
chr1	228540000	228660000	rDNA5S
BED
cat "$WORK/dj.GRCh38.loci.bed"; } > "$OUT/build_inputs/GRCh38.class_loci.bed"
{ cat <<'BED'
chr13	5770548	9348041	rDNA_13_1
chr14	2099537	2817811	rDNA_14_1
chr15	2506442	4707485	rDNA_15_1
chr21	3108298	5612715	rDNA_21_1
chr22	4793794	5720650	rDNA_22_1
chr1	227730000	228040000	rDNA5S
BED
cat "$DJ5"; } > "$OUT/build_inputs/CHM13v2.class_loci.bed"

printf 'rDNA45S\tpositional\tunits/rDNA45S.KY962518.1.fa\t1\nrDNA5S\tpositional\tunits/rDNA5S.X12811.1.fa\t1\nDJ\tpositional\tunits/DJ.CHM13v2_chr21_2708299_3108298.fa\t0\t1\tbuild_inputs/DJ.core.bed\n' \
  > "$OUT/build_inputs/classes.tsv"
# manifest paths are relative to the manifest's own directory
sed -e "s|\tunits/|\t$OUT/units/|" -e "s|\tbuild_inputs/|\t$OUT/build_inputs/|" "$OUT/build_inputs/classes.tsv" > "$WORK/classes.abs.tsv"
# the panel header records each mask by its file name, and counts record the panel's sha256: the
# masks go in under the names the shipped panel carries, so that a rebuild reproduces it exactly
mkdir -p "$WORK/masks"
cp "$OUT/build_inputs/GRCh38.class_loci.bed" "$WORK/masks/grch38.class_loci.bed"
cp "$OUT/build_inputs/CHM13v2.class_loci.bed" "$WORK/masks/chm13.class_loci.bed"
"$BIN" panel -m "$WORK/classes.abs.tsv" -k 31 \
  -b "$GRCH38:$WORK/masks/grch38.class_loci.bed" \
  -b "$CHM13:$WORK/masks/chm13.class_loci.bed" \
  --report "$OUT/build_inputs/panel.k31.report.tsv.gz" -o "$OUT/panel.k31.tsv.gz"

# ---- 5. controls and known-truth regions -------------------------------------------------------
python3 "$HERE/select_controls.py" --reference "$GRCH38" --exclude "$EXCLUDE" --out "$OUT/controls.bed"
"$BIN" controls -b "$OUT/controls.bed" -T "$GRCH38" --flank 1000 -o "$OUT/controls.fa.gz"

# ---- 6. data-dependent pieces ------------------------------------------------------------------
# sinks.bed   two runs of `ngsdose sinks` on whole-file scans (ngs-dose count -m scan) of the target pipeline:
#             positional classes: ngsdose sinks NA12878.json.gz HG02258.json.gz > pos.bed  (compositional classes get no sinks by default)
#             TEL: scans made with -p resources/experimental/telomere.k31.panel.tsv.gz (372 cohort scans), then
#                  ngsdose sinks <scans> --classes TEL | awk '$4=="TEL"' > tel.bed
#             cat pos.bed tel.bed | sort -k1,1 -k2,2n > resources/GRCh38/sinks.bed   (bundle.json sinks_learned_from, sinks_history)
# anchors.json       NGS-DOSE-1000G: python pilot/evaluate_pilot.py --write-anchors   (needs cross-chemistry replicate pairs)
# efficiencies       optional, not part of the bundle: ngsdose cohort estimates/*.json.gz --save-efficiencies eff.json,
#                    applied later with ngsdose cohort ... --efficiencies eff.json

# ---- 7. build manifest: what was built from what, and how it compares with the shipped bundle ----
for f in units/rDNA45S.KY962518.1.fa units/rDNA5S.X12811.1.fa units/DJ.CHM13v2_chr21_2708299_3108298.fa \
         build_inputs/CHM13v2.DJ.bed build_inputs/DJ.core.bed build_inputs/GRCh38.class_loci.bed \
         build_inputs/CHM13v2.class_loci.bed build_inputs/classes.tsv panel.k31.tsv.gz controls.bed controls.fa.gz; do
  got=$(sha256 "$OUT/$f")
  if [ ! -e "$SHIPPED/$f" ]; then cmp="not in resources/GRCh38"
  elif [ "$got" = "$(sha256 "$SHIPPED/$f")" ]; then cmp="identical to resources/GRCh38"
  else cmp="differs from resources/GRCh38"
  fi
  printf 'output\t%s\t%s\t%s\n' "$f" "$got" "$cmp" >> "$MANIFEST.part"
done
{ printf '# %s  NGS-DOSE %s  %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(git -C "$ROOT" describe --always --dirty 2>/dev/null || echo unknown)" "$("$BIN" --version)"
  cat "$MANIFEST.part"; } > "$MANIFEST" && rm "$MANIFEST.part"
awk -F'\t' '$1 == "output"' "$MANIFEST"
echo "bundle written to $OUT (manifest: $MANIFEST)"
