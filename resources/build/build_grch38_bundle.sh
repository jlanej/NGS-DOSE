#!/usr/bin/env bash
# Rebuild the GRCh38 resource bundle from public inputs. Everything in resources/GRCh38/ is the
# output of this script (plus `ngsdose sinks` and `ngsdose cohort --save-efficiencies`, which
# need sequencing data; see the end of this file).
#
# Needs: the ngs-dose binary, python3 + numpy, curl, samtools, minimap2. ~8 GB of downloads, ~15 min.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="${OUT:-$ROOT/resources/GRCh38}"
WORK="${WORK:-$ROOT/work/ref}"
BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"
# NGS-PCA's exclusion set (10x SV blacklist + GEM 100-mer mappability < 1 + DGV + segmental duplications)
EXCLUDE="${EXCLUDE:-$HOME/git/NGS-PCA/resources/GRCh38/ngs_pca_exclude.sv_blacklist.map.kmer.100.1.0.dgv.gsd.sorted.merge.bed.gz}"
mkdir -p "$OUT/units" "$OUT/build_inputs" "$WORK"

# ---- 1. references -----------------------------------------------------------------------
GRCH38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CHM13="$WORK/chm13v2.0.fa"
[ -s "$GRCH38" ] || curl -sSL -o "$GRCH38" \
  https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
[ -s "$GRCH38.fai" ] || samtools faidx "$GRCH38"
[ -s "$CHM13" ] || { curl -sSL https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz | gzip -dc > "$CHM13"; }

# ---- 2. unit sequences of the positional classes (GenBank) ---------------------------------
efetch() { curl -sSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$1&rettype=fasta&retmode=text"; }
[ -s "$OUT/units/rDNA45S.KY962518.1.fa" ] || efetch KY962518.1 > "$OUT/units/rDNA45S.KY962518.1.fa"   # 44,838 bp 45S unit
[ -s "$OUT/units/rDNA5S.X12811.1.fa" ]   || efetch X12811.1   > "$OUT/units/rDNA5S.X12811.1.fa"       # 2,231 bp 5S unit

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
[ -s "$CHM13.fai" ] || samtools faidx "$CHM13"
[ -s "$DJFA" ] || samtools faidx "$CHM13" chr21:2708299-3108298 | sed 's/^>.*/>DJ_chr21/' > "$DJFA"
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
"$BIN" panel -m "$WORK/classes.abs.tsv" -k 31 \
  -b "$GRCH38:$OUT/build_inputs/GRCh38.class_loci.bed" \
  -b "$CHM13:$OUT/build_inputs/CHM13v2.class_loci.bed" \
  --report "$OUT/build_inputs/panel.k31.report.tsv.gz" -o "$OUT/panel.k31.tsv.gz"

# ---- 5. controls and known-truth regions -------------------------------------------------------
python3 "$HERE/select_controls.py" --reference "$GRCH38" --exclude "$EXCLUDE" --out "$OUT/controls.bed"
"$BIN" controls -b "$OUT/controls.bed" -T "$GRCH38" --flank 1000 -o "$OUT/controls.fa.gz"

# ---- 6. data-dependent pieces ------------------------------------------------------------------
# sinks.bed          ngs-dose count -m scan on a few whole CRAMs of the target pipeline, then
#                    ngsdose sinks scan1.json.gz [scan2.json.gz ...] -o resources/GRCh38/sinks.bed
# anchors.json       python example/1000G/pilot/evaluate_pilot.py --write-anchors   (needs cross-chemistry replicate pairs)
# efficiencies.json  ngsdose cohort estimates/*.json.gz --save-efficiencies resources/GRCh38/efficiencies.json
echo "bundle written to $OUT"
