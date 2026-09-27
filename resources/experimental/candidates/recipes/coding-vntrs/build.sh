#!/usr/bin/env bash
# Candidate panel "coding-vntrs": long-unit coding VNTRs from Mukamel et al. 2021 (Science 373:1499) and
# 2023 (Cell 186:3659), one COMPOSITIONAL class per locus (VNTR_<gene>), built from every allele the
# two reference assemblies hold: the GRCh38 primary array, the GRCh38 alt-contig arrays (MUC4: 7 alts,
# MUC6: 2 alts) and the CHM13v2.0 array.
#
# Why compositional, not positional: a positional class keeps only k-mers seen once in its unit
# (src/panel.rs), and the positional estimator measures copies of one consensus in 250-bp windows. A
# VNTR array is the opposite case: units 64-99% identical to each other (a single consensus unit
# misses the divergent ones), units of 29-99 bp for most loci (shorter than a window), and the
# quantity of interest is the array length. The compositional estimator gives it directly: diploid
# array mass in bp (ngsdose estimate: mass_bp = sum over read GC of reads / expected reads per bp),
# i.e. L1 + L2 plus a fixed edge term per allele (reads that overlap an array end by >= 34 bp carry
# the 4 k-mers a read needs), which recall.py measures per allele (eff_len_bp).
#
# Loci (sources/loci.GRCh38.tsv): Mukamel's 734 candidate coding VNTR regions (Zenodo 4776804,
# 734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt) with an array of >= 1 kb in GRCh38 and
# IBD2R >= 0.7 (ACAN, 0.69, as required), minus the segmental-duplication gene families (NBPF,
# ANKRD36, GOLGA6L, NPIPB, PRR20, SPDYE, PRB, KRAB-ZNF clusters, 22q11 LCR) and LPA (its own
# class). Array coordinates are the UCSC simpleRepeat (TRF) span, or for units beyond TRF's range
# (NEB, DMBT1, C2orf78) the span of a minimap2 self-alignment. The CHM13 array of each locus lies
# between its two 3-kb GRCh38 flanks mapped to CHM13 (minimap2 asm5, MAPQ 60, sources/loci.CHM13.tsv);
# alt arrays were found by mapping the flanks and the arrays to all 3,341 non-primary contigs
# (sources/loci.GRCh38_alts.tsv). Two hs38d1 decoys are pieces of a class's own array
# (sources/loci.GRCh38_decoys.tsv: chrUn_JTFH01000840v1_decoy is MUC19, chrUn_JTFH01000899v1_decoy MUC6;
# found by counting every source's 31-mers in all 3,341 non-primary contigs, where no other decoy,
# alt, random or chrUn contig is mostly made of a class's sequence). They are masked as the class's
# own copies, so the k-mers they share with the class are kept, but they are not class sources: a
# decoy is a fragment of a person's allele, not a reference allele of the array.
#
# A k-mer is kept if it is in no other class of this panel, nowhere in GRCh38 (analysis set: primary,
# alts, HLA, hs38d1 decoys, EBV) or CHM13 outside the arrays in masks/ (+-200 bp), and in none of the
# shipped panels (bundle panel, satellite panel, TEL panel): a k-mer shared with a loaded class is
# dropped by the engine from both classes, which would change the counts of the shipped ones.
#
# The loci in masks/ are where the reference assemblies hold the arrays. They are NOT sinks: where a
# pipeline puts these reads is learned from whole-file scans that carry this panel (ngsdose sinks
# --classes ...) and checked on held-out scans before any fetch relies on it.
#
# Needs: ngs-dose, samtools, python3. Inputs: the GRCh38 analysis set and CHM13v2.0 in WORK
# (build_grch38_bundle.sh fetches both) and NGS-DOSE's shipped panels. Writes OUT/panel.tsv.gz,
# OUT/units/, OUT/masks/, OUT/panel_stats.tsv; RECALL=1 also OUT/recall.tsv. ~4 min, ~2 GB.
# Called by resources/build/build_candidate_panels.sh.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/coding-vntrs}"
GRCH38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"; CHM13="$WORK/chm13v2.0.fa"
SHIPPED=("$ROOT/resources/GRCh38/panel.k31.tsv.gz" "$ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz" "$ROOT/resources/experimental/telomere.k31.panel.tsv.gz")
PAD=200
T="$OUT/tmp"; mkdir -p "$OUT/units" "$OUT/masks" "$T"
SRC="$HERE/sources"
cd "$OUT"

# ---- 1. class sources: every reference allele of each array ----------------------------------------
# loci.GRCh38.tsv: gene chrom start end unit_bp how (0-based, half-open); loci.CHM13.tsv: gene chrom
# start end strand; loci.GRCh38_alts.tsv: gene contig start end how
: > manifest.tsv
while IFS=$'\t' read -r g c s e u how; do
  f="units/VNTR_$g.fa"
  samtools faidx "$GRCH38" "$c:$((s + 1))-$e" | sed "1s/.*/>${g}_GRCh38_${c}_$((s + 1))_$e/" > "$f"
  awk -F'\t' -v g="$g" '$1 == g' "$SRC"/loci.GRCh38_alts.tsv | while IFS=$'\t' read -r _ ac as ae _; do
    samtools faidx "$GRCH38" "$ac:$((as + 1))-$ae" | sed "1s/.*/>${g}_GRCh38_${ac}_$((as + 1))_$ae/"
  done >> "$f"
  awk -F'\t' -v g="$g" '$1 == g' "$SRC"/loci.CHM13.tsv | while IFS=$'\t' read -r _ hc hs he _; do
    samtools faidx "$CHM13" "$hc:$((hs + 1))-$he" | sed "1s/.*/>${g}_CHM13_${hc}_$((hs + 1))_$he/"
  done >> "$f"
  [ "$(grep -c '^>' "$f")" -ge 2 ] || { echo "error: $f lacks its CHM13 allele" >&2; exit 1; }
  if grep -v '^>' "$f" | grep -qi n; then echo "error: $f contains N" >&2; exit 1; fi
  printf 'VNTR_%s\tcompositional\t%s\t0\t1\n' "$g" "$f" >> manifest.tsv
done < "$SRC"/loci.GRCh38.tsv

# ---- 2. masks: the arrays (+-PAD), every copy in each assembly; the own-copy decoys whole -----------
{ awk -F'\t' -v OFS='\t' -v p=$PAD '{print $2, ($3 > p ? $3 - p : 0), $4 + p, "VNTR_" $1}' "$SRC"/loci.GRCh38.tsv
  awk -F'\t' -v OFS='\t' -v p=$PAD '{print $2, ($3 > p ? $3 - p : 0), $4 + p, "VNTR_" $1 "_alt"}' "$SRC"/loci.GRCh38_alts.tsv
  awk -F'\t' -v OFS='\t' '{print $2, $3, $4, "VNTR_" $1 "_decoy"}' "$SRC"/loci.GRCh38_decoys.tsv
} | sort -k1,1 -k2,2n > masks/GRCh38.coding_vntr_loci.bed
awk -F'\t' -v OFS='\t' -v p=$PAD '{print $2, $3 - p, $4 + p, "VNTR_" $1}' "$SRC"/loci.CHM13.tsv | sort -k1,1 -k2,2n > masks/CHM13.coding_vntr_loci.bed

# ---- 3. the shipped panels' k-mers as one more background: any shared k-mer is dropped here -------
for p in "${SHIPPED[@]}"; do gzip -dc "$p" | awk '!/^#/ {print $1}'; done \
  | awk 'BEGIN {print ">shipped_panel_kmers"} {printf "%s%s", (NR > 1 ? "N" : ""), $1} END {print ""}' > "$T/shipped_kmers.fa"

# ---- 4. panels: A without, B with the shipped k-mers (B is the candidate; A only counts what B removed)
"$BIN" panel -m manifest.tsv -k 31 --max-bg 0 \
  -b "$GRCH38:masks/GRCh38.coding_vntr_loci.bed" -b "$CHM13:masks/CHM13.coding_vntr_loci.bed" \
  --report "$T/rep_A.tsv.gz" -o "$T/panel.noshipped.tsv.gz" 2> "$T/panel_A.log"
"$BIN" panel -m manifest.tsv -k 31 --max-bg 0 \
  -b "$GRCH38:masks/GRCh38.coding_vntr_loci.bed" -b "$CHM13:masks/CHM13.coding_vntr_loci.bed" -b "$T/shipped_kmers.fa" \
  --report "$T/rep.tsv.gz" -o panel.tsv.gz 2> "$T/panel_B.log"
gzip -dc "$T/rep.tsv.gz" > rep.tsv
"$PY" "$HERE/panel_stats.py" panel.tsv.gz "$T/panel.noshipped.tsv.gz" rep.tsv "${SHIPPED[@]}" > panel_stats.tsv
cat panel_stats.tsv

[ "${RECALL:-0}" = 1 ] || exit 0
# ---- 5. read recall: 150-bp reads every 10 bp over each source allele, share with >= 4 class k-mers,
#         with the full panel and with a panel built from the GRCh38 primary allele alone (held-out
#         alleles: CHM13, alts) - how well the k-mers of one allele see another ----------------------
: > "$T/manifest.grch38only.tsv"
for f in units/VNTR_*.fa; do
  n=$(basename "$f" .fa)
  awk '/^>/ {p = (NR == 1)} p' "$f" > "$T/$n.grch38only.fa"
  printf '%s\tcompositional\t%s\t0\t1\n' "$n" "$T/$n.grch38only.fa" >> "$T/manifest.grch38only.tsv"
done
"$BIN" panel -m "$T/manifest.grch38only.tsv" -k 31 --max-bg 0 \
  -b "$GRCH38:masks/GRCh38.coding_vntr_loci.bed" -b "$CHM13:masks/CHM13.coding_vntr_loci.bed" -b "$T/shipped_kmers.fa" \
  -o "$T/panel.grch38only.tsv.gz" 2> "$T/panel_L.log"
export GRCH38 CHM13
"$PY" "$HERE/recall.py" panel.tsv.gz "$T/panel.grch38only.tsv.gz" > recall.tsv
cat recall.tsv
