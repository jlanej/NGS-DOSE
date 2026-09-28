#!/usr/bin/env bash
# Candidate panel "macrosatellites": tandem macrosatellites and RNA-gene arrays, plus the D4Z4 distal
# haplotypes. Twelve positional classes, one unit record each:
#   DXZ4 (Xq23), CT47 (Xq24), RS447 (USP17L, 4p16.1), MSR5p (TAF11-like, 5p15.1),
#   FLJ40296 (PRR20, 13q21.1), RNU2 (U2 snRNA, 17q21.31), D4Z4 (4q35 + 10q26 units, one class),
#   ZAV (9q32), REXO1L1 (8q21.2): circular tandem units;
#   D4Z4_4qA, D4Z4_4qB, D4Z4_10q: sequence directly distal to the last D4Z4 unit (4qA and 10q: 14 kb,
#   pLAM + beta satellite + distal; 4qB: the qB segment up to the telomere repeat, 29 kb), not circular.
# A k-mer is kept if it occurs once in its unit, in no other class of this panel, nowhere in GRCh38
# (analysis set: primary, alts, HLA, hs38d1 decoys, EBV) or CHM13 outside the class's own copies
# (masks/), and in none of the shipped panels (bundle, satellite and TEL panels): a k-mer shared with
# a loaded class is dropped by the engine from both classes, which would change the counts of the shipped ones.
# The masks list the class's own copies only (arrays, partial units beside them, the inverted D4Z4 copy
# and its chr4 alt): diverged paralogs (chr8 USP17L, CT47C1, the DXZ4-like fragment, D4Z4-like decoys
# and acrocentric copies, the 4p copy of the 4qB segment) stay in the background, so their k-mers go.
#
# The loci in masks/ are where the reference assemblies hold copies. They are NOT sinks: where a
# pipeline puts these reads is learned from whole-file scans that carry this panel (ngsdose sinks) and
# checked on held-out scans before any fetch relies on it.
#
# Needs: ngs-dose, samtools, curl, python3 (minimap2 for RECALL=1). Inputs: the GRCh38 analysis set and
# CHM13v2.0 in WORK (build_grch38_bundle.sh fetches both), the shipped panels, and three GenBank records
# (E-utilities). Writes OUT/panel.tsv.gz, OUT/units/, OUT/masks/, OUT/panel_stats.tsv. ~11 min (five
# panel builds, three at a time), ~1 GB per build. Called by resources/build/build_candidate_panels.sh.
# RECALL=1 also measures read recall (OUT/recall.tsv); the HPRC r2 haplotype rows of recall_sources.tsv
# need extracts in OUT/sources/hprc/ (made by the literature assessments; skipped when absent).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/macrosatellites}"
GRCH38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"; CHM13="$WORK/chm13v2.0.fa"
SHIPPED=("$ROOT/resources/GRCh38/panel.k31.tsv.gz" "$ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz" "$ROOT/resources/experimental/telomere.k31.panel.tsv.gz")
T="$OUT/tmp"; mkdir -p "$OUT/units" "$OUT/masks" "$OUT/sources" "$T"
cd "$OUT"

# ---- 1. units -----------------------------------------------------------------------------------
efetch() { curl -fsSL --retry 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$1&rettype=fasta&retmode=text"; }
[ -s sources/U57614.1.fa ] || efetch U57614.1 > sources/U57614.1.fa
[ -s sources/AF117653.3.fa ] || efetch AF117653.3 > sources/AF117653.3.fa
# unit NAME RECORD-NAME FASTA REGION EXPECTED-LENGTH
unit() {
  samtools faidx "$3" "$4" | sed "1s/.*/>$2/" > "units/$1.fa"
  local n; n=$(awk '!/^>/ {n += length($0)} END {print n + 0}' "units/$1.fa")
  [ "$n" -eq "$5" ] || { echo "error: units/$1.fa is $n bp, expected $5" >&2; exit 1; }
}
# tandem units: one period of the array, the next unit starting exactly where this one ends (checked:
# the unit's first 31-mer recurs exactly one unit length on, and every 31-mer across the circular join
# occurs in the array; RNU2 is the GenBank clone less its last 6 bp, which repeat its first 6, the
# HindIII site AAGCTT at both ends of the cloned fragment. REXO1L1 is the normal 12,194-bp period that
# ends where the atypical 15,275-bp unit at chr8:85,772,425 begins: that unit carries a 2.9-kb internal
# duplication which about 8 of CHM13's 70 units have)
unit DXZ4     DXZ4_CHM13_chrX_114191115_114194086     "$CHM13"  chrX:114191115-114194086   2972
unit CT47     CT47_CHM13_chrX_119264861_119269720     "$CHM13"  chrX:119264861-119269720   4860
unit RS447    RS447_GRCh38_chr4_9215000_9219746       "$GRCH38" chr4:9215000-9219746       4747
unit MSR5p    MSR5p_GRCh38_chr5_17518000_17521433     "$GRCH38" chr5:17518000-17521433     3434
unit FLJ40296 FLJ40296_GRCh38_chr13_57143000_57149569 "$GRCH38" chr13:57143000-57149569    6570
unit ZAV      ZAV_GRCh38_chr9_113064500_113069947     "$GRCH38" chr9:113064500-113069947   5448
unit REXO1L1  REXO1L1_GRCh38_chr8_85760231_85772424   "$GRCH38" chr8:85760231-85772424     12194
unit RNU2     U57614.1_1_6126                         sources/U57614.1.fa  U57614.1:1-6126 6126
unit D4Z4     AF117653.3_5743_9038                    sources/AF117653.3.fa AF117653.3:5743-9038 3296
# distal segments (not circular), from the base after the last D4Z4 unit (4qA: pLAM, beta satellite, the
# first ~7 kb of the 4qA subtelomere; 10q: the same on chr10; 4qB: the qB segment and subtelomere up to
# the telomere repeat). Every chr10 copy in both assemblies carries the 10qA polyadenylation signal
# ATCAAA, the chr4 4qA copies ATTAAA (checked in this build's scratch).
unit D4Z4_4qA D4Z4_4qA_GRCh38_chr4_190175604_190189603   "$GRCH38" chr4:190175604-190189603   14000
unit D4Z4_4qB D4Z4_4qB_GRCh38_chr4_190093265_190122572   "$GRCH38" chr4:190093265-190122572   29308
unit D4Z4_10q D4Z4_10q_GRCh38_chr10_133761681_133775680  "$GRCH38" chr10:133761681-133775680  14000

cat > manifest.tsv <<'EOF'
DXZ4	positional	units/DXZ4.fa	1
CT47	positional	units/CT47.fa	1
RS447	positional	units/RS447.fa	1
MSR5p	positional	units/MSR5p.fa	1
FLJ40296	positional	units/FLJ40296.fa	1
RNU2	positional	units/RNU2.fa	1
D4Z4	positional	units/D4Z4.fa	1
ZAV	positional	units/ZAV.fa	1
REXO1L1	positional	units/REXO1L1.fa	1
D4Z4_4qA	positional	units/D4Z4_4qA.fa	0
D4Z4_4qB	positional	units/D4Z4_4qB.fa	0
D4Z4_10q	positional	units/D4Z4_10q.fa	0
EOF

# ---- 2. every copy of each class in the two assemblies (BED, 0-based, name = class). Found with minimap2
#         asm20 of 1-kb unit chunks against 0.4-2.7 Mb windows around each locus; copies of >= 97% identity
#         within 50 kb of the array count as the class's own (partial units, the chr5:17.63 Mb MSR5p and
#         chr17:43.32 Mb RNU2 fragments, the inverted D4Z4 copy on chr4 and on chr4_KI270896v1_alt); 2 kb
#         padding at array ends. No alt, HLA or decoy contig holds any other copy (the panel report would
#         show a class losing most of its k-mers). -----------------------------------------------------
cat > masks/GRCh38.class_loci.bed <<'EOF'
chr4	9208250	9371058	RS447
chr4	190021140	190023276	D4Z4
chr4	190063228	190093264	D4Z4
chr4	190093264	190122572	D4Z4_4qB
chr4	190171121	190175603	D4Z4
chr4	190175603	190189603	D4Z4_4qA
chr4_KI270896v1_alt	358671	360807	D4Z4
chr5	17515431	17601621	MSR5p
chr5	17629296	17634720	MSR5p
chr8	85638339	85830705	REXO1L1
chr9	113057718	113092511	ZAV
chr10	133662429	133685491	D4Z4
chr10	133685491	133690466	D4Z4_10q
chr10	133738608	133761680	D4Z4
chr10	133761680	133775680	D4Z4_10q
chr13	57137162	57175450	FLJ40296
chr17	43229443	43310000	RNU2
chr17	43320321	43325825	RNU2
chrX	115723929	115891415	DXZ4
chrX	120870714	120989017	CT47
EOF
cat > masks/CHM13.class_loci.bed <<'EOF'
chr4	9121235	9391581	RS447
chr4	193389358	193391494	D4Z4
chr4	193431357	193543408	D4Z4
chr4	193543408	193557408	D4Z4_4qA
chr5	17454872	17601835	MSR5p
chr5	17629506	17635184	MSR5p
chr8	86077107	86948157	REXO1L1
chr9	125229476	125291380	ZAV
chr10	134613846	134726331	D4Z4
chr10	134726331	134743331	D4Z4_10q
chr13	56354034	56392314	FLJ40296
chr17	44088297	44095163	RNU2
chr17	44112173	44171769	RNU2
chrX	114131487	114301042	DXZ4
chrX	119246389	119288877	CT47
EOF

# ---- 3. panels. The engine takes one mask per background genome and exempts it for every class, so the
#         D4Z4 unit and the three distal classes, whose loci adjoin and are 92-99% alike, are each built with
#         only their own copies masked (a 4qA k-mer found in the 10qA distal is then dropped, and the other
#         way round); the eight unrelated tandem arrays are built together. Every build takes the whole
#         manifest, so a k-mer in two units is dropped from both whichever build a class comes from. ------
TANDEM="DXZ4 CT47 RS447 MSR5p FLJ40296 RNU2 ZAV REXO1L1"
build() {  # build NAME CLASS...
  local name=$1; shift
  for a in GRCh38 CHM13; do
    awk -v c=" $* " 'index(c, " " $4 " ")' "masks/$a.class_loci.bed" > "masks/$a.$name.bed"
  done
  "$BIN" panel -m manifest.tsv -k 31 --max-bg 0 \
    -b "$GRCH38:masks/GRCh38.$name.bed" -b "$CHM13:masks/CHM13.$name.bed" \
    --report "$T/rep.$name.tsv.gz" -o "$T/panel.$name.tsv.gz" 2> "$T/panel.$name.log"
}
build tandem $TANDEM & build D4Z4 D4Z4 & build D4Z4_4qA D4Z4_4qA & wait
build D4Z4_4qB D4Z4_4qB & build D4Z4_10q D4Z4_10q & wait

# ---- 4. one panel: each class from its build, minus every k-mer of the shipped panels (a k-mer shared
#         with a loaded class is dropped by the engine from both, which would change the shipped counts) --
"$PY" "$HERE/merge_panels.py" manifest.tsv panel.tsv.gz rep.tsv \
  "$T/panel.tandem.tsv.gz:$T/rep.tandem.tsv.gz:$TANDEM" \
  "$T/panel.D4Z4.tsv.gz:$T/rep.D4Z4.tsv.gz:D4Z4" \
  "$T/panel.D4Z4_4qA.tsv.gz:$T/rep.D4Z4_4qA.tsv.gz:D4Z4_4qA" \
  "$T/panel.D4Z4_4qB.tsv.gz:$T/rep.D4Z4_4qB.tsv.gz:D4Z4_4qB" \
  "$T/panel.D4Z4_10q.tsv.gz:$T/rep.D4Z4_10q.tsv.gz:D4Z4_10q" \
  -- "${SHIPPED[@]}" > panel_stats.tsv
cat panel_stats.tsv

[ "${RECALL:-0}" = 1 ] || exit 0
# the 37 GenBank 4q/10q distal haplotypes of Lemmers et al. 2010 (Science 329:1650)
[ -s sources/HM190160-196.fa ] || efetch "$(seq 160 196 | sed 's/^/HM190/; s/$/.1/' | paste -sd, -)" > sources/HM190160-196.fa
# ---- 5. read recall: 150-bp reads every 10 bp over each source copy, share with >= 4 class k-mers ----
# recall_sources.tsv lists the copies: CHM13 and GRCh38 arrays, GenBank records, HPRC r2 haplotype
# extracts in sources/hprc, and paralogs that must stay unassigned (the D4Z4-like loci of both
# assemblies: 81-92% identity pieces found by minimap2 in the D4Z4 assessment, regions in sources/)
samtools faidx "$CHM13" -r "$HERE/sources/D4Z4like.CHM13.regions.txt" > "$T/D4Z4like.CHM13.fa"
samtools faidx "$GRCH38" -r "$HERE/sources/D4Z4like.GRCh38.regions.txt" > "$T/D4Z4like.GRCh38.fa"
export GRCH38 CHM13
# rows whose FASTA is missing (the HPRC extracts) are left out
awk -F'\t' '/^#/ || $3 ~ /^\$/ || system("test -s \"" $3 "\"") == 0' "$HERE/recall_sources.tsv" > recall_sources.tsv
"$PY" "$HERE/recall.py" panel.tsv.gz recall_sources.tsv > recall.tsv
cat recall.tsv
