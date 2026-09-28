#!/usr/bin/env bash
# CANDIDATE panel: non-human sequence in a human WGS file - the reads that align nowhere.
#
#   HHV6A, HHV6B  roseoloviruses that integrate into a telomere; 0/1/2 germline copies in
#                 iciHHV-6 carriers, about 1% of people. Positional, so the unit profile
#                 separates a whole integrated genome from a solo direct repeat.
#   HHV7          the third roseolovirus: no integration, but it shares its direct-repeat
#                 architecture with HHV-6 and must compete for those reads.
#   SMRV          squirrel monkey retrovirus, a B95-8/LCL culture contaminant with a published
#                 per-sample table for this cohort.
#   MYCO          the five Mycoplasma species of cell culture, as one class.
#   PHIX          the Illumina spike-in: a run/lane property and an index-hopping gauge.
#   EBV2          EBV type 2, from the type-divergent latency genes of AG876; the reference's
#                 chrEBV is type 1, so a type-2 virus in a donor is invisible today.
#
# All of these live in the unmapped bin of a GRCh38 CRAM (plus five hs38d1 decoys: one HHV-6A
# direct-repeat decoy and four Mycoplasmopsis arginini ones). A scan already k-mer-classifies
# unmapped reads, so loading this panel into a scan costs nothing but memory, and where the
# reads land has to be learned from those scans like any other class - nothing here assumes a
# read sits at the class's reference coordinates, because for six of the seven classes there
# are no reference coordinates at all.
#
# Backgrounds. Every class k-mer is dropped if it occurs even once in: the GRCh38 analysis set
# (primary, alt, HLA, the 2,385 hs38d1 decoys, chrM and chrEBV), CHM13v2.0, thirteen bacteria
# (lab, skin and reagent flora plus four Mollicutes that are not culture contaminants), sixteen
# other human viruses, and five EBV type-1 genomes. NO MASK is passed: the classes' only copies
# in a human reference are five decoy contigs, and dropping the k-mers they hold is cheaper
# than risking a human-derived k-mer in the panel. work/expected_loci.bed records those copies.
#
# Two filters are applied afterwards by tools/postfilter.py:
#   - every k-mer that the shipped panels (bundle, satellites, telomere) already claim, because
#     the engine drops a k-mer held by two panels from BOTH, which would change the counts of
#     the cohort now running;
#   - low-complexity k-mers (at most two distinct bases, or a homopolymer of 18 bp or more),
#     because the unmapped bin is full of two-colour poly-G reads.
# recipes/assemble.py then removes periodic k-mers (period 1-6 bp, <= 2 mismatches): the direct repeats
# of HHV-6 and HHV-7 hold telomere-repeat variants and microsatellites that GRCh38 and CHM13 lack but
# human reads carry (they gave NA12878, which carries neither virus, 283 HHV7 and 21 HHV6A reads).
# Needs: ngs-dose, curl, python3. Inputs: the GRCh38 analysis set and CHM13v2.0 in WORK, the shipped
# panels, and GenBank records (E-utilities; kept in SRC, default OUT/src, and not fetched again).
# Writes OUT/panel.tsv.gz (all seven classes; build_candidate_panels.sh splits MYCO into a file of its
# own), OUT/units/, OUT/rep.tsv; RECALL=1 also OUT/recall.tsv. Called by
# resources/build/build_candidate_panels.sh.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/nonhuman}"
SRC="${SRC:-$OUT/src}"; UNITS="$OUT/units"; W="$OUT/work"
REF="$WORK"; SHIPPED="$ROOT/resources"
mkdir -p "$SRC" "$UNITS" "$W"

# 1. sources (GenBank, via E-utilities)
CLASS_ACC="NC_001664.4 NC_000898.1 NC_001716.2 NC_001514.1 NC_001422.1 NC_009334.1"
# recall: other strains of each class, and the iciHHV-6 genomes of this cohort's own carriers
RECALL_ACC="KC465951.1 KP257584.1 MG894370.1 MG894374.1 KY316047.1 KY316049.1 \
AB021506.1 MG894368.1 MG894369.1 MG894372.1 MG894373.1 MH698401.1 MH698402.1 KY316045.1 \
U43400.1 M23385.1 MK561030.1 J02482.1 LN827800.1 \
NZ_CP064323.1 NZ_CP102738.1 NC_021002.1 NC_014552.1 NZ_CP143577.1 NZ_CP165990.1 NZ_CP055144.1 NZ_CP026341.1"
# background: EBV type 1 strains, and specificity checks for EBV2
EBV1_ACC="KC207813.1 KC207814.1 KF373730.1 AY961628.3 KF717093.1"
SPEC_ACC="OR134767.1 OR134769.1 OR134766.1 LN827799.1 NC_007605.1"
BG_ACC="NC_000913.3 NC_004461.1 NC_006085.1 NC_000964.3 NC_007795.1 NC_003028.3 NC_010682.1 \
NC_012660.1 NZ_CP103951.1 NC_000912.1 NC_000908.2 NC_002162.1 NC_006360.1 \
NC_006273.2 NC_001806.2 NC_001798.2 NC_001348.1 NC_009333.1 NC_001436.1 NC_001802.1 \
NC_000883.2 NC_001405.1 NC_003977.2 NC_001538.1 NC_001699.1 NC_010277.2 NC_001526.4 \
NC_001357.1 NC_001550.1"
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
fetch() { local d=$1; shift; mkdir -p "$d"; for a in "$@"; do
  [ -s "$d/$a.fa" ] || { curl -sf "$E?db=nuccore&id=$a&rettype=fasta&retmode=text" -o "$d/$a.fa.tmp" && mv "$d/$a.fa.tmp" "$d/$a.fa"; sleep 0.4; }; done; }
ft() { for a in "$@"; do [ -s "$SRC/$a.ft" ] || { curl -sf "$E?db=nuccore&id=$a&rettype=ft&retmode=text" -o "$SRC/$a.ft.tmp" && mv "$SRC/$a.ft.tmp" "$SRC/$a.ft"; sleep 0.4; }; done; }
# shellcheck disable=SC2086
fetch "$SRC" $CLASS_ACC $RECALL_ACC $EBV1_ACC $SPEC_ACC
# shellcheck disable=SC2086
fetch "$SRC/bg" $BG_ACC
# the roseolovirus repeat coordinates and the Mycoplasma rRNA/tRNA features come from the
# GenBank feature tables, not from a guess
ft NC_001664.4 NC_000898.1 NC_001716.2 NC_009334.1 \
   NC_014448.1 NZ_LR214940.1 NC_013511.1 NC_014921.1 NZ_AP014657.1
fetch "$SRC" NC_014448.1 NZ_LR214940.1 NC_013511.1 NC_014921.1 NZ_AP014657.1

# 2. class units
"$PY" "$HERE/tools/mkclasses.py" "$SRC" "$UNITS"

# 3. backgrounds
cat "$SRC"/bg/NC_000913.3.fa "$SRC"/bg/NC_004461.1.fa "$SRC"/bg/NC_006085.1.fa "$SRC"/bg/NC_000964.3.fa \
    "$SRC"/bg/NC_007795.1.fa "$SRC"/bg/NC_003028.3.fa "$SRC"/bg/NC_010682.1.fa "$SRC"/bg/NC_012660.1.fa \
    "$SRC"/bg/NZ_CP103951.1.fa "$SRC"/bg/NC_000912.1.fa "$SRC"/bg/NC_000908.2.fa "$SRC"/bg/NC_002162.1.fa \
    "$SRC"/bg/NC_006360.1.fa > "$W/bg_bacteria.fa"
cat "$SRC"/bg/NC_006273.2.fa "$SRC"/bg/NC_001806.2.fa "$SRC"/bg/NC_001798.2.fa "$SRC"/bg/NC_001348.1.fa \
    "$SRC"/bg/NC_009333.1.fa "$SRC"/bg/NC_001436.1.fa "$SRC"/bg/NC_001802.1.fa "$SRC"/bg/NC_000883.2.fa \
    "$SRC"/bg/NC_001405.1.fa "$SRC"/bg/NC_003977.2.fa "$SRC"/bg/NC_001538.1.fa "$SRC"/bg/NC_001699.1.fa \
    "$SRC"/bg/NC_010277.2.fa "$SRC"/bg/NC_001526.4.fa "$SRC"/bg/NC_001357.1.fa "$SRC"/bg/NC_001550.1.fa > "$W/bg_viruses.fa"
# shellcheck disable=SC2086
cat $(for a in $EBV1_ACC; do echo "$SRC/$a.fa"; done) > "$W/bg_ebv1.fa"

# 4. manifest and panel
cat > "$W/manifest.tsv" <<EOF
HHV6A	positional	$UNITS/HHV6A.fa	0
HHV6B	positional	$UNITS/HHV6B.fa	0
HHV7	positional	$UNITS/HHV7.fa	0
SMRV	positional	$UNITS/SMRV.fa	0
PHIX	positional	$UNITS/PHIX.fa	1
MYCO	compositional	$UNITS/MYCO.fa	0	1
EBV2	compositional	$UNITS/EBV2.fa	0	1
EOF
cp "$W/manifest.tsv" "$OUT/manifest.tsv"
"$BIN" panel -m "$W/manifest.tsv" -k 31 \
  -b "$REF/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  -b "$REF/chm13v2.0.fa" \
  -b "$W/bg_bacteria.fa" -b "$W/bg_viruses.fa" -b "$W/bg_ebv1.fa" \
  --max-bg 0 -o "$W/raw.panel.tsv.gz" --report "$W/rep.raw.tsv.gz"

# 5. shipped-panel overlap and low complexity
"$PY" "$HERE/tools/postfilter.py" "$W/raw.panel.tsv.gz" "$OUT/panel.tsv.gz" "$OUT/rep.tsv" \
  "$SHIPPED/GRCh38/panel.k31.tsv.gz" \
  "$SHIPPED/experimental/satellites.CHM13v2.k31.panel.tsv.gz" \
  "$SHIPPED/experimental/telomere.k31.panel.tsv.gz"

[ "${RECALL:-0}" = 1 ] || exit 0
# 7. recall: 150-bp reads tiled every 10 bp over every source copy that could be found, plus
#    EBV type-1 genomes (EBV2 must score 0 on them) and four bacteria left out of the background
"$PY" "$HERE/tools/mkrecall.py" "$SRC" "$W"
# shellcheck disable=SC2086
fetch "$SRC/bg" NC_016845.1 NC_004668.1 NC_004663.1 NC_021490.2
cat "$SRC"/bg/NC_016845.1.fa "$SRC"/bg/NC_004668.1.fa "$SRC"/bg/NC_004663.1.fa "$SRC"/bg/NC_021490.2.fa > "$W/offpanel_bacteria.fa"
"$PY" "$HERE/tools/recall.py" "$OUT/panel.tsv.gz" \
  "HHV6A=$W/recall_HHV6A.fa" "HHV6B=$W/recall_HHV6B.fa" "HHV7=$W/recall_HHV7.fa" \
  "SMRV=$W/recall_SMRV.fa" "PHIX=$W/recall_PHIX.fa" \
  "EBV2=$W/recall_EBV2.fa" "EBV2=$W/spec_EBV1.fa" "MYCO=$W/recall_MYCO.fa" > "$OUT/recall.tsv"
"$PY" "$HERE/tools/recall.py" "$OUT/panel.tsv.gz" --step 50 \
  "MYCO=$W/offpanel_bacteria.fa" "HHV6B=$W/offpanel_bacteria.fa" "PHIX=$W/offpanel_bacteria.fa" \
  > "$W/specificity.tsv"
