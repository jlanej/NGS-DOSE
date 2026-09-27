#!/usr/bin/env bash
# CANDIDATE panels: classes to be loaded in whole-file scans, whose sinks are not learned yet
# (resources/experimental/candidates/README.md). Rebuilds them from their sources:
#   1. each group's recipe (resources/experimental/candidates/recipes/<group>/build.sh) builds its
#      panel in OUT/<group>/ with the engine's panel builder, against the GRCh38 analysis set and
#      CHM13v2.0 in WORK, with the shipped panels' k-mers removed;
#   2. recipes/assemble.py checks that no candidate k-mer is in a shipped panel and no class name is
#      used twice, removes k-mers found in two candidate panels (from both, so that any combination
#      of panels loads with the same k-mers per class) and simple-repeat k-mers (periodic with a period
#      of 1-6 bp at <= 2 mismatches: telomere-repeat variants and microsatellites absent from the
#      references but present in every genome's reads), splits MYCO from the other non-human classes,
#      and writes DEST/<group>.k31.panel.tsv.gz (gzip without a timestamp: a rebuild gives the same
#      bytes, and the estimator matches experimental panels by their sha256);
#   3. the units go to DEST/units/<class>.fa, where `ngsdose estimate` looks for them.
# CANDIDATE_GROUPS picks the recipes to run (default all six); a group not rebuilt is taken from DEST as it is.
# RECALL=1 makes each recipe also measure read recall (more downloads and time).
# Needs: ngs-dose (cargo build --release), samtools, minimap2, curl, python3 with numpy.
# ~30 min for all six groups, peak ~3.2 GB (rna-arrays).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates}"; DEST="${DEST:-$ROOT/resources/experimental/candidates}"
REC="$ROOT/resources/experimental/candidates/recipes"
ALL_GROUPS="macrosatellites multicopy-genes sex-chromosome-arrays rna-arrays nonhuman coding-vntrs"
CANDIDATE_GROUPS="${CANDIDATE_GROUPS:-$ALL_GROUPS}"
for f in "$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" "$WORK/chm13v2.0.fa.fai" "$BIN"; do
  [ -s "$f" ] || { echo "error: $f missing (build_grch38_bundle.sh fetches the assemblies)" >&2; exit 1; }
done
export ROOT WORK NGSDOSE_BIN="$BIN" PY
mkdir -p "$OUT" "$DEST/units"

for g in $CANDIDATE_GROUPS; do
  echo "== $g" >&2
  OUT="$OUT/$g" bash "$REC/$g/build.sh" > "$OUT/$g.log" 2>&1 || { echo "error: $g failed, see $OUT/$g.log" >&2; exit 1; }
done

# the group panels: rebuilt ones from OUT, the others as they are in DEST
src() {  # src GROUP [CLASSES]: NAME=PANEL[:CLASSES] for assemble.py
  local g=$1 name=${3:-$1}
  if [[ " $CANDIDATE_GROUPS " == *" $g "* ]]; then echo "$name=$OUT/$g/panel.tsv.gz${2:+:$2}"; else echo "$name=$DEST/$name.k31.panel.tsv.gz"; fi
}
VIRAL=HHV6A,HHV6B,HHV7,SMRV,PHIX,EBV2
tmp="$OUT/assembled"; mkdir -p "$tmp"
"$PY" "$REC/assemble.py" \
  --shipped "$ROOT/resources/GRCh38/panel.k31.tsv.gz" "$ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz" \
            "$ROOT/resources/experimental/telomere.k31.panel.tsv.gz" \
  --out "$tmp" --shared "$tmp/shared_between_panels.tsv" \
  "$(src macrosatellites)" "$(src multicopy-genes)" "$(src sex-chromosome-arrays)" "$(src rna-arrays)" \
  "$(src nonhuman $VIRAL)" "$(src nonhuman MYCO nonhuman-myco)" "$(src coding-vntrs)" > "$OUT/assemble.tsv"
cp "$tmp"/*.k31.panel.tsv.gz "$tmp/shared_between_panels.tsv" "$DEST/"

# units: the sequence each class's k-mer positions refer to (the estimator computes each window's
# expected count from it). OPN1 is the real sequence: its panel was built with OPN1MW exon 5 set to N.
for g in $CANDIDATE_GROUPS; do
  for f in "$OUT/$g"/units/*.fa; do
    n=$(basename "$f" .fa)
    case $n in
      MYCO) gzip -9nc "$f" > "$DEST/units/MYCO.fa.gz" ;;
      OPN1) { echo ">OPN1 GRCh38 chrX:154182596-154219733 (OPN1MW gene start to the next unit); the panel was built with unit 10513-11352 (OPN1MW exon 5 +-300 bp) set to N"
              samtools faidx "$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa" chrX:154182596-154219733 | tail -n +2 | tr -d '\n' | fold -w 60; echo; } > "$DEST/units/OPN1.fa" ;;
      *) cp "$f" "$DEST/units/$n.fa" ;;
    esac
  done
done
cat "$OUT/assemble.tsv"
