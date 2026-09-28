#!/usr/bin/env bash
# Candidate panel "multicopy-genes": multi-copy and multi-allelic gene families as positional
# k-mer classes, for whole-file scans first (where the reads land is learned from the scans, not
# assumed from these coordinates). From the sources:
#   classes.tsv    the units (GRCh38 coordinates, extracted with samtools faidx) and build modes
#   masks/         every copy of each class in GRCh38 and CHM13 (masks/<class>.<assembly>.bed), the
#                  background masks: made by masks.py from a minimap2 search (DERIVE_MASKS=1 redoes it)
#   GRCh38 analysis set + CHM13v2.0 (backgrounds), the shipped panels (k-mers that must not recur)
# Steps:
#   1. OUT/units/<class>.fa (and OUT/aux/ paralogs) from classes.tsv           (build.py --units-only)
#   2. DERIVE_MASKS=1 only: every copy of every unit in both assemblies, minimap2 with the UNITS as
#      the index and 1-Mb windows of each genome as queries (no genome index is built, ~20 min per
#      assembly), then OUT/masks/ by per-class identity/length rules (masks.py). A rerun is expected,
#      not verified, to give the shipped masks/ (they were made from search_units.fa plus a few
#      paralog genes).
#   3. stage 1: one engine build per class against both assemblies with its own masks, --max-bg 0;
#      stage 2: all classes together (keep BEDs = stage-1 k-mers), with the shipped panels' k-mers
#      as the only background: OUT/panel.tsv.gz, OUT/rep.tsv, OUT/manifest.tsv            (build.py)
#   4. RECALL=1 only: recall per source copy and leakage from paralogs              (recall.sh)
# Needs: ngs-dose, samtools, python3 (minimap2 for DERIVE_MASKS=1; bedtools for RECALL=1). ~10 min.
# Called by resources/build/build_candidate_panels.sh.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; BIN="${NGSDOSE_BIN:-$ROOT/target/release/ngs-dose}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/multicopy-genes}"
G38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"; CHM="$WORK/chm13v2.0.fa"
SHIPPED=("$ROOT/resources/GRCh38/panel.k31.tsv.gz" "$ROOT/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz" "$ROOT/resources/experimental/telomere.k31.panel.tsv.gz")
mkdir -p "$OUT"
cd "$OUT"

"$PY" "$HERE/build.py" --grch38 "$G38" --chm13 "$CHM" --engine "$BIN" --shipped "${SHIPPED[@]}" --out "$OUT" --units-only

MASKS="$HERE/masks"
if [ "${DERIVE_MASKS:-0}" = 1 ]; then
  mkdir -p copies
  # the search units: every unit of classes.tsv except SMN1 (the same sequence as SMN) and the three
  # sub-units searched in their loci only (HERVC4 inside C4; HBA, HPdup: masks.py REGIONAL); CCL3L is
  # searched under the name CCLblk
  for f in units/*.fa aux/*.fa; do
    case $(basename "$f" .fa) in SMN1|HERVC4|HBA|HPdup) continue;; esac
    sed 's/^>CCL3L[[:space:]].*/>CCLblk/' "$f"
  done > copies/search_units.fa
  for asm in GRCh38 CHM13; do
    fa=$G38; [ $asm = CHM13 ] && fa=$CHM
    [ -s copies/$asm.paf ] && continue
    awk -v OFS='\t' '{print $1,$2}' "$fa.fai" \
      | awk '{for (s = 0; s < $2; s += 800000) {e = s + 1000000; if (e > $2) e = $2; printf "%s:%d-%d\n", $1, s + 1, e; if (e == $2) break}}' > copies/$asm.win.txt
    samtools faidx "$fa" -r copies/$asm.win.txt \
      | minimap2 -t 4 -x asm20 -c --cs -N 20 copies/search_units.fa - 2> copies/$asm.mm.log \
      | awk '$11 >= 500' > copies/$asm.paf
  done
  "$PY" "$HERE/masks.py" copies/GRCh38.paf copies/CHM13.paf "$G38" "$CHM" masks
  MASKS="$OUT/masks"
fi

"$PY" "$HERE/build.py" --grch38 "$G38" --chm13 "$CHM" --engine "$BIN" --shipped "${SHIPPED[@]}" --out "$OUT" --masks "$MASKS" -j 5
[ "${RECALL:-0}" = 1 ] && MASKS="$MASKS" OUT="$OUT" ROOT="$ROOT" WORK="$WORK" PY="$PY" bash "$HERE/recall.sh"
echo "panel: $OUT/panel.tsv.gz"
