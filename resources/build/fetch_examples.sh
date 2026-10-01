#!/usr/bin/env bash
# Worked examples of `ngsdose fetchplan`, written to docs/fetch_examples.md: the command line of each,
# the plan table it prints, what it writes to PREFIX.count_flags.txt, its notes, and its total MB per
# genome (median over the CRAM indexes given). The READMEs and docs/DESIGN.md quote plan figures from
# that file, so rerun this after changing the menu, a sinks or statistics file, the control regions or
# fetchplan.
#
# usage: resources/build/fetch_examples.sh CRAI...
#        resources/build/fetch_examples.sh --get DIR
#   CRAI      CRAM indexes to cost on. The docs' file uses the 13 below: NYGC 30x CRAMs of the 1000 Genomes
#             high-coverage release (bwa-mem alignments to GRCh38_full_analysis_set_plus_decoy_hla), whose
#             .crai sit beside the CRAMs named in the release's two sequence indexes (the manifest of
#             NGS-DOSE-1000G's pipeline/00_setup.sh). --get DIR downloads those 13 (about 19 MB) and the
#             reference's .fai into DIR, keeping files already there, and uses them.
# environment:
#   CONTIGS   the contigs of the CRAMs in header order (default work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai,
#             or DIR's copy with --get)
#   ENGINE    the ngs-dose that would run the fetch (default target/release/ngs-dose); it decides what
#             count_flags.txt holds (--classes for an engine that takes it)
#   PYTHON    the interpreter that runs this checkout's ngsdose (default python3). The script runs it from
#             the checkout's root with the checkout first on PYTHONPATH, so neither an installed ngsdose nor
#             one in the caller's working directory is imported instead
#   OUT       the Markdown file to write (default docs/fetch_examples.md)
# Needs: python3 with the ngsdose dependencies, and curl for --get. About 30 s.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
usage() { sed -n '2,25p' "$0" >&2; exit 2; }
[ $# -ge 1 ] || usage
S3=https://1000genomes.s3.amazonaws.com
URLS=(  # the 13 indexes: CRAM URL + .crai
  "$S3/1000G_2504_high_coverage/data/ERR3240114/HG00096.final.cram.crai"
  "$S3/1000G_2504_high_coverage/additional_698_related/data/ERR3988821/HG00706.final.cram.crai"
  "$S3/1000G_2504_high_coverage/additional_698_related/data/ERR3988835/HG01084.final.cram.crai"
  "$S3/1000G_2504_high_coverage/additional_698_related/data/ERR3988903/HG01552.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3242096/HG01871.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3242274/HG02383.final.cram.crai"
  "$S3/1000G_2504_high_coverage/additional_698_related/data/ERR3989016/HG02451.final.cram.crai"
  "$S3/1000G_2504_high_coverage/additional_698_related/data/ERR3989018/HG02466.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3243025/HG02792.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3242643/HG03007.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3242561/HG03072.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3242588/HG03267.final.cram.crai"
  "$S3/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram.crai"
)
FAI_URL="$S3/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
get() {  # URL FILE: download unless present; the name is taken only once whole
  [ -s "$2" ] || { curl -fsSL --retry 5 -o "$2.part" "$1" && mv "$2.part" "$2"; } || { echo "error: could not download $1" >&2; exit 1; }
}
if [ "$1" = "--get" ]; then
  [ $# -eq 2 ] || usage
  mkdir -p "$2"; set -- "$(cd "$2" && pwd)"; d="$1"; shift
  for u in "${URLS[@]}"; do f="${u##*/}"; get "$u" "$d/${f%.final.cram.crai}.crai"; set -- "$@" "$d/${f%.final.cram.crai}.crai"; done
  get "$FAI_URL" "$d/${FAI_URL##*/}"; CONTIGS="${CONTIGS:-$d/${FAI_URL##*/}}"
fi
CONTIGS="${CONTIGS:-$ROOT/work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai}"
ENGINE="${ENGINE:-$ROOT/target/release/ngs-dose}"
OUT="${OUT:-$ROOT/docs/fetch_examples.md}"
MENU="$ROOT/resources/fetch_menu.tsv"
abs() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$PWD" "$1" ;; esac; }
CONTIGS="$(abs "$CONTIGS")"; ENGINE="$(abs "$ENGINE")"; OUT="$(abs "$OUT")"
CRAI=(); for c in "$@"; do CRAI+=("$(abs "$c")"); done
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
NGSDOSE="${PYTHON:-python3} -m ngsdose.cli"
for f in "$CONTIGS" "$ENGINE" "$MENU" "${CRAI[@]}"; do [ -s "$f" ] || { echo "error: $f missing" >&2; exit 1; }; done
# python -m puts the working directory ahead of PYTHONPATH: run from the checkout so its ngsdose is used
cd "$ROOT"
if command -v sha256sum >/dev/null; then SHA=(sha256sum); else SHA=(shasum -a 256); fi
L="$ROOT/resources/GRCh38/controls.lite200.bed"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# name|arguments (word-split; no argument holds a space)
EXAMPLES=(
  "core: the bundle's positional classes|--preset core"
  "core_tel: core and the telomeric repeat|--preset core_tel"
  "biobank_lite: core_tel and the four cheapest satellites|--preset biobank_lite"
  "core and the cheap satellites SST1, CER, SATR, ACRO at a capture target of 0.995|--preset core --classes SST1 CER SATR ACRO --capture 0.995"
  "core with the 200 lite control regions|--preset core --controls $L"
  "core_tel with a capture target for TEL alone|--preset core_tel --capture-class TEL=0.99"
  "the per-array options DXZ1 and DYZ3|--classes DXZ1 DYZ3"
  "core_tel and the tier-A candidates, for the whole-file scans|--preset core_tel candidates_A"
  "a budget of 600 MB, filled|--budget-mb 600 --fill"
  "a budget of 1000 MB, filled|--budget-mb 1000 --fill"
  "the ten satellite families and the unmapped bin at a capture target of 0.995|--preset satellites --capture 0.995"
  "every option with sinks (each row's mb_median is that option alone)|--preset core_tel satellites xy_arrays"
  "core_tel with the lite control regions and a capture target of 0.995|--preset core_tel --controls $L --capture 0.995"
  "biobank_lite at a capture target of 0.995|--preset biobank_lite --capture 0.995"
  "the four per-array options (preset xy_arrays)|--preset xy_arrays"
  "core_tel and the four per-array options|--preset core_tel xy_arrays"
  "core with the karyotype set: every chromosome at full precision|--preset core --controls karyotype"
  "core with the screen set: every chromosome, the sex chromosomes at full precision|--preset core --controls screen"
  "core_tel with the karyotype set|--preset core_tel --controls karyotype"
  "core with every region of the bundle (the controls and the windows)|--preset core --controls all"
)

strip() {  # repository and engine paths out of the text
  sed -e "s#$ENGINE#ENGINE#g" -e "s#$ROOT/##g"
}
table() {  # plan TSV -> Markdown, the columns a reader needs
  awk -F'\t' '
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; n = split("option status tier intervals dropped capture_target expected_capture capture_lost mb_median mb_saved cum_mb_median cum_mb_floor cum_pct_median order note", k, " ")
              h = "|"; s = "|"; for (j = 1; j <= n; j++) { h = h " " k[j] " |"; s = s " --- |" } print h; print s; next }
    { r = "|"; for (j = 1; j <= n; j++) { v = $(c[k[j]]); gsub(/\|/, "\\|", v); r = r " " v " |" } print r }'
}
total() {  # the last cumulative row: MB as the engine reads it, percent of the CRAM, and the floor with every slice once
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $(c["cum_mb_median"]) != "NA" { m = $(c["cum_mb_median"]); p = $(c["cum_pct_median"]); f = $(c["cum_mb_floor"]) } END { print m "\t" p "\t" f }'
}

{
  echo "# fetchplan worked examples"
  echo
  echo "Written by \`resources/build/fetch_examples.sh\` on $(date +%Y-%m-%d); do not edit by hand."
  echo "The READMEs and \`docs/DESIGN.md\` quote plan figures (MB per genome, intervals kept, expected capture)"
  echo "from this file: rerun the script after a change to the menu, a sinks or statistics file, the control"
  echo "regions or \`ngsdose fetchplan\`. The descriptions in \`resources/fetch_menu.tsv\` and the headers and"
  echo "table of \`resources/experimental/subsets/\` carry rounded costs written when those files were made;"
  echo "examples 7 and 14 to 16 give those options' own costs (mb_median) to check them against."
  echo "Medians of cumulative totals are not additive: the difference between two rows, or between two plans'"
  echo "totals, is not the median of what an option adds per CRAM."
  echo
  echo "Costs are what \`ngs-dose count -m fetch\` reads: the CRAM slices its fetches decode, in MB (1e6 bytes), median over"
  echo "${#CRAI[@]} CRAM indexes:"
  echo "$(for c in "${CRAI[@]}"; do basename "$c" .crai; done | paste -sd, - | sed 's/,/, /g')"
  echo "(\`resources/build/fetch_examples.sh --get DIR\` downloads them; sha256 of their concatenation"
  echo "\`$(cat "${CRAI[@]}" | "${SHA[@]}" | cut -c1-16)\`)."
  echo "These are NYGC 30x CRAMs of the 1000 Genomes high-coverage release (bwa-mem, GRCh38 with decoys"
  echo "and HLA; contigs from \`$(basename "$CONTIGS")\`). The sinks, their statistics and so these costs belong to"
  echo "that aligner and reference: another pipeline learns its own sinks from its own scans. Percentages are of"
  echo "each whole CRAM. The engine joins a plan's intervals into one indexed fetch where less than 50,000 bp separate"
  echo "them (16,384 in a BAM: \`count --group-gap\`), and each fetch decodes every slice that overlaps it (with its"
  echo "container's compression header), so a slice under two fetches is decoded twice. In each table, mb_median is the"
  echo "option's own intervals alone and cum_mb_median the plan up to and including that row, both priced so; cum_mb_floor"
  echo "is the plan with every slice decoded once (the engine comes within a few percent of it: a slice is read twice only"
  echo "where it spans two fetches more than 50,000 bp apart, as the few sparse chrY slices of a woman do). Over HTTPS"
  echo "htslib moves more bytes than it decodes: each fetch opens a range request without an end, and what is in flight"
  echo "at the next seek is thrown away (measured on the windows' candidates: 1,517 MB moved for 410 MB of slices); a"
  echo "copy of the containers by exact byte range moves what is priced here. Candidate classes have no sinks yet: a plan"
  echo "lists them for the whole-file scans only."
  echo
  if "$ENGINE" count --help 2>/dev/null | grep -q -- '--classes'; then takes="takes"; else takes="does not take"; fi
  echo "Engine for count_flags.txt: \`$(basename "$ENGINE")\` ($("$ENGINE" --version 2>/dev/null || echo unknown); the binary's sha256 begins"
  echo "\`$("${SHA[@]}" < "$ENGINE" | cut -c1-16)\`, as \`--version\` does not tell builds apart), which $takes \`count --classes\`."
  echo "\`count --classes\` is in engines from the fetch-menu change of 2026-09-26 on; fae1124 lacks it. The engine changes only count_flags.txt, and whether fetchplan accepts a plan whose"
  echo "panels define classes it does not select (\`ngsdose fetchplan --help\`, --engine)."
  echo
  echo "| example | MB per genome | % of the CRAM | floor: every slice once |"
  echo "| --- | ---: | ---: | ---: |"
} > "$TMP/head.md"
: > "$TMP/body.md"

i=0
for e in "${EXAMPLES[@]}"; do
  i=$((i + 1)); name="${e%%|*}"; args="${e#*|}"
  # shellcheck disable=SC2086
  $NGSDOSE fetchplan --menu "$MENU" $args --crai "${CRAI[@]}" --contigs "$CONTIGS" --engine "$ENGINE" -o "$TMP/p$i" \
    > "$TMP/p$i.tsv" 2> "$TMP/p$i.err" || { cat "$TMP/p$i.err" >&2; echo "error: example $i ($args) failed" >&2; exit 1; }
  IFS=$'\t' read -r mb pct floor < <(total < "$TMP/p$i.tsv")
  shown="$(echo "ngsdose fetchplan $args" | strip) --crai CRAI... --contigs $(basename "$CONTIGS") --engine ENGINE"
  echo "| [$i. $name](#example-$i) | $mb | $pct | $floor |" >> "$TMP/head.md"
  {
    echo
    echo "<a id=\"example-$i\"></a>"
    echo "## $i. $name"
    echo
    echo '```'
    echo "$shown"
    echo '```'
    echo
    echo "**$mb MB per genome ($pct% of the CRAM; $floor MB with every slice decoded once).**"
    echo
    table < "$TMP/p$i.tsv" | strip
    echo
    echo "controls.txt (the -c FASTA): $(xargs -n1 basename < "$TMP/p$i.controls.txt")"
    echo
    flags="$(cat "$TMP/p$i.count_flags.txt")"
    if [ -n "$flags" ]; then echo "count_flags.txt: \`$(echo "$flags" | paste -sd' ' -)\`"
    else echo "count_flags.txt: (empty: the fetch needs no count flags)"; fi
    if ! cmp -s "$TMP/p$i.panels.txt" "$TMP/p$i.scan_panels.txt"; then
      echo
      echo "scan_panels.txt (panels for the whole-file scans): $(xargs -n1 basename < "$TMP/p$i.scan_panels.txt" | paste -sd' ' - | sed 's/ /, /g')"
    fi
    echo
    echo "Notes:"
    echo
    # one line for the candidates' identical notes, none for the 'wrote' line
    { grep -v '^\[fetchplan\] wrote ' "$TMP/p$i.err" || true; } | strip | sed 's/^\[fetchplan\] //' | awk '
      / is a candidate: no sinks have been learned for it/ { n++; if (n == 1) first = $0; next }
      { print "- " $0 }
      END { if (n) print "- (each of the " n " candidates has this note; the first:) " first
            if (NR == 0) print "- (none)" }'
  } >> "$TMP/body.md"
done
mkdir -p "$(dirname "$OUT")"
cat "$TMP/head.md" "$TMP/body.md" > "$OUT"
echo "wrote $OUT" >&2
