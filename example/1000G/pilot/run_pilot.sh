#!/usr/bin/env bash
# The pilot, exactly as run: four 1000 Genomes trios (NYGC 30x, AWS Open Data mirror) and an
# independent, older library of every sample (HGSVC for the nine non-CEU samples, Illumina Platinum
# for the CEU trio; EBI), counted in fetch mode straight from the public CRAMs. One to two minutes
# per NYGC sample and three to four per 75x HGSVC sample on a home connection; nothing but the
# .crai index files is downloaded.
#
#   REF=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa bash run_pilot.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"
BIN="${NGSDOSE_BIN:-$REPO/target/release/ngs-dose}"; B="$REPO/resources/GRCh38"
REF="${REF:?set REF to the GRCh38 analysis-set FASTA (CRAM decoding)}"
IDX="${IDX:-$HERE/crai}"; mkdir -p "$IDX" "$HERE/counts_nygc" "$HERE/counts_replicates"
NYGC=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage
HGSVC=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data

count() {  # sample url outdir
  local out="$3/$1.json.gz" crai="$IDX/$(basename "$2").crai" attempt
  [ -s "$out" ] && return 0
  [ -s "$crai" ] || { curl -sSL --fail --retry 5 -o "$crai.part" "$2.crai" && mv "$crai.part" "$crai"; } || { echo "FAILED $1 (index)" >&2; return 1; }
  # fetch mode is bound by request latency, not bandwidth (EBI can take seconds to answer a range
  # request), so many connections; a connection that goes silent ends the engine with status 75
  for attempt in 1 2 3 4; do
    "$BIN" count -m fetch -@ 32 --stall-timeout 120 -i "$2" --index "$crai" -T "$REF" -s "$1" \
        -p "$B/panel.k31.tsv.gz" -c "$B/controls.fa.gz" --sinks "$B/sinks.bed" -o "$out" && return 0
    echo "attempt $attempt failed: $1" >&2; sleep 20
  done
  echo "FAILED $1" >&2; return 1
}

while read -r s path; do count "$s" "$NYGC/$path" "$HERE/counts_nygc"; done <<'LIST'
NA12878 data/ERR3239334/NA12878.final.cram
NA12891 additional_698_related/data/ERR3989341/NA12891.final.cram
NA12892 additional_698_related/data/ERR3989342/NA12892.final.cram
NA19240 additional_698_related/data/ERR3989410/NA19240.final.cram
NA19238 data/ERR3239453/NA19238.final.cram
NA19239 data/ERR3239454/NA19239.final.cram
HG00514 additional_698_related/data/ERR3988781/HG00514.final.cram
HG00512 additional_698_related/data/ERR3988780/HG00512.final.cram
HG00513 data/ERR3241684/HG00513.final.cram
HG00733 additional_698_related/data/ERR3988823/HG00733.final.cram
HG00731 data/ERR3241754/HG00731.final.cram
HG00732 data/ERR3241755/HG00732.final.cram
LIST

while read -r s pop; do
  count "$s" "$HGSVC/$pop/$s/high_cov_alignment/$s.alt_bwamem_GRCh38DH.20150715.$pop.high_coverage.cram" "$HERE/counts_replicates"
done <<'LIST'
HG00512 CHS
HG00513 CHS
HG00514 CHS
HG00731 PUR
HG00732 PUR
HG00733 PUR
NA19238 YRI
NA19239 YRI
NA19240 YRI
LIST

# a third library generation for the CEU trio: Illumina Platinum pedigree (HiSeq 2000, 2x100, ~50x)
PLAT=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/illumina_platinum_pedigree/data/CEU
for s in NA12878 NA12891 NA12892; do
  count "$s" "$PLAT/$s/alignment/$s.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram" "$HERE/counts_replicates"
done

# Hall, Turner & Queitsch 2021, Supplementary Data 1: per-sample estimates from the same NYGC CRAMs
[ -s "$HERE/hall2021_MOESM1.txt" ] || curl -sSL -o "$HERE/hall2021_MOESM1.txt" \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-80049-y/MediaObjects/41598_2020_80049_MOESM1_ESM.txt"
python3 "$HERE/evaluate_pilot.py"
