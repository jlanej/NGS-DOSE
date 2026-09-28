#!/usr/bin/env bash
# Recall of panel.tsv.gz per source copy: 150-bp reads tiled every 10 bp over every copy of each
# class in GRCh38 (primary and alts) and CHM13 - the mask intervals, merged - and the share with
# >= 4 of the class's k-mers. Then leakage: the same over paralogs that are background for a class
# (the share of a paralog's reads the class would wrongly claim). Output recall.tsv, leak.tsv.
# Called by build.sh with RECALL=1; needs bedtools and samtools.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="${ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
WORK="${WORK:-$ROOT/work/ref}"; PY="${PY:-python3}"
OUT="${OUT:-$ROOT/work/candidates/multicopy-genes}"; MASKS="${MASKS:-$HERE/masks}"
G38="$WORK/GRCh38_full_analysis_set_plus_decoy_hla.fa"; CHM="$WORK/chm13v2.0.fa"
cd "$OUT"; mkdir -p recall
classes=$(awk -F'\t' '!/^#/ && $7 == "class" {print $1}' "$HERE/classes.tsv")
args=()
for c in $classes; do
  for asm in GRCh38 CHM13; do
    fa=$G38; [ $asm = CHM13 ] && fa=$CHM
    [ -s "$MASKS/$c.$asm.bed" ] || continue
    sort -k1,1 -k2,2n "$MASKS/$c.$asm.bed" | bedtools merge -i - | awk '{printf "%s:%d-%d\n", $1, $2 + 1, $3}' > recall/$c.$asm.regions
    samtools faidx "$fa" -r recall/$c.$asm.regions > recall/$c.$asm.fa
    args+=("$c=$asm=recall/$c.$asm.fa")
  done
done
"$PY" "$HERE/recall.py" panel.tsv.gz "${args[@]}" > recall.tsv

# leakage sources: paralog rows of classes.tsv, plus a few loci named here
leak=()
pl() { leak+=("$1=$2=recall/leak.$2.fa"); }
unit_leak() { cp aux/$2.fa recall/leak.$2.fa; pl "$1" "$2"; }
reg_leak() { local c=$1 name=$2 fa=$3; shift 3; samtools faidx "$fa" "$@" > recall/leak.$name.fa; pl "$c" "$name"; }
unit_leak KIV2 KIVflank5; unit_leak KIV2 KIVflank3
unit_leak SMN1 SMN2g
reg_leak SMN1 SMN2like_alt_CHM13 "$G38" chr5_KI270897v1_alt:274862-302933
samtools faidx "$CHM" chr5:70809747-70837821 >> recall/leak.SMN2like_alt_CHM13.fa
unit_leak CYP21A2 CYP21A1P
unit_leak RHD RHCE
unit_leak HPR HPgene; unit_leak HPdup HPRgene
unit_leak ORM1 ORM2
unit_leak AMY1 AMY2Ap; cp aux/AMY2A.fa recall/leak.AMY2Afull.fa; pl AMY1 AMY2Afull
unit_leak CCL3L CCL3blk
unit_leak UGT2B17 UGT2B15
unit_leak GSTM1 GSTM2to5
reg_leak HERVC4 HERVK_elsewhere "$G38" chrY:12995015-13001093 chr6:27187544-27194236 chr10:39107643-39112352 chr6:122504854-122512077 chr19:52461504-52466478 chr1:19926927-19932710
reg_leak HBA HBAP1_HBM_HBQ1 "$G38" chr16:169368-170788 chr16:165977-166764 chr16:180458-181179
reg_leak RASA4 RASA4_paralogs "$G38" chr7:43961331-44041894 chr7:102341679-102357966
reg_leak TBC1D3 TBC1D3P1_P2 "$G38" chr17:60008125-60019098 chr17:62264705-62275701
reg_leak LCE3BC LCE_other "$G38" chr1:152570000-152583066 chr1:152615264-152640000
reg_leak APOBEC3B APOBEC3A_side "$G38" chr22:38940000-38961989
"$PY" "$HERE/recall.py" panel.tsv.gz "${leak[@]}" > leak.tsv

# C4 and HERVC4 on their own parts of each C4 copy: the insertion of each long copy (HERVC4), and
# each copy without it (C4)
for asm in GRCh38 CHM13; do
  "$PY" "$HERE/subseq.py" inside units/HERVC4.fa recall/C4.$asm.fa > recall/HERVC4.$asm.insert.fa
  "$PY" "$HERE/subseq.py" outside units/HERVC4.fa recall/C4.$asm.fa > recall/C4.$asm.noinsert.fa
done
"$PY" "$HERE/recall.py" panel.tsv.gz C4=GRCh38=recall/C4.GRCh38.noinsert.fa C4=CHM13=recall/C4.CHM13.noinsert.fa \
  HERVC4=GRCh38=recall/HERVC4.GRCh38.insert.fa HERVC4=CHM13=recall/HERVC4.CHM13.insert.fa > recall.c4parts.tsv

# SMN1: reads that cover at least one SMN1/SMN2 site of their copy (sites of each SMN1-like copy
# against GRCh38 SMN2), and reads that cover c.840
"$PY" - "$HERE" <<'PYEOF' > recall.smn1sites.tsv
import subprocess, sys
sys.path.insert(0, sys.argv[1])
from recall import load_panel, canonical, read_fasta
from diffsites import sites as diffsites
k, kmers = load_panel("panel.tsv.gz")
print("copy\treads_over_a_site\trecall\treads_over_c840\trecall_c840")
for asm in ("GRCh38", "CHM13"):
    for name, seq in read_fasta(f"recall/SMN1.{asm}.fa"):
        open("recall/one.fa", "w").write(f">{name}\n{seq}\n")
        st = sorted(diffsites("recall/one.fa", "aux/SMN2g.fa"))
        # c.840: the site whose 20-bp context matches the GRCh38 SMN1 exon 7 context
        c840 = [x for x in st if seq[x - 10:x + 11] == "TACAGGGTTTCAGACAAAATC"]
        n = h = n8 = h8 = 0
        for i in range(0, len(seq) - 150 + 1, 10):
            hits = sum(1 for j in range(i, i + 150 - k + 1) if kmers.get(canonical(seq[j:j + k])) == "SMN1")
            if any(i <= x < i + 150 for x in st):
                n += 1; h += hits >= 4
            if any(i <= x < i + 150 for x in c840):
                n8 += 1; h8 += hits >= 4
        print(f"{asm}:{name}\t{n}\t{h / n:.3f}\t{n8}\t{(h8 / n8 if n8 else float('nan')):.3f}")
PYEOF
