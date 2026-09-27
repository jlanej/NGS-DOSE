#!/usr/bin/env python3
"""Assemble the source copies each class is scored against.

A class is built from one or two reference strains, so recall has to be measured on the other
genomes of the same organism - including, for HHV-6, the integrated genomes sequenced from this
cohort's own carriers. The EBV type-2 class is scored on the homologous block of a second type-2
genome (Jijoye), located by aligning the AG876 regions to it, and its specificity is scored on
type-1 genomes that were not in the background.

    mkrecall.py SRCDIR OUTDIR
"""
import subprocess
import sys
from pathlib import Path

SRC, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

SETS = {
    # class: (label, accession) - the first is the strain the class was built from
    "HHV6A": [("U1102_source", "NC_001664.4"), ("GS", "KC465951.1"), ("AJ", "KP257584.1"),
              ("ici_HG00657", "MG894370.1"), ("ici_NA18999", "MG894374.1"),
              ("NA18999_KY316047", "KY316047.1"), ("ici_3A_10q26", "KY316049.1")],
    "HHV6B": [("Z29_source", "NC_000898.1"), ("HST", "AB021506.1"), ("ici_HG00245", "MG894368.1"),
              ("ici_HG00362", "MG894369.1"), ("ici_HG02016", "MG894372.1"), ("ici_HG02301", "MG894373.1"),
              ("ici_GTEX_YF7O", "MH698401.1"), ("ici_GTEX_OXRO", "MH698402.1"), ("ici_2B_9q34", "KY316045.1")],
    "HHV7": [("RK_source", "NC_001716.2"), ("JI", "U43400.1")],
    "SMRV": [("HLB_source", "NC_001514.1"), ("SMRV_H", "M23385.1"), ("MK561030", "MK561030.1")],
    "PHIX": [("NC_001422_source", "NC_001422.1"), ("J02482", "J02482.1")],
    "MYCO": [("hyorhinis_HUB1_source", "NC_014448.1"), ("hyorhinis_IMT49388", "NZ_CP064323.1"),
             ("hyorhinis_N36N21c", "NZ_CP102738.1"), ("orale_NCTC10112_source", "NZ_LR214940.1"),
             ("hominis_ATCC23114_source", "NC_013511.1"), ("hominis_SP2565", "NZ_CP055144.1"),
             ("hominis_2539", "NZ_CP026341.1"), ("fermentans_M64_source", "NC_014921.1"),
             ("fermentans_PG18", "NC_021002.1"), ("fermentans_JER", "NC_014552.1"),
             ("arginini_HAZ145_source", "NZ_AP014657.1"), ("arginini_NCTC10129", "NZ_CP143577.1"),
             ("arginini_LGL1", "NZ_CP165990.1")],
}
# EBV type 1 genomes that were NOT in the panel's background: a false-positive check for EBV2
SPEC_EBV1 = [("C666_1", "OR134767.1"), ("NPC43", "OR134769.1"), ("C17", "OR134766.1"),
             ("sLCL_IM1_16", "LN827799.1"), ("chrEBV_B95_8_Raji", "NC_007605.1")]


def read(acc):
    return "".join(l.strip() for l in (SRC / f"{acc}.fa").read_text().splitlines() if not l.startswith(">")).upper()


def write(path, recs):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


for cls, items in SETS.items():
    write(OUT / f"recall_{cls}.fa", [(f"{cls}_{lab}_{acc}", read(acc)) for lab, acc in items])

# EBV2: the AG876 source regions, plus the same block of Jijoye located by alignment
ag = OUT.parent / "units" / "EBV2.fa"
jij = SRC / "LN827800.1.fa"
paf = subprocess.run(["minimap2", "-x", "asm20", "-c", "--secondary=no", str(jij), str(ag)],
                     capture_output=True, text=True).stdout
recs = []
for rec in ag.read_text().split(">")[1:]:
    lines = rec.split("\n")
    recs.append((lines[0].split()[0], "".join(lines[1:])))
jseq = read("LN827800.1")
for line in paf.splitlines():
    f = line.split("\t")
    q, ts, te = f[0], int(f[7]), int(f[8])
    if te - ts >= 500:
        recs.append((f"EBV2_Jijoye_LN827800.1_{q.split('_')[3]}_{ts}_{te}", jseq[ts:te]))
write(OUT / "recall_EBV2.fa", recs)
write(OUT / "spec_EBV1.fa", [(f"EBV1_{lab}_{acc}", read(acc)) for lab, acc in SPEC_EBV1])

for f in sorted(OUT.glob("recall_*.fa")) + [OUT / "spec_EBV1.fa"]:
    n = sum(1 for l in f.read_text().splitlines() if l.startswith(">"))
    bp = sum(len(l) for l in f.read_text().splitlines() if not l.startswith(">"))
    print(f"{f.name}\t{n} records\t{bp} bp")
