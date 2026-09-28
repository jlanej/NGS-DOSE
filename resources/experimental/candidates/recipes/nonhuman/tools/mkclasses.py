#!/usr/bin/env python3
"""Cut the class units from the downloaded GenBank records.

Positional units of the roseoloviruses keep the direct repeat once (DRL + U, the right-hand
copy removed), so that every DR k-mer occurs once in the unit and the engine places it: a
positional class drops k-mers seen twice in its unit, which would otherwise delete the whole DR.
The Mycoplasma genomes have their rRNA and tRNA features replaced by N, because those operons
are shared with bacteria that are not in any background list.
"""
import re
import sys
from pathlib import Path

SRC = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)


def read(acc):
    seq = []
    for line in (SRC / f"{acc}.fa").read_text().splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def write(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def rna_intervals(acc, pad=50):
    """rRNA/tRNA/ncRNA feature intervals of a GenBank feature table, 0-based half-open."""
    iv = []
    for line in (SRC / f"{acc}.ft").read_text().splitlines():
        m = re.match(r"^[<>]?(\d+)\t[<>]?(\d+)\t(rRNA|tRNA|ncRNA|misc_RNA|tmRNA)\s*$", line)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            s, e = min(a, b) - 1, max(a, b)
            iv.append((max(0, s - pad), e + pad))
    return iv


def mask(seq, iv):
    b = bytearray(seq, "ascii")
    n = 0
    for s, e in iv:
        e = min(e, len(b))
        n += e - s
        b[s:e] = b"N" * (e - s)
    return b.decode(), n


# --- roseoloviruses: unit = the genome with the right-hand terminal direct repeat removed
# NC_001664.4 (HHV-6A U1102): repeat_region 1-8089 DRL, 151289-159378 DRR
# NC_000898.1 (HHV-6B Z29):   repeat_region 1-8793 DRL, 153322-162114 DRR
# NC_001716.2 (HHV-7 RK):     repeat_region 1-10034 DR,  143047-153080 DR
for cls, acc, end in [("HHV6A", "NC_001664.4", 151288), ("HHV6B", "NC_000898.1", 153321), ("HHV7", "NC_001716.2", 143046)]:
    s = read(acc)
    assert len(s) in (159378, 162114, 153080), (acc, len(s))
    write(OUT / f"{cls}.fa", [(f"{cls}_{acc}_1_{end}_DRL_plus_U", s[:end])])

# --- SMRV: the whole proviral genome
write(OUT / "SMRV.fa", [("SMRV_NC_001514.1", read("NC_001514.1"))])

# --- phiX174: the whole circular genome
write(OUT / "PHIX.fa", [("PHIX_NC_001422.1", read("NC_001422.1"))])

# --- EBV type 2: AG876 EBNA-2 and the EBNA-3A/3B/3C block (the type-divergent latency genes)
ag = read("NC_009334.1")
write(OUT / "EBV2.fa", [("EBV2_NC_009334.1_EBNA2_36201_37565", ag[36200:37565]),
                        ("EBV2_NC_009334.1_EBNA3A_3C_80026_89937", ag[80025:89937])])

# --- Mycoplasma: five culture species, rRNA/tRNA operons masked
myco = [("Mesomycoplasma_hyorhinis", "NC_014448.1"), ("Metamycoplasma_orale", "NZ_LR214940.1"),
        ("Metamycoplasma_hominis", "NC_013511.1"), ("Mycoplasmopsis_fermentans", "NC_014921.1"),
        ("Mycoplasmopsis_arginini", "NZ_AP014657.1")]
recs, log = [], []
for name, acc in myco:
    s = read(acc)
    iv = rna_intervals(acc)
    s, n = mask(s, iv)
    recs.append((f"{name}_{acc}", s))
    log.append(f"{name}\t{acc}\t{len(s)}\t{len(iv)}\t{n}")
write(OUT / "MYCO.fa", recs)
(OUT / "MYCO.rna_masked.tsv").write_text("#species\taccession\tbp\trna_features\tbp_masked\n" + "\n".join(log) + "\n")
print("\n".join(log))
for f in sorted(OUT.glob("*.fa")):
    n = sum(len(x) for x in f.read_text().split("\n") if not x.startswith(">"))
    print(f"{f.name}\t{n}")
