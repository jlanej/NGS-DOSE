#!/usr/bin/env python3
"""Re-derive the copies table (copies/copies.tsv) from the two assemblies, so that build.sh can
check the table it builds from. Only local extracts of each region are aligned (minimap2), never a
whole genome.

    locate.py GRCh38.fa CHM13.fa UNITS_DIR WORKDIR > derived.tsv

Rules (0-based half-open output: class assembly chrom start end strand):
  TDNA1Q23  the unit (asm20, -p 0.1) against GRCh38 chr1:161.40-161.52 Mb and CHM13
            chr1:160.45-160.85 Mb; the array is the first to the last unit hit.
  RNU1      the 164-bp RNU1-1 gene (sr, <= 2 mismatches) against GRCh38 chr1:15.5-17.5 Mb plus
            chr1_KI270713v1_random and CHM13 chr1:15.0-17.2 Mb; each copy is the gene -3 kb to
            +9 kb along the gene's strand (12 kb).
  SNORD3    the 4,542-bp U3 module (asm20, -p 0.05, >= 4,000 bp aligned) against GRCh38
            chr17:18.9-19.4 Mb and CHM13 chr17:18.85-19.35 Mb.
  SNORD116, SNORD115  GRCh38: the unit itself; CHM13: the unit (asm5) against chr15:22.6-23.2 Mb.
"""
import subprocess
import sys
from pathlib import Path

g38, c13, units, work = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
work.mkdir(parents=True, exist_ok=True)


def extract(fa, regions, out):
    with open(out, "w") as fh:
        subprocess.run(["samtools", "faidx", fa] + regions, stdout=fh, check=True)
    return out


def paf(target, query, preset, extra=()):
    r = subprocess.run(["minimap2", "-c", "-x", preset, *extra, str(target), str(query)],
                       capture_output=True, text=True, check=True)
    for line in r.stdout.splitlines():
        p = line.split("\t")
        tname = p[5]
        if ":" in tname and "-" in tname.rsplit(":", 1)[1]:
            chrom, iv = tname.rsplit(":", 1)
            off = int(iv.split("-")[0]) - 1
        else:
            chrom, off = tname, 0
        nm = next((int(x[5:]) for x in p[12:] if x.startswith("NM:i:")), None)
        yield dict(qs=int(p[2]), qe=int(p[3]), strand=p[4], chrom=chrom, ts=off + int(p[7]),
                   te=off + int(p[8]), match=int(p[9]), alen=int(p[10]), nm=nm)


out = []
# TDNA1Q23
for asm, fa, reg in (("GRCh38", g38, "chr1:161400001-161520000"), ("CHM13", c13, "chr1:160450001-160850000")):
    t = extract(fa, [reg], work / f"tdna.{asm}.fa")
    h = list(paf(t, units / "TDNA1Q23.fa", "asm20", ("-N", "100", "-p", "0.1", "--secondary=yes")))
    out.append(("TDNA1Q23", asm, h[0]["chrom"], min(x["ts"] for x in h), max(x["te"] for x in h), "+"))
# RNU1: the gene, from GRCh38 RNU1-1 (chr1:16514122-16514285, minus strand)
gene = work / "RNU1-1.fa"
with open(gene, "w") as fh:
    subprocess.run(["samtools", "faidx", "-i", g38, "chr1:16514122-16514285"], stdout=fh, check=True)
for asm, fa, regs in (("GRCh38", g38, ["chr1:15500001-17500000", "chr1_KI270713v1_random"]), ("CHM13", c13, ["chr1:15000001-17200000"])):
    t = extract(fa, regs, work / f"u1.{asm}.fa")
    for x in sorted(paf(t, gene, "sr", ("-N", "50", "-p", "0.3", "--secondary=yes")), key=lambda x: (x["chrom"], x["ts"])):
        if x["nm"] is not None and x["nm"] <= 2 and x["alen"] >= 160:
            if x["strand"] == "+":
                out.append(("RNU1", asm, x["chrom"], x["ts"] - 3000, x["ts"] + 9000, "+"))
            else:
                out.append(("RNU1", asm, x["chrom"], x["te"] - 9000, x["te"] + 3000, "-"))
# SNORD3
for asm, fa, reg in (("GRCh38", g38, "chr17:18900001-19400000"), ("CHM13", c13, "chr17:18850001-19350000")):
    t = extract(fa, [reg], work / f"u3.{asm}.fa")
    for x in sorted(paf(t, units / "SNORD3.fa", "asm20", ("-N", "100", "-p", "0.05", "--secondary=yes")), key=lambda x: x["ts"]):
        if x["qe"] - x["qs"] >= 4000:
            out.append(("SNORD3", asm, x["chrom"], x["ts"], x["te"], x["strand"]))
# SNORD116, SNORD115
t = extract(c13, ["chr15:22600001-23200000"], work / "snord.CHM13.fa")
for cls, reg in (("SNORD116", (25051000, 25109000)), ("SNORD115", (25170000, 25270500))):
    out.append((cls, "GRCh38", "chr15", reg[0], reg[1], "+"))
    h = [x for x in paf(t, units / f"{cls}.fa", "asm5") if x["alen"] > 0.9 * (reg[1] - reg[0])]
    out.append((cls, "CHM13", h[0]["chrom"], h[0]["ts"], h[0]["te"], h[0]["strand"]))
for r in out:
    print("\t".join(map(str, r)))
