#!/usr/bin/env python3
"""recall.py PANEL [PANEL_GRCH38_ONLY]

Per class and per source allele (the arrays in sources/loci.*.tsv): 150-bp reads tiled every 10 bp
over the array, share carrying >= 4 of the class's k-mers (the engine's rule), and median hits per
read; with the full panel and, if given, with the panel built from the GRCh38 primary allele alone
(for the other alleles that is a held-out test). eff_len_bp: read start positions per strand, over the
array +-5 kb tiled every 1 bp, whose read carries >= 4 class k-mers - what the compositional mass
counts for that allele (array length plus the edge term edge_bp; negative where the panel is blind to
part of the array). flank_reads_with_hits: of those, reads wholly outside the array (partial units
just beyond the TRF span, inside the +-200 bp mask). No class k-mer occurs in either assembly outside
the masked arrays, so no read from elsewhere in the two references is assigned."""
import gzip
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMP = str.maketrans("ACGT", "TGCA")
R, K, MIN = 150, 31, 4


def canonical(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def load(path):
    names, km = {}, {}
    for line in gzip.open(path, "rt"):
        if line.startswith("##class"):
            f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
            names[f["id"]] = f["name"]
        elif not line.startswith("#"):
            p = line.split("\t")
            km[p[0]] = names[p[1]]
    return km


def faidx(fa, reg):
    out = subprocess.run(["samtools", "faidx", fa, reg], check=True, capture_output=True, text=True).stdout
    return "".join(out.split("\n")[1:]).upper()


def hits_per_start(seq, cls, km):
    """hits[i] = class k-mers in seq[i:i+R] (both strands give the same count: k-mers are canonical)"""
    isk = [1 if km.get(canonical(seq[j:j + K])) == cls else 0 for j in range(len(seq) - K + 1)]
    cs = [0]
    for v in isk:
        cs.append(cs[-1] + v)
    n = R - K + 1
    return [cs[i + n] - cs[i] for i in range(len(seq) - R + 1)]


def main():
    full = load(sys.argv[1])
    loo = load(sys.argv[2]) if len(sys.argv) > 2 else None
    g38, chm = os.environ["GRCH38"], os.environ["CHM13"]
    alleles = []
    for line in open(HERE / "sources/loci.GRCh38.tsv"):
        g, c, s, e = line.split("\t")[:4]
        alleles.append((g, "GRCh38", g38, c, int(s), int(e)))
    for line in open(HERE / "sources/loci.GRCh38_alts.tsv"):
        g, c, s, e = line.split("\t")[:4]
        alleles.append((g, "GRCh38_alt", g38, c, int(s), int(e)))
    for line in open(HERE / "sources/loci.CHM13.tsv"):
        g, c, s, e = line.split("\t")[:4]
        alleles.append((g, "CHM13", chm, c, int(s), int(e)))
    alleles.sort(key=lambda a: (a[0], a[1] != "GRCh38", a[1]))
    print("class\tallele\tarray_bp\treads\trecall\tmedian_hits\trecall_grch38only_panel\tmedian_hits_grch38only\teff_len_bp\tedge_bp\tflank_reads_with_hits")
    for g, src, fa, c, s, e in alleles:
        cls = "VNTR_" + g
        L = e - s
        wide = faidx(fa, f"{c}:{max(1, s - 5000 + 1)}-{e + 5000}")
        off = s - max(0, s - 5000)
        arr = wide[off:off + L]
        row = [cls, f"{src}:{c}:{s}-{e}", str(L)]
        for panel in ([full, loo] if loo is not None else [full]):
            if L >= R:
                h = hits_per_start(arr, cls, panel)[::10]
                rec = sum(x >= MIN for x in h) / len(h)
                med = sorted(h)[len(h) // 2]
                if panel is full:
                    row += [str(len(h))]
                row += [f"{rec:.3f}", str(med)]
            else:
                if panel is full:
                    row += ["0"]
                row += ["NA", "NA"]
        if loo is None:
            row += ["NA", "NA"]
        # effective length: every read start in the array +-5 kb whose read is assigned to the class
        h = hits_per_start(wide, cls, full)
        eff = sum(x >= MIN for x in h)
        # of those, reads that start more than 150 bp before the array or after its end
        fl = sum(x >= MIN for i, x in enumerate(h) if i < off - R or i >= off + L)
        row += [str(eff), str(eff - L), str(fl)]
        print("\t".join(row))


if __name__ == "__main__":
    main()
