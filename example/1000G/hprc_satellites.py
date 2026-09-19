#!/usr/bin/env python3
"""Satellite array mass: NGS-DOSE estimates from short reads against HPRC release-2 assemblies.

The HPRC CenSat annotation of a diploid assembly gives, per haplotype, the span of every
satellite array. Summed over both haplotypes that is the sample's array mass per class - an
assembly-based truth for the compositional classes of the experimental satellite panel, for the
~200 HPRC samples that are also in the 1000 Genomes 30x cohort.

    hprc_satellites.py --estimates estimates.tsv --censat DIR --out hprc_satellites.tsv

Caveats worth keeping next to the numbers: an assembly is not a truth for the arrays it failed
to span (GAP annotations are counted and reported); and the panel's k-mers come from CHM13, so a
class whose sequence differs between people is under-recovered in proportion to that.
"""
import argparse
import collections
import csv
import glob
import os

import numpy as np

# CenSat annotation name prefix -> panel class
CLASS_OF = {"HSat1A": "HSat1A", "HSat1B": "HSat1B", "HSat2": "HSat2", "HSat3": "HSat3", "bSat": "bSat", "hor": "aSatHOR", "active": "aSatHOR"}
CLASSES = ("HSat1A", "HSat1B", "HSat2", "HSat3", "bSat", "aSatHOR")


def assembly_mass(sample: str, censat_dir: str):
    """bp per class over all haplotype annotations of one sample, the number of files and of GAP records."""
    files = sorted(glob.glob(os.path.join(censat_dir, f"{sample}_*cenSat.bed")))
    mass, gaps = collections.Counter(), 0
    for path in files:
        with open(path) as fh:
            for line in fh:
                if line.startswith(("track", "#")) or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                name = p[3]
                gaps += "GAP" in name
                key = name.split("(")[0].split("_")[0]
                if key in CLASS_OF:
                    mass[CLASS_OF[key]] += int(p[2]) - int(p[1])
    return mass, len(files), gaps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--estimates", required=True, help="table from `ngsdose estimate` or `ngsdose cohort` on scan-mode counts")
    ap.add_argument("--censat", required=True, help="directory of <sample>_<hap>_hprc_r2_v1*.cenSat.bed files")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    with open(a.estimates) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mass, n_files, gaps = assembly_mass(r["sample"], a.censat)
            if n_files != 2:
                continue
            for cls in CLASSES:
                v = r.get(f"{cls}.mass_Mb", "NA")
                if v not in ("NA", "", None):
                    rows.append(dict(sample=r["sample"], cls=cls, assembly_Mb=round(mass[cls] / 1e6, 3), ngsdose_Mb=float(v), assembly_gaps=gaps))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "cls", "assembly_Mb", "ngsdose_Mb", "assembly_gaps"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = len({r["sample"] for r in rows})
    print(f"{n} samples with both haplotype annotations and satellite estimates -> {a.out}")
    print(f"{'class':8s} {'n':>4s} {'assembly Mb (median, range)':>30s} {'est/assembly (median)':>22s} {'pair SD of log ratio':>21s} {'Pearson r':>10s} {'Spearman':>9s}")
    for cls in CLASSES:
        x = np.array([r["assembly_Mb"] for r in rows if r["cls"] == cls])
        y = np.array([r["ngsdose_Mb"] for r in rows if r["cls"] == cls])
        ok = (x > 0) & (y > 0)
        x, y = x[ok], y[ok]
        if len(x) < 3:
            continue
        rank = lambda v: np.argsort(np.argsort(v))
        lr = np.log(y / x)
        print(f"{cls:8s} {len(x):4d} {np.median(x):12.1f} ({x.min():.1f}-{x.max():.1f}){'':6s} {np.exp(np.median(lr)):22.2f} {lr.std(ddof=1):21.3f} "
              f"{np.corrcoef(x, y)[0, 1]:10.3f} {np.corrcoef(rank(x), rank(y))[0, 1]:9.3f}")


if __name__ == "__main__":
    main()
