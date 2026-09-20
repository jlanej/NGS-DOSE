#!/usr/bin/env python3
"""Satellite array mass: NGS-DOSE estimates from short reads against HPRC release-2 assemblies.

The HPRC CenSat annotation of a diploid assembly gives, per haplotype, the span of every
satellite array. Summed over both haplotypes that is the sample's array mass per class - an
assembly-based truth for the compositional classes of the experimental satellite panel, for the
~200 HPRC samples that are also in the 1000 Genomes 30x cohort.

    hprc_satellites.py --estimates estimates.tsv --censat DIR --out hprc_satellites.tsv

Caveats worth keeping next to the numbers: an assembly is not a truth for an array it failed to
span - such arrays are annotated together with their gap ("GAP,HSat2"), are tallied separately,
and a sample that has any in a class is left out of that class's comparison; and the panel's
k-mers come from CHM13, so a class whose sequence differs between people, or too few of whose
reads carry the four k-mers a read needs, is under-recovered in proportion (the recall of each
class on CHM13 itself is in resources/experimental/README.md).
"""
import argparse
import collections
import csv
import glob
import os

import numpy as np

# CenSat annotation label -> panel class. A label is the category before "(" or "_" ("HSat3",
# "active_hor(...)" -> "active"); inside the catch-all "cenSat(...)" category it is the family named
# first in the parentheses ("cenSat(SST1,SST1v)" -> SST1, "cenSat(ACRO1,COMP-...)" -> ACRO1).
CLASS_OF = {"HSat1A": "HSat1A", "HSat1B": "HSat1B", "HSat2": "HSat2", "HSat3": "HSat3", "bSat": "bSat", "hor": "aSatHOR", "active": "aSatHOR",
            "ACRO1": "ACRO", "SST1": "SST1", "CER": "CER", "SATR1": "SATR", "SATR2": "SATR"}
CLASSES = ("HSat1A", "HSat1B", "HSat2", "HSat3", "bSat", "aSatHOR", "ACRO", "SST1", "CER", "SATR")
MAX_GAPPED = 0.02


def labels_of(name: str) -> list[str]:
    """The labels of one annotation record. An array that the assembly did not span is labelled
    together with its gap ("GAP,HSat3"): the record is that class, of unknown true size."""
    out = []
    depth, part = 0, ""
    for ch in name + ",":                                   # split on commas outside parentheses
        if ch == "," and depth == 0:
            out.append(part)
            part = ""
        else:
            depth += ch == "("
            depth -= ch == ")"
            part += ch
    labels = []
    for part in out:
        cat = part.split("(")[0].split("_")[0]
        if cat == "cenSat" and "(" in part:
            cat = part.split("(", 1)[1].split(",")[0].rstrip(")")
        labels.append(cat)
    return labels


def assembly_mass(sample: str, censat_dir: str):
    """Per class: bp in arrays the assembly spans, and bp in arrays annotated with a gap (a lower
    bound on something larger), over all haplotype annotations of one sample; and the number of files."""
    files = sorted(glob.glob(os.path.join(censat_dir, f"{sample}_*cenSat.bed")))
    mass, gapped = collections.Counter(), collections.Counter()
    for path in files:
        with open(path) as fh:
            for line in fh:
                if line.startswith(("track", "#")) or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                labels = labels_of(p[3])
                cls = next((CLASS_OF[x] for x in labels if x in CLASS_OF), None)
                if cls:
                    (gapped if "GAP" in labels else mass)[cls] += int(p[2]) - int(p[1])
    return mass, gapped, len(files)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--estimates", required=True, help="table from `ngsdose estimate` or `ngsdose cohort` on scan-mode counts")
    ap.add_argument("--censat", required=True, help="directory of <sample>_<hap>_hprc_r2_v1*.cenSat.bed files")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    with open(a.estimates) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mass, gapped, n_files = assembly_mass(r["sample"], a.censat)
            if n_files != 2:
                continue
            for cls in CLASSES:
                v = r.get(f"{cls}.mass_Mb", "NA")
                if v not in ("NA", "", None):
                    rows.append(dict(sample=r["sample"], cls=cls, assembly_Mb=round(mass[cls] / 1e6, 3),
                                     assembly_gapped_Mb=round(gapped[cls] / 1e6, 3), ngsdose_Mb=float(v)))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "cls", "assembly_Mb", "assembly_gapped_Mb", "ngsdose_Mb"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = len({r["sample"] for r in rows})
    print(f"{n} samples with both haplotype annotations and satellite estimates -> {a.out}")
    print(f"An assembly is a truth only for the arrays it spans: a sample with more than {100 * MAX_GAPPED:.0f}% of a class's annotated sequence in "
          "gap-containing arrays is left out of that class's comparison (column 'gapped' counts them).")
    print(f"{'class':8s} {'n':>4s} {'gapped':>7s} {'assembly Mb (median, range)':>30s} {'est/assembly (median)':>22s} {'SD of log ratio':>16s} {'Pearson r':>10s} {'Spearman':>9s}")
    for cls in CLASSES:
        sub = [r for r in rows if r["cls"] == cls]
        # a gap matters when it could hide a material part of the class: more than 2% of what is annotated
        clean = [r for r in sub if r["assembly_Mb"] > 0 and r["ngsdose_Mb"] > 0
                 and r["assembly_gapped_Mb"] <= MAX_GAPPED * (r["assembly_Mb"] + r["assembly_gapped_Mb"])]
        x = np.array([r["assembly_Mb"] for r in clean])
        y = np.array([r["ngsdose_Mb"] for r in clean])
        if len(x) < 2:
            print(f"{cls:8s} {len(x):4d} {len(sub) - len(clean):7d}   too few samples without gaps")
            continue
        rank = lambda v: np.argsort(np.argsort(v))
        lr = np.log(y / x)
        corr = (f"{np.corrcoef(x, y)[0, 1]:10.3f} {np.corrcoef(rank(x), rank(y))[0, 1]:9.3f}" if len(x) >= 3 else f"{'':>10s} {'':>9s}")
        print(f"{cls:8s} {len(x):4d} {len(sub) - len(clean):7d} {np.median(x):12.1f} ({x.min():.1f}-{x.max():.1f}){'':6s} {np.exp(np.median(lr)):22.2f} {lr.std(ddof=1):16.3f} {corr}")


if __name__ == "__main__":
    main()
