"""Satellite array mass from HPRC release-2 assemblies, as a truth for the compositional classes.

The HPRC CenSat annotation of a diploid assembly gives, per haplotype, the span of every
satellite array. Summed over both haplotypes that is the sample's array mass per class. An
assembly is a truth only for the arrays it spans: an array the assembly did not close is
annotated together with its gap ("GAP,HSat2"), is a lower bound on something larger, and is
tallied separately here - a sample with a material share of a class in such arrays is left out
of that class's comparison.

The same annotation labels the rDNA the assemblies hold. They are no truth for it (the arrays are
not closed), but how much of it they hold is the answer to why the rDNA is measured from reads.
"""
from __future__ import annotations

import collections
import glob
import os

import numpy as np

# CenSat annotation label -> panel class. A label is the category before "(" or "_" ("HSat3",
# "active_hor(...)" -> "active"); inside the catch-all "cenSat(...)" category it is the family named
# first in the parentheses ("cenSat(SST1,SST1v)" -> SST1, "cenSat(ACRO1,COMP-...)" -> ACRO1).
CLASS_OF = {"HSat1A": "HSat1A", "HSat1B": "HSat1B", "HSat2": "HSat2", "HSat3": "HSat3", "bSat": "bSat", "hor": "aSatHOR", "active": "aSatHOR",
            "ACRO1": "ACRO", "SST1": "SST1", "CER": "CER", "SATR1": "SATR", "SATR2": "SATR"}
CLASSES = ("HSat1A", "HSat1B", "HSat2", "HSat3", "bSat", "aSatHOR", "ACRO", "SST1", "CER", "SATR")
MAX_GAPPED = 0.02          # a gap matters when it could hide more than this share of what is annotated


def labels_of(name: str) -> list[str]:
    """The labels of one annotation record ("GAP,HSat3" -> ["GAP", "HSat3"])."""
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
    """Per class: bp in arrays the assembly spans, and bp in arrays annotated with a gap, over
    all haplotype annotations of one sample; and the number of annotation files found."""
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


def rdna_in_assembly(sample: str, censat_dir: str, merge_bp: int = 1000):
    """Sequence annotated as rDNA in the haplotype assemblies of one sample: total bp, the longest
    stretch (records on one contig less than `merge_bp` apart are one stretch), and the number of
    annotation files. Assemblies do not close the rDNA arrays, so this is what an assembly holds of
    them, not a copy number; records annotated with a gap are not sequence and are not counted."""
    files = sorted(glob.glob(os.path.join(censat_dir, f"{sample}_*cenSat.bed")))
    total = longest = 0
    for path in files:
        by_contig = collections.defaultdict(list)
        with open(path) as fh:
            for line in fh:
                if line.startswith(("track", "#")) or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                labels = labels_of(p[3])
                if "rDNA" in labels and "GAP" not in labels:
                    by_contig[p[0]].append((int(p[1]), int(p[2])))
        for ivs in by_contig.values():
            ivs.sort()
            s, e = ivs[0]
            for a, b in ivs[1:] + [(1 << 62, 1 << 62)]:
                if a - e < merge_bp:
                    e = max(e, b)
                else:
                    total += e - s
                    longest = max(longest, e - s)
                    s, e = a, b
    return total, longest, len(files)


def compare(estimates: list[dict], censat_dir: str) -> tuple[list[dict], dict[str, dict]]:
    """`estimates`: rows with `sample` and `CLASS.mass_Mb` columns. Returns the per-sample rows
    (assembly, gapped and estimated Mb per class) and per-class statistics over the samples
    whose gaps are immaterial."""
    rows = []
    for r in estimates:
        mass, gapped, n_files = assembly_mass(r["sample"], censat_dir)
        if n_files != 2:
            continue
        for cls in CLASSES:
            v = r.get(f"{cls}.mass_Mb")
            if v in (None, "", "NA"):
                continue
            rows.append(dict(sample=r["sample"], cls=cls, assembly_Mb=round(mass[cls] / 1e6, 3), assembly_gapped_Mb=round(gapped[cls] / 1e6, 3),
                             ngsdose_Mb=float(v)))
    stats = {}
    for cls in CLASSES:
        sub = [x for x in rows if x["cls"] == cls]
        clean = [x for x in sub if x["assembly_Mb"] > 0 and x["ngsdose_Mb"] > 0
                 and x["assembly_gapped_Mb"] <= MAX_GAPPED * (x["assembly_Mb"] + x["assembly_gapped_Mb"])]
        st = dict(n=len(clean), n_gapped=len(sub) - len(clean))
        if len(clean) >= 2:
            x = np.array([c["assembly_Mb"] for c in clean])
            y = np.array([c["ngsdose_Mb"] for c in clean])
            lr = np.log(y / x)
            st.update(ratio_median=float(np.exp(np.median(lr))), sd_log=float(lr.std(ddof=1)), assembly_median=float(np.median(x)),
                      assembly_min=float(x.min()), assembly_max=float(x.max()))
            if len(clean) >= 3:
                rank = lambda v: np.argsort(np.argsort(v))
                st.update(pearson=float(np.corrcoef(x, y)[0, 1]), spearman=float(np.corrcoef(rank(x), rank(y))[0, 1]))
        stats[cls] = st
    return rows, stats
