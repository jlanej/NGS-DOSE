"""Satellite array mass from HPRC release-2 assemblies, as a truth for the compositional classes.

The HPRC CenSat annotation of a diploid assembly gives, per haplotype, the span of every
satellite array. Summed over both haplotypes that is the sample's array mass per class. An
assembly is a truth only for the arrays it spans, and the annotation marks where it did not
close one in two ways: an HSat2 or HSat3 array is annotated together with its gap ("GAP,HSat2",
the record spanning the whole array), while for every other family the break is a record of its
own, labelled "GAP", next to or between the array's records. Both are tallied as gapped here: the
"GAP,X" record whole, a standalone GAP by its span, given to the class of the records it touches
(split when it touches two). A sample is left out of a class's comparison when gapped sequence
exceeds MAX_GAPPED of the class. A standalone GAP of exactly 100 bp is a placeholder of unknown
size, so its span says nothing about what it hides; such gaps are counted (`n_unsized_gap`), and
`compare(exclude_unsized=True)` leaves out every sample that has one next to the class.

The same annotation labels the rDNA the assemblies hold. They are no truth for it (the arrays are
not closed), but how much of it they hold is the answer to why the rDNA is measured from reads.
"""
from __future__ import annotations

import collections
import glob
import gzip
import os

import numpy as np

# CenSat annotation label -> panel class. A label is the category before "(" or "_" ("HSat3",
# "active_hor(...)" -> "active"); inside the catch-all "cenSat(...)" category it is the family named
# first in the parentheses ("cenSat(SST1,SST1v)" -> SST1, "cenSat(ACRO1,COMP-...)" -> ACRO1). As in
# the panel build (resources/build/build_satellite_panel.sh), a family maps to ACRO, SST1 or SATR
# by its prefix ("SATR1v" -> SATR) and to CER exactly; a composite named first after another family
# ("cenSat(HSAT5v1,SST1,SST1v)", the chrY RBMY/TSPY composites) is left out, as the build leaves it.
PREFIX_CLASSES = (("ACRO", "ACRO"), ("SST1", "SST1"), ("SATR", "SATR"))


class _ClassMap(dict):
    def __missing__(self, label):
        cls = next((c for pre, c in PREFIX_CLASSES if label.startswith(pre)), None)
        if cls is None:
            raise KeyError(label)
        return cls

    def __contains__(self, label):
        return dict.__contains__(self, label) or any(label.startswith(pre) for pre, _ in PREFIX_CLASSES)

    def get(self, label, default=None):
        return self[label] if label in self else default


CLASS_OF = _ClassMap({"HSat1A": "HSat1A", "HSat1B": "HSat1B", "HSat2": "HSat2", "HSat3": "HSat3", "bSat": "bSat", "hor": "aSatHOR",
                      "active": "aSatHOR", "ACRO1": "ACRO", "SST1": "SST1", "CER": "CER", "SATR1": "SATR", "SATR2": "SATR"})
CLASSES = ("HSat1A", "HSat1B", "HSat2", "HSat3", "bSat", "aSatHOR", "ACRO", "SST1", "CER", "SATR")
MAX_GAPPED = 0.02          # a gap matters when it could hide more than this share of what is annotated
PLACEHOLDER_GAP = 100      # the span of a gap record whose size the assembly does not know
# labels that name no panel class, as expected (not reported by `unrecognised`)
NOT_CLASSES = ("ct", "mon", "dhor", "gSat", "rDNA", "GAP", "mixedAlpha")


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


def annotation_files(sample: str, censat_dir: str) -> list[str]:
    """The sample's haplotype annotations, plain (`.cenSat.bed`) or gzipped (`.cenSat.bed.gz`), one
    file per haplotype if both forms are present."""
    found: dict[str, str] = {}
    for path in sorted(glob.glob(os.path.join(censat_dir, f"{sample}_*cenSat.bed*"))):
        if path.endswith((".bed", ".bed.gz")):
            found.setdefault(path.removesuffix(".gz"), path)
    return sorted(found.values())


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _records(path: str):
    with _open(path) as fh:
        for line in fh:
            if line.startswith(("track", "#")) or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            yield p[0], int(p[1]), int(p[2]), labels_of(p[3])


def _class_of(labels: list[str]) -> str | None:
    return next((CLASS_OF[x] for x in labels if x in CLASS_OF), None)


def _tally(sample: str, censat_dir: str, unrecognised: collections.Counter | None = None):
    """(mass, gapped, placeholder gaps touching each class, number of files)."""
    files = annotation_files(sample, censat_dir)
    mass, gapped, unsized = collections.Counter(), collections.Counter(), collections.Counter()
    for path in files:
        by_contig = collections.defaultdict(list)
        for contig, s, e, labels in _records(path):
            cls = _class_of(labels)
            if cls:
                (gapped if "GAP" in labels else mass)[cls] += e - s
            elif unrecognised is not None:
                for x in labels:
                    if x not in NOT_CLASSES and not x.startswith(("HSAT", "COMP")):
                        unrecognised[x] += e - s
            by_contig[contig].append((s, e, labels, cls))
        # a standalone GAP breaks the array it sits in or next to: its span is gapped sequence of the
        # class(es) of the records that touch it
        for recs in by_contig.values():
            for s, e, labels, _ in recs:
                if labels != ["GAP"]:
                    continue
                near = sorted({c for s2, e2, lab2, c in recs if c and lab2 != ["GAP"] and s2 <= e and e2 >= s})
                for c in near:
                    gapped[c] += (e - s) if len(near) == 1 else (e - s) / len(near)
                    if e - s == PLACEHOLDER_GAP:
                        unsized[c] += 1
    return mass, gapped, unsized, len(files)


def assembly_mass(sample: str, censat_dir: str, unrecognised: collections.Counter | None = None):
    """Per class: bp in arrays the assembly spans, and bp of gaps (a "GAP,X" array whole, a
    standalone GAP record by its span) over all haplotype annotations of one sample; and the
    number of annotation files found. `unrecognised`, when given, collects the bp of labels that
    map to no class and are not known to be other sequence (NOT_CLASSES, HSAT4/5, composites), so
    that a new spelling of a family in a later release is seen rather than dropped."""
    mass, gapped, _, n = _tally(sample, censat_dir, unrecognised)
    return mass, gapped, n


def rdna_in_assembly(sample: str, censat_dir: str, merge_bp: int = 1000):
    """Sequence annotated as rDNA in the haplotype assemblies of one sample: total bp (of the rDNA
    records themselves), the longest stretch (records on one contig less than `merge_bp` apart are
    one stretch), and the number of annotation files. Assemblies do not close the rDNA arrays, so
    this is what an assembly holds of them, not a copy number. A record annotated with a gap is not
    sequence: it is not counted, and it ends a stretch however close the next record is."""
    files = annotation_files(sample, censat_dir)
    total = longest = 0
    for path in files:
        by_contig = collections.defaultdict(list)
        for contig, s, e, labels in _records(path):
            if "GAP" in labels:
                by_contig[contig].append((s, e, True))
            elif "rDNA" in labels:
                by_contig[contig].append((s, e, False))
        for recs in by_contig.values():
            recs.sort()
            cur = None                                      # [start, end, rDNA bp] of the open stretch
            for a, b, gap in recs + [(1 << 62, 1 << 62, True)]:
                if cur and not gap and a - cur[1] < merge_bp:
                    cur[2] += max(0, b - max(a, cur[1]))
                    cur[1] = max(cur[1], b)
                    continue
                if cur:
                    total += cur[2]
                    longest = max(longest, cur[1] - cur[0])
                cur = None if gap else [a, b, b - a]
    return total, longest, len(files)


def compare(estimates: list[dict], censat_dir: str, unrecognised: collections.Counter | None = None,
            exclude_unsized: bool = False) -> tuple[list[dict], dict[str, dict]]:
    """`estimates`: rows with `sample` and `CLASS.mass_Mb` columns. Returns the per-sample rows
    (assembly, gapped and estimated Mb per class) and per-class statistics over the samples
    whose gaps are immaterial (gapped sequence at most MAX_GAPPED of the class); `n_gapped` counts
    the samples left out. `n_unsized_gap` counts the samples otherwise compared that have a
    placeholder gap (unknown size) next to the class; `exclude_unsized` leaves those out as well
    (for aSatHOR that is most samples: 172 of the 200 HPRC release-2 samples have one inside an
    array)."""
    rows, unsized = [], {}
    for r in estimates:
        mass, gapped, uns, n_files = _tally(r["sample"], censat_dir, unrecognised)
        if n_files != 2:
            continue
        for cls in CLASSES:
            v = r.get(f"{cls}.mass_Mb")
            if v in (None, "", "NA"):
                continue
            rows.append(dict(sample=r["sample"], cls=cls, assembly_Mb=round(mass[cls] / 1e6, 3), assembly_gapped_Mb=round(gapped[cls] / 1e6, 3),
                             ngsdose_Mb=float(v)))
            unsized[(r["sample"], cls)] = uns[cls]
    stats = {}
    for cls in CLASSES:
        sub = [x for x in rows if x["cls"] == cls]
        small = [x for x in sub if x["assembly_Mb"] > 0 and x["ngsdose_Mb"] > 0
                 and x["assembly_gapped_Mb"] <= MAX_GAPPED * (x["assembly_Mb"] + x["assembly_gapped_Mb"])]
        has = [x for x in small if unsized[(x["sample"], cls)]]
        clean = [x for x in small if not (exclude_unsized and unsized[(x["sample"], cls)])]
        st = dict(n=len(clean), n_gapped=len(sub) - len(clean), n_unsized_gap=len(has))
        if len(clean) >= 2:
            x = np.array([c["assembly_Mb"] for c in clean])
            y = np.array([c["ngsdose_Mb"] for c in clean])
            lr = np.log(y / x)
            med = float(np.median(lr))
            st.update(ratio_median=float(np.exp(med)), sd_log=float(lr.std(ddof=1)), assembly_median=float(np.median(x)),
                      assembly_min=float(x.min()), assembly_max=float(x.max()),
                      # per-genome agreement without the few outliers, and how much people differ: r means little when the two are alike
                      sd_log_robust=float(1.4826 * np.median(np.abs(lr - med))), cv_assembly=float(x.std(ddof=1) / x.mean()))
            # the genomes far from the family's median, and how many of them have less in the assembly than in the reads
            far = lr[np.abs(lr - med) > 3 * st["sd_log_robust"]] if st["sd_log_robust"] > 0 else lr[:0]
            st.update(n_far=int(len(far)), n_far_assembly_short=int((far > med).sum()))
            if len(clean) >= 3:
                rank = lambda v: np.argsort(np.argsort(v))
                st.update(pearson=float(np.corrcoef(x, y)[0, 1]), spearman=float(np.corrcoef(rank(x), rank(y))[0, 1]))
        stats[cls] = st
    return rows, stats
