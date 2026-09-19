"""Learn class sinks from scan-mode runs.

A multi-copy class has no fixed address in a linear reference: the aligner scatters its reads
over whatever paralogs, decoys and look-alike loci the reference happens to contain. Scan mode
classifies every read by k-mers regardless of placement and records where class reads *were*
placed; the sinks are the (few) intervals that hold essentially all of them. Fetch mode then
retrieves only those intervals. Sinks depend on the reference build and on the aligner, so they
are re-learned, not assumed, whenever either changes.
"""
from __future__ import annotations

from collections import defaultdict


def learn(scan_counts, min_frac: float = 1e-5, pad: int = 1000, min_reads: int = 25) -> tuple[list[tuple[str, int, int, str]], dict]:
    """Intervals holding >= min_frac of a class's reads (and >= min_reads reads, so that a
    small class does not collect every stray placement) in any sample, merged and padded.

    `scan_counts` is an iterable of counts dicts or of paths; files are read one at a time, so a
    whole cohort of scans can be passed. Only positional classes get sinks: a satellite family's
    reads are spread over too much of the alignment for targeted retrieval to make sense.
    Returns (BED rows, per-class capture: the fraction of each sample's scan-mode class reads
    that fall inside the learned sinks)."""
    from .io import load_counts
    keep: dict[str, set[tuple[str, int]]] = defaultdict(set)
    summaries = []                                         # per sample: {class: {(contig, start): reads}} is too big; keep totals + kept bins
    width = None
    for item in scan_counts:
        c = item if isinstance(item, dict) else load_counts(item)
        if c["mode"] != "scan":
            raise ValueError(f"{c['sample']}: sinks must be learned from scan-mode counts")
        width = c["placement_bin"]
        positional = {x["name"] for x in c["classes"] if x["kind"] == "positional"}
        total: dict[str, int] = defaultdict(int)
        for p in c["placements"]:
            total[p["class"]] += p["reads"]
        bins: dict[str, dict[tuple[str, int], int]] = defaultdict(dict)
        for p in c["placements"]:
            if p["class"] in positional and p["contig"] != "*":
                bins[p["class"]][(p["contig"], p["start"])] = p["reads"]
                if p["reads"] >= max(min_reads, min_frac * total[p["class"]]):
                    keep[p["class"]].add((p["contig"], p["start"]))
        summaries.append((c["sample"], {k: total[k] for k in positional}, bins))
    rows = []
    for cls, kept in keep.items():
        by = defaultdict(list)
        for contig, start in kept:
            by[contig].append(start)
        for contig, starts in by.items():
            starts.sort()
            s0 = e0 = None
            for st in starts:
                if s0 is None:
                    s0, e0 = st, st + width
                elif st <= e0 + pad:
                    e0 = st + width
                else:
                    rows.append((contig, max(0, s0 - pad), e0 + pad, cls))
                    s0, e0 = st, st + width
            rows.append((contig, max(0, s0 - pad), e0 + pad, cls))
    rows.sort()
    stats: dict[str, dict[str, float]] = defaultdict(dict)
    for sample, total, bins in summaries:
        for cls, tot in total.items():
            iv = [(r[0], r[1], r[2]) for r in rows if r[3] == cls]
            inside = sum(n for (contig, st), n in bins[cls].items() if any(contig == c and a <= st < b for c, a, b in iv))
            stats[cls][sample] = inside / max(tot, 1)
    return rows, dict(stats)


def capture(scan: dict, bed: list[tuple[str, int, int, str]]) -> dict[str, tuple[int, int]]:
    """Per class: (class reads in the scan, of which placed inside that class's sink intervals).
    Unmapped reads count against capture unless fetch mode is run with --unmapped."""
    by: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for contig, s0, e0, cls in bed:
        by[cls].append((contig, s0, e0))
    out: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for p in scan["placements"]:
        out[p["class"]][0] += p["reads"]
        if any(p["contig"] == c and s0 <= p["start"] < e0 for c, s0, e0 in by.get(p["class"], ())):
            out[p["class"]][1] += p["reads"]
    return {k: (v[0], v[1]) for k, v in out.items()}
