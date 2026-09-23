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

DECISION_WIDTH = 10_000


def learn(scan_counts, min_frac: float = 1e-5, pad: int = 1000, min_reads: int = 25, classes=None) -> tuple[list[tuple[str, int, int, str]], dict]:
    """Intervals holding >= min_frac of a class's reads (and >= min_reads reads, so that a
    small class does not collect every stray placement) in any sample, merged, padded and
    clipped to the contig.

    `scan_counts` is an iterable of counts dicts or of paths; files are read one at a time, so a
    whole cohort of scans can be passed. Positional classes always get sinks. A compositional
    class gets them only when named in `classes`: a satellite family's reads are spread over too
    much of the alignment for targeted retrieval to make sense, but the telomeric repeat is not
    dispersed - the aligner concentrates its reads at the chromosome ends (in 372 NYGC scans,
    92% within 25 kb of an end and 60% in one 10-kb bin of chr5p), so it is fetched like a
    positional class. Returns (BED rows, per-class capture: the fraction of each sample's
    scan-mode class reads that fall inside the learned sinks)."""
    from .io import load_counts
    keep: dict[str, set[tuple[str, int, int]]] = defaultdict(set)      # class -> {(contig, start, end)}
    summaries = []                                         # per sample: class totals and the reads of every bin
    lengths: dict[str, int] = {}
    for item in scan_counts:
        c = item if isinstance(item, dict) else load_counts(item)
        if c["mode"] != "scan":
            raise ValueError(f"{c['sample']}: sinks must be learned from scan-mode counts")
        kinds = {x["name"]: x["kind"] for x in c["classes"]}
        chosen = {n for n, k in kinds.items() if k == "positional"} | ({n for n in (classes or ()) if n in kinds})
        # scans made with different bin widths can be pooled; a compositional class is placed on its own, coarser grid
        width = {n: (c["placement_bin"] if kinds[n] == "positional" else c.get("placement_bin_compositional", c["placement_bin"])) for n in chosen}
        lengths.update({x["name"]: x["len"] for x in c.get("contigs", ()) if isinstance(x, dict) and "len" in x})
        total: dict[str, int] = defaultdict(int)
        for p in c["placements"]:
            total[p["class"]] += p["reads"]
        bins: dict[str, dict[tuple[str, int], int]] = defaultdict(dict)
        coarse: dict[tuple[str, str, int], int] = defaultdict(int)
        for p in c["placements"]:
            if p["class"] in chosen and p["contig"] != "*":
                bins[p["class"]][(p["contig"], p["start"])] = p["reads"]
                coarse[(p["class"], p["contig"], p["start"] // DECISION_WIDTH)] += p["reads"]
        # whether a neighbourhood is a sink is decided per DECISION_WIDTH, whatever the grid of the
        # scan (a finer grid must not promote three stray reads in one kilobase); a finer grid then
        # trims the sink to the bins that actually hold reads
        for cls, per in bins.items():
            for (contig, start), reads in per.items():
                if coarse[(cls, contig, start // DECISION_WIDTH)] >= max(min_reads, min_frac * total[cls]):
                    keep[cls].add((contig, start, start + width[cls]))
        summaries.append((c["sample"], {k: total[k] for k in chosen}, bins))
    rows = []
    for cls, kept in keep.items():
        by = defaultdict(list)
        for contig, start, end in kept:
            by[contig].append((start, end))
        for contig, ivs in by.items():
            top = lengths.get(contig)                      # a padded interval does not run past the contig
            ivs.sort()
            s0, e0 = ivs[0]
            for st, en in ivs[1:]:
                if st <= e0 + pad:
                    e0 = max(e0, en)
                else:
                    rows.append((contig, max(0, s0 - pad), min(e0 + pad, top) if top else e0 + pad, cls))
                    s0, e0 = st, en
            rows.append((contig, max(0, s0 - pad), min(e0 + pad, top) if top else e0 + pad, cls))
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
