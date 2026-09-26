"""Learn class sinks from scan-mode runs.

A multi-copy class has no fixed address in a linear reference: the aligner scatters its reads
over whatever paralogs, decoys and look-alike loci the reference happens to contain. Scan mode
classifies every read by k-mers regardless of placement and records where class reads *were*
placed; the sinks are the (few) intervals that hold essentially all of them. Fetch mode then
retrieves only those intervals. Sinks depend on the reference build and on the aligner, so they
are re-learned, not assumed, whenever either changes.
"""
from __future__ import annotations

import re
from bisect import bisect_right
from collections import defaultdict

DECISION_WIDTH = 10_000
CUT_SHARE = 0.05                                           # control 5' ends / primary reads: about 0.3% in a whole 30x genome


def cut_share(counts: dict) -> float:
    """Control-region 5' ends per mapped primary read: 0.33-0.34% in 50 of the 1000 Genomes 30x
    scans, 26% in the test fixture (a subsample restricted to the controls and sinks)."""
    if not counts.get("primary") or not counts.get("ctrl_reads"):
        return 0.0
    return counts["ctrl_reads"] / counts["primary"]


def looks_cut(counts: dict) -> bool:
    """A scan-mode counts file made from reads cut out along a fetch plan (`ngs-dose plan`) holds
    little but control and sink reads. It says `mode: scan`, and to a sink learner or evaluator it
    would look like a whole file in which every class read lands inside the plan."""
    return cut_share(counts) > CUT_SHARE


def check_scan(counts: dict, allow_cut: bool = False):
    """Raise ValueError unless the counts are a scan of a whole file: what learning sinks and
    measuring their capture both need. A fetch holds only what its plan retrieved, so every class
    read in it lies inside the plan; a cut (`looks_cut`) likewise lacks the reads sinks are for."""
    s = counts.get("sample")
    if counts.get("mode") != "scan":
        raise ValueError(f"{s}: {counts.get('mode')}-mode counts; sinks are learned and evaluated from scan-mode counts only "
                         "(a fetch holds only the reads its plan retrieved)")
    if looks_cut(counts) and not allow_cut:
        raise ValueError(f"{s}: {100 * cut_share(counts):.0f}% of its primary reads lie in the control regions, where a whole-genome file "
                         "has well under 1%: this is a cut along a fetch plan, not a whole-file scan, and it can teach nothing about where reads land")


def read_bed(path) -> list[tuple[str, int, int, str]]:
    """(contig, start, end, class) of a sinks BED, read as the engine reads it (fasta::read_bed):
    plain or gzipped; empty lines and lines starting with '#' or 'track' are skipped, and every
    other line (a 'browser' line or one of spaces too) must have at least three tab-separated
    columns, integer start and end, and end > start. The class is the fourth column as written,
    '' when there is none: such an interval serves every class."""
    from .io import _open
    rows = []
    with _open(path) as fh:
        for i, line in enumerate(fh, 1):
            line = line[:-1] if line.endswith("\n") else line
            if not line or line.startswith(("#", "track")):
                continue
            p = line.split("\t")
            if len(p) < 3:
                raise ValueError(f"{path}:{i}: expected at least 3 tab-separated BED columns")
            if not (_INT.fullmatch(p[1]) and _INT.fullmatch(p[2])):
                raise ValueError(f"{path}:{i}: start and end must be integers")
            s0, e0 = int(p[1]), int(p[2])
            if e0 <= s0:
                raise ValueError(f"{path}:{i}: end <= start")
            rows.append((p[0], s0, e0, p[3] if len(p) > 3 else ""))
    return rows


_INT = re.compile(r"[+-]?[0-9]+")


def _index(bed, classes=()) -> dict[tuple[str, str], tuple[list[int], list[int]]]:
    """(class, contig) -> the class's intervals there, merged where they overlap or touch, as
    sorted start and end lists. A row without a class is added to every name in `classes`."""
    by = defaultdict(list)
    for contig, s0, e0, cls in bed:
        for c in ((cls,) if cls else classes):
            by[(c, contig)].append((s0, e0))
    out = {}
    for key, ivs in by.items():
        ivs.sort()
        st, en = [ivs[0][0]], [ivs[0][1]]
        for a, b in ivs[1:]:
            if a <= en[-1]:
                en[-1] = max(en[-1], b)
            else:
                st.append(a)
                en.append(b)
        out[key] = (st, en)
    return out


def _inside(index, cls: str, contig: str, a: int, b: int) -> bool:
    iv = index.get((cls, contig))
    if not iv:
        return False
    i = bisect_right(iv[0], a) - 1
    return i >= 0 and b <= iv[1][i]


def _bin_end(start: int, width: int, top: int | None) -> int:
    """End of a placement bin, clipped at the contig's end (the last bin of a contig is shorter)."""
    return min(start + width, top) if top else start + width


def _widths(counts: dict, names) -> dict[str, int]:
    """Placement-bin width per class: scans made with different bin widths can be pooled; a
    compositional class is placed on its own, coarser grid."""
    kinds = {x["name"]: x.get("kind") for x in counts.get("classes", [])}
    return {n: (counts["placement_bin"] if kinds.get(n) == "positional" else counts.get("placement_bin_compositional", counts["placement_bin"]))
            for n in names}


def _placed(counts: dict, chosen) -> tuple[dict[str, int], dict[str, dict[tuple[str, int], int]], dict[tuple[str, str, int], int]]:
    """Per class: the reads of the scan (unmapped included), the reads of every placement bin of
    the chosen classes, and their sums per DECISION_WIDTH neighbourhood."""
    total: dict[str, int] = defaultdict(int)
    bins: dict[str, dict[tuple[str, int], int]] = defaultdict(dict)
    coarse: dict[tuple[str, str, int], int] = defaultdict(int)
    for p in counts["placements"]:
        total[p["class"]] += p["reads"]
        if p["class"] in chosen and p["contig"] != "*":
            bins[p["class"]][(p["contig"], p["start"])] = p["reads"]
            coarse[(p["class"], p["contig"], p["start"] // DECISION_WIDTH)] += p["reads"]
    return total, bins, coarse


def _panel_set(counts: dict) -> frozenset:
    h = counts.get("panel_sha256")
    return frozenset([h] if isinstance(h, str) else (h or ()))


def learn(scan_counts, min_frac: float = 1e-5, pad: int = 1000, min_reads: int = 25, classes=None, allow_cut: bool = False,
          allow_mixed_panels: bool = False, log=None) -> tuple[list[tuple[str, int, int, str]], dict]:
    """Intervals holding >= min_frac of a class's reads (and >= min_reads reads, so that a
    small class does not collect every stray placement) in any sample, merged, padded and
    clipped to the contig.

    `scan_counts` is a list of counts dicts or of paths. Files are read one at a time, twice (once
    to learn, once to measure capture), so only the learned intervals are kept between them and a
    whole cohort of scans can be passed as paths; dicts, and any one-shot iterator, are held in
    memory for the second pass. Positional classes always get sinks. A compositional
    class gets them when named in `classes` (a name no scan has is an error). In the NYGC bwa-mem
    1000 Genomes scans the aligner concentrates each satellite family on a small, stable set of
    intervals, some of them on decoy and unplaced contigs: sinks learned from 30 scans held
    >= 0.998 of HSat1A, HSat2, HSat3, aSatHOR, bSat, ACRO, SST1, CER and SATR in each of 200 other
    genomes (>= 0.9985 in two of three random draws) (aSatHOR in 59.6 Mb of intervals, the others in 0.2-3.5 Mb each), and >= 0.966 of
    HSat1B, part of which is left fully unmapped in the 698-genome batch (so HSat1B needs a
    calibrated capture, or `count --unmapped`). The telomeric repeat is concentrated at the
    chromosome ends (in 372 NYGC scans, 92% within 25 kb of an end and 60% in one 10-kb bin of
    chr5p); its sinks held >= 0.9956 in 0.69 Mb. So these classes can be fetched too, but no
    satellite sinks ship with the bundle. Sinks are specific to the aligner and the reference: with
    DRAGEN 4.x and an alt-masked reference, most of the 45S rDNA and DJ reads of the one genome
    checked were left unmapped, and 64-90% of its HSat1A, HSat1B, bSat, ACRO and TEL reads (almost
    none of its HSat2 or aSatHOR). A pipeline must therefore learn its own sinks from a scanned subset
    and check their capture on held-out scans (`ngsdose sinks --evaluate`).

    Returns (BED rows, per-class capture: the fraction of each sample's scan-mode class reads
    that fall inside the learned sinks, scored as in `capture`). A file that holds only a region
    subset (`looks_cut`) is refused unless `allow_cut`: the reads it lacks are the ones sinks are
    for, so it can only confirm what the subset holds. Scans aligned to different references
    (a contig of two lengths) are refused; scans counted with panel files of which neither set
    holds the other are refused unless `allow_mixed_panels`. `log` receives notes (a named class
    that only some scans have)."""
    from .io import load_counts
    items = scan_counts if isinstance(scan_counts, (list, tuple)) else list(scan_counts)
    keep: dict[str, set[tuple[str, int, int]]] = defaultdict(set)      # class -> {(contig, start, end)}
    lengths: dict[str, tuple[int, str]] = {}               # contig -> (length, the first scan that gave it)
    panels: dict[frozenset, str] = {}                      # panel file hashes -> the first scan counted with them
    seen: dict[str, int] = defaultdict(int)
    for item in items:
        c = item if isinstance(item, dict) else load_counts(item)
        check_scan(c, allow_cut)
        for x in c.get("contigs", ()):
            if isinstance(x, dict) and "len" in x:
                n0, s0 = lengths.setdefault(x["name"], (x["len"], c["sample"]))
                if n0 != x["len"]:
                    raise ValueError(f"{c['sample']}: contig {x['name']} is {x['len']:,} bp long, in {s0} {n0:,} bp: the scans were aligned to "
                                     "different references. Sinks are specific to the reference and the aligner: learn one sinks BED per pipeline")
        ps = _panel_set(c)
        if ps and ps not in panels:
            clash = next((s for q, s in panels.items() if not (q <= ps or ps <= q)), None)
            if clash and not allow_mixed_panels:
                raise ValueError(f"{c['sample']} and {clash} were counted with different panel files (neither set of panel_sha256 holds the "
                                 "other): a class counted with other k-mers is found in other places. Learn from scans of one panel set, or "
                                 "pass --allow-mixed-panels")
            panels[ps] = c["sample"]
        kinds = {x["name"]: x["kind"] for x in c["classes"]}
        for n in kinds:
            seen[n] += 1
        chosen = {n for n, k in kinds.items() if k == "positional"} | ({n for n in (classes or ()) if n in kinds})
        width = _widths(c, chosen)
        total, bins, coarse = _placed(c, chosen)
        # whether a neighbourhood is a sink is decided per DECISION_WIDTH, whatever the grid of the
        # scan (a finer grid must not promote three stray reads in one kilobase); a finer grid then
        # trims the sink to the bins that actually hold reads
        for cls, per in bins.items():
            for (contig, start), reads in per.items():
                if coarse[(cls, contig, start // DECISION_WIDTH)] >= max(min_reads, min_frac * total[cls]):
                    keep[cls].add((contig, start, start + width[cls]))
    unknown = [n for n in (classes or ()) if n not in seen]
    if unknown:
        raise ValueError(f"no scan has a class {', '.join(unknown)} (they have {', '.join(sorted(seen))})")
    for n in classes or ():
        if seen[n] < len(items) and log:
            log(f"[sinks] {n} is in {seen[n]} of {len(items)} scans: its sinks are learned from those")
    rows = []
    for cls, kept in keep.items():
        by = defaultdict(list)
        for contig, start, end in kept:
            by[contig].append((start, end))
        for contig, ivs in by.items():
            top = lengths.get(contig, (None,))[0]          # a padded interval does not run past the contig
            ivs.sort()
            if top and ivs[-1][0] >= top:
                raise ValueError(f"{cls} reads placed at {contig}:{ivs[-1][0]:,}, beyond the contig's end ({top:,} bp): the placements "
                                 "and the contig lengths do not belong together")
            s0, e0 = ivs[0]
            for st, en in ivs[1:]:
                if st <= e0 + pad:
                    e0 = max(e0, en)
                else:
                    rows.append((contig, max(0, s0 - pad), min(e0 + pad, top) if top else e0 + pad, cls))
                    s0, e0 = st, en
            rows.append((contig, max(0, s0 - pad), min(e0 + pad, top) if top else e0 + pad, cls))
    rows.sort()
    index = _index(rows)
    stats: dict[str, dict[str, float]] = defaultdict(dict)
    for item in items:                                     # second pass: capture, one file at a time
        c = item if isinstance(item, dict) else load_counts(item)
        kinds = {x["name"]: x["kind"] for x in c["classes"]}
        chosen = sorted({n for n, k in kinds.items() if k == "positional"} | ({n for n in (classes or ()) if n in kinds}))
        width = _widths(c, chosen)
        total, bins, _ = _placed(c, chosen)
        for cls in chosen:
            inside = sum(n for (contig, st), n in bins[cls].items()
                         if _inside(index, cls, contig, st, _bin_end(st, width[cls], lengths.get(contig, (None,))[0])))
            stats[cls][c["sample"]] = inside / max(total[cls], 1)
    return rows, dict(stats)


def capture(scan: dict, bed: list[tuple[str, int, int, str]]) -> dict[str, tuple[int, int]]:
    """Per class: (class reads in the scan, of which placed inside that class's sink intervals).
    A placement bin counts as inside only when all of it, [start, min(start + bin width, contig
    length)], lies inside the class's intervals (merged; a BED row without a class serves every
    class): a bin that merely begins inside a sink also holds reads beyond it. Unmapped reads
    (contig '*') always count against capture here; a fetch reads them only with `count --unmapped`."""
    top = {x["name"]: x["len"] for x in scan.get("contigs", ()) if isinstance(x, dict) and "len" in x}
    names = {p["class"] for p in scan["placements"]}
    width = _widths(scan, names)
    index = _index(bed, names)
    out: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for p in scan["placements"]:
        o = out[p["class"]]
        o[0] += p["reads"]
        if p["contig"] != "*" and _inside(index, p["class"], p["contig"], p["start"], _bin_end(p["start"], width[p["class"]], top.get(p["contig"]))):
            o[1] += p["reads"]
    return {k: (v[0], v[1]) for k, v in out.items()}
