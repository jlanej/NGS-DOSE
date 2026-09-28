"""What a counts file must satisfy before its numbers are used, checked in one place.

A counts file records what it was made with (panel, controls, sinks, what the fetch could not
read); the estimator must not turn a class into a copy number or a mass when that record says
the class was not fully counted, or was counted with a different panel than the bundle's. Keys
that later engine versions add (`sinks_skipped`, `pipeline`, `pad`) are optional: files written
before them carry none and pass.
"""
from __future__ import annotations

import gzip
import hashlib
from bisect import bisect_right
from collections import defaultdict

from . import io, sinks
from .io import CountsError  # noqa: F401  (raised by load)


def load(path) -> dict:
    """One counts file; any failure to read or decode it is a CountsError naming the path."""
    return io.load_counts(path)


def missing_from_bed(counts: dict, known_sinks: dict | None = None) -> list[str]:
    """Classes a fetch counted that its sinks BED gives no interval, when that BED is known here by its
    SHA-256 (`known_sinks`, sha256 -> rows, as estimate.sinks_by_hash makes it). Engines since 645ae55
    refuse such a fetch or list the classes in `sinks_missing_classes`; fae1124 counts them without a
    word, and only the BED tells. A BED that names no class serves every class."""
    if counts.get("mode") != "fetch" or not known_sinks:
        return []
    bed = known_sinks.get(counts.get("sinks_sha256"))
    named = {r[3] for r in bed or () if r[3]}
    if not named:
        return []
    return [c["name"] for c in counts.get("classes", []) if c["name"] not in named]


def incomplete_sinks(counts: dict, known_sinks: dict | None = None) -> dict[str, str]:
    """Classes a fetch did not read in full: no sink interval at all (`--allow-missing-sinks`, or
    missing_from_bed when the fetch's sinks BED is known), or sink intervals on contigs the input
    header lacks, which the engine drops (`sinks_skipped`). The engine keys dropped rows without a
    class as '': such a row serves every class, so every class of the file is marked."""
    out: dict[str, str] = {}
    skipped = dict(counts.get("sinks_skipped") or {})
    anyclass = skipped.pop("", None)
    counted = {c["name"] for c in counts.get("classes", [])}
    for name, v in skipped.items():
        if counted and name not in counted:                # a class of the sinks BED this fetch did not load
            continue
        out[name] = (f"{v.get('intervals', '?')} sink intervals ({v.get('bp', '?')} bp) are on contigs absent from the "
                     "input header and were not fetched")
    if anyclass is not None:
        why = (f"{anyclass.get('intervals', '?')} sink intervals without a class ({anyclass.get('bp', '?')} bp), which "
               "serve every class, are on contigs absent from the input header and were not fetched")
        for name in dict.fromkeys(c["name"] for c in counts.get("classes", [])):
            out[name] = f"{out[name]}; {why}" if name in out else why
    for name in missing_from_bed(counts, known_sinks):
        out[name] = (f"no sink intervals in the fetch: its sinks BED ({counts.get('sinks')}) gives {name} none, and the engine did "
                     "not record it (fae1124 does not)")
    for name in counts.get("sinks_missing_classes") or []:
        out[name] = "no sink intervals in the fetch"
    return out


def _overlaps(merged, a: int, b: int) -> bool:
    """Whether [a, b) overlaps any of `merged` (sorted, merged start and end lists, as sinks._index keeps them)."""
    if not merged:
        return False
    j = bisect_right(merged[1], a)                         # the first interval that ends after a
    return j < len(merged[1]) and merged[0][j] < b


def subset_only(counts: dict, sub_options: dict, known_sinks: dict | None = None) -> dict[str, str]:
    """Classes a fetch read only inside the intervals of their sub-options (estimate.SubOption: named
    subsets of a class's sinks, such as DXZ1 of aSatHOR), not their whole sinks: the class's count is
    then the sub-options' and not the class's. When the fetch's sinks BED is known (`known_sinks`,
    sha256 -> rows) it says so; otherwise the placements do: every placement bin of the class overlaps
    a sub-option's interval, where a fetch of the class's whole sinks has placements across its many
    intervals. Scans are whole and never subset-only."""
    if counts.get("mode") != "fetch" or not sub_options or "placement_bin" not in counts:
        return {}
    parts = defaultdict(list)
    for o in sub_options.values():
        parts[o.parent].append(o)
    counted = {c["name"] for c in counts.get("classes", [])}
    bed = (known_sinks or {}).get(counts.get("sinks_sha256"))
    out = {}
    for parent, subs in parts.items():
        if parent not in counted:
            continue
        index = sinks._index([(c, s, e, parent) for o in subs for c, s, e in o.intervals])
        width = sinks._widths(counts, [parent])[parent]
        names = ", ".join(sorted(o.name for o in subs))
        if bed is not None:
            top = {x["name"]: x["len"] for x in counts.get("contigs", ()) if isinstance(x, dict) and "len" in x}
            # a row without a class serves every class (the engine, incomplete_sinks and estimate.sub_option_values
            # read it so): it is one of the parent's intervals here too
            own = [(c, s, min(e, top.get(c, e))) for c, s, e, n in bed if n in (parent, "") and s < top.get(c, e)]
            if own and all(sinks._inside(index, parent, c, s, e) for c, s, e in own):
                out[parent] = (f"the fetch's sinks BED holds no {parent} interval outside those of its sub-option(s) {names}: {parent} was "
                               "read only there (the sub-options' values are reported, the class's is not)")
            continue
        placed = [p for p in counts.get("placements", []) if p["class"] == parent and p["contig"] != "*"]
        # no placed bin says nothing (all([]) would): the class is then unverified, not subset-only
        if placed and all(_overlaps(index.get((parent, p["contig"])), p["start"], p["start"] + width) for p in placed):
            out[parent] = (f"every {parent} placement of this fetch ({len(placed)} bins) lies at the intervals of its sub-option(s) {names}: "
                           f"the fetch read those, not {parent}'s sinks (the sub-options' values are reported, the class's is not)")
    return out


def unverified_parents(counts: dict, sub_options: dict, known_sinks: dict | None = None, only: dict | None = None) -> dict[str, str]:
    """Classes with sub-options (estimate.SubOption) that a fetch counted with a sinks BED not known here
    (`known_sinks`, sha256 -> rows): whether it read their whole sinks or only the sub-options'
    intervals cannot be told, and a count at the X and Y arrays alone would read as the whole family's.
    The placement check of subset_only does not settle it: reads of the class land in other fetched
    intervals (controls, other classes' sinks) too. Classes subset_only already names (`only`) are left
    to it. Scans, and fetches whose BED is known, give none."""
    if counts.get("mode") != "fetch" or not sub_options:
        return {}
    if (known_sinks or {}).get(counts.get("sinks_sha256")) is not None:
        return {}
    parents = {o.parent: [] for o in sub_options.values()}
    for o in sub_options.values():
        parents[o.parent].append(o.name)
    counted = {c["name"] for c in counts.get("classes", [])}
    return {p: (f"a fetch made with a sinks BED not known here ({counts.get('sinks')}, sha256 {str(counts.get('sinks_sha256'))[:12]}...): "
                f"it may have read {p}'s whole sinks or only the intervals of its sub-option(s) {', '.join(sorted(subs))}, so its value is "
                "not reported; pass that BED to `ngsdose estimate --fetch-sinks`")
            for p, subs in parents.items() if p in counted and p not in (only or {})}


def pipeline_check(counts: dict, learned: dict | None) -> str | None:
    """Whether a fetch's input came from the pipeline its sinks were learned under. `learned` is the
    bundle's `sinks_learned_from.pipeline`: the aligner's @PG program name and the engine's `sq_sha256`
    of the reference's @SQ lines (names, lengths, M5). Sinks are where one aligner put a class's reads
    on one reference; another aligner, or another reference, puts them elsewhere, and a fetch through
    these sinks then undercounts with no other sign (DESIGN.md section 12: under DRAGEN 4.x most 45S
    reads are unmapped). None when either side lacks the record (counts written by engines before it,
    a bundle without one, a scan); 'ok' when they agree; otherwise what differs."""
    if not learned or counts.get("mode") != "fetch":
        return None
    rec = counts.get("pipeline") or {}
    if not rec:
        return None
    bad, checked = [], False
    want = str(learned.get("aligner") or "").lower()
    pns = sorted({str(p.get("pn") or p.get("id") or "").lower() for p in rec.get("pg", [])} - {""})
    if want and pns:                                       # a header without @PG lines (a stripped test file) names no aligner
        checked = True
        if not any(want in pn for pn in pns):
            bad.append(f"aligned by {', '.join(pns)}, not {learned['aligner']}, which the sinks were learned under")
    # the @SQ hash is comparable only between whole headers: a file whose header keeps a subset of the contigs
    # (a cut, the test fixture) hashes differently without being from another reference
    if learned.get("sq_sha256") and rec.get("sq_sha256") and (learned.get("sq_n") is None or rec.get("sq_n") == learned.get("sq_n")):
        checked = True
        if rec["sq_sha256"] != learned["sq_sha256"]:
            bad.append("its reference's @SQ lines (names, lengths, M5) are not those the sinks were learned under")
    if not checked:
        return None
    return "ok" if not bad else "; ".join(bad)


def _panel(x) -> io.Panel:
    if isinstance(x, io.Panel):
        return x
    if getattr(x, "_panel_obj", None) is None:
        x._panel_obj = io.load_panel(x.panel)
    return x._panel_obj


def panel_mismatch(counts: dict, bundle) -> dict[str, str]:
    """Positional classes counted with k-mers other than those of the bundle panel (`bundle`: a
    Bundle or a loaded Panel). The estimator's callable mask comes from the bundle panel, so a
    class counted with fewer k-mers (an older panel, or k-mers a co-loaded panel shared and the
    engine dropped) would read low with no other sign. Fields a counts file lacks are not checked."""
    panel = _panel(bundle)
    out = {}
    k = counts.get("k")
    for c in counts.get("classes", []):
        pc = panel.classes.get(c["name"])
        if c.get("kind") != "positional" or pc is None:
            continue
        bad = []
        if pc.kind != "positional":
            bad.append(f"the bundle panel has it as {pc.kind}")
        if k is not None and k != panel.k:
            bad.append(f"k={k}, the bundle panel k={panel.k}")
        if c.get("length") is not None and c["length"] != pc.length:
            bad.append(f"unit length {c['length']}, the bundle panel {pc.length}")
        if c.get("panel_kmers") is not None and c["panel_kmers"] != len(pc.kmer_pos):
            bad.append(f"{c['panel_kmers']} panel k-mers, the bundle panel {len(pc.kmer_pos)}")
        if bad:
            out[c["name"]] = "counted with a different panel (" + "; ".join(bad) + ")"
    return out


def panel_hashes(path) -> set[str]:
    """sha256 of a panel file as stored and of its decompressed content: the engine records the
    first; the second still matches after the file is recompressed."""
    raw = open(path, "rb").read()
    out = {hashlib.sha256(raw).hexdigest()}
    if raw[:2] == b"\x1f\x8b":
        out.add(hashlib.sha256(gzip.decompress(raw)).hexdigest())
    return out


def issues(counts: dict, bundle=None, known_sinks: dict | None = None) -> list[tuple[str, str]]:
    """(level, message) for everything that makes a counts file unusable ('fatal') or parts of it
    unmeasured or doubtful ('warn'). With a bundle, also what ties the file to that bundle.
    `known_sinks`: sinks BEDs by SHA-256 (incomplete_sinks)."""
    s = counts.get("sample", "?")
    out: list[tuple[str, str]] = []
    if counts.get("format") != io.COUNTS_FORMAT:
        return [("fatal", f"{s}: not an {io.COUNTS_FORMAT} file (format={counts.get('format')!r})")]
    names = [c["name"] for c in counts.get("classes", [])]
    for n in sorted({n for n in names if names.count(n) > 1}):
        out.append(("fatal", f"{s}: class {n} appears more than once"))
    eof = counts.get("eof_marker")
    if eof == "absent":
        out.append(("warn", f"{s}: counted from a file without an end-of-file marker (truncated input)"))
    elif eof == "unchecked":
        out.append(("warn", f"{s}: the input's end-of-file marker could not be checked (a stream or a server without "
                            "range requests): a truncated input would not have been noticed"))
    for name, why in incomplete_sinks(counts, known_sinks).items():
        out.append(("warn", f"{s}: {name} is not estimated: {why}"))
    if bundle is None:
        return out
    from .estimate import check_build
    try:
        check_build(counts, bundle.contig_lengths())
    except ValueError as e:
        out.append(("fatal", str(e)))
    regions = bundle.regions()
    ctrl = {n for n, role in regions if role == "control"}
    have = counts.get("regions", [])
    have_ctrl = {r["name"] for r in have if r.get("role", "control") == "control"}
    # a lighter controls file the bundle names (Bundle.control_subsets) is accepted; no other set is
    subsets = bundle.control_subsets() if hasattr(bundle, "control_subsets") else {}
    if have_ctrl != ctrl and not (have_ctrl < ctrl and any(set(v) == have_ctrl for v in subsets.values())):
        why = ("a subset of them that no named subset of the bundle declares" if have_ctrl and have_ctrl < ctrl
               else "the control regions differ from the bundle's")
        out.append(("fatal", f"{s}: made with a different controls file ({why})"))
    known = {n for n, _ in regions}
    unknown = [r["name"] for r in have if r["name"] not in known and r.get("role", "control") != "control"]
    if unknown:
        out.append(("warn", f"{s}: {len(unknown)} non-control regions unknown to the bundle are left out ({', '.join(unknown[:3])}"
                            f"{', ...' if len(unknown) > 3 else ''})"))
    for name, why in panel_mismatch(counts, bundle).items():
        out.append(("warn", f"{s}: {name} is not estimated: {why}"))
    panel, units = _panel(bundle), set(bundle.meta.get("units", {}))
    ex = bundle.experimental() if hasattr(bundle, "experimental") else None
    for c in counts.get("classes", []):
        if c.get("kind") == "positional" and (c["name"] not in panel.classes or c["name"] not in units):
            # an experimental unit estimates it (status `experimental`); otherwise say what was missing
            found = ex.lookup(c["name"], counts, panel.classes.get(c["name"])) if ex is not None else None
            if isinstance(found, tuple):
                continue
            if found is None and ex is not None:
                found = ex.missing(c["name"])
            out.append(("warn", f"{s}: positional class {c['name']} is not estimated: the bundle has no "
                                f"{'panel entry' if c['name'] not in panel.classes else 'unit sequence'} for it"
                                + (f"; {found}" if found else "")))
    if getattr(bundle, "_panel_hashes", None) is None:
        bundle._panel_hashes = panel_hashes(bundle.panel)
    recorded = counts.get("panel_sha256")
    recorded = [recorded] if isinstance(recorded, str) else (recorded or [])
    if recorded and not bundle._panel_hashes & set(recorded):
        out.append(("warn", f"{s}: none of the panels it was counted with is the bundle's panel file ({bundle.panel.name})"))
    learned = (bundle.meta.get("sinks_learned_from") or {}).get("pipeline") if hasattr(bundle, "meta") else None
    why = pipeline_check(counts, learned)
    if why and why != "ok":
        out.append(("warn", f"{s}: this fetch's input is not from the pipeline the sinks were learned under ({why}). The sinks may not "
                            "hold its class reads; learn sinks from whole-file scans of this pipeline (ngsdose sinks) before trusting its "
                            "fetches (the table records the verdict in sinks_pipeline)"))
    return out
