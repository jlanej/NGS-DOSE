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

from . import io
from .io import CountsError  # noqa: F401  (raised by load)


def load(path) -> dict:
    """One counts file; any failure to read or decode it is a CountsError naming the path."""
    return io.load_counts(path)


def incomplete_sinks(counts: dict) -> dict[str, str]:
    """Classes a fetch did not read in full: no sink interval at all (`--allow-missing-sinks`), or
    sink intervals on contigs the input header lacks, which the engine drops (`sinks_skipped`).
    The engine keys dropped rows without a class as '': such a row serves every class, so every
    class of the file is marked."""
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
    for name in counts.get("sinks_missing_classes") or []:
        out[name] = "no sink intervals in the fetch"
    return out


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


def issues(counts: dict, bundle=None) -> list[tuple[str, str]]:
    """(level, message) for everything that makes a counts file unusable ('fatal') or parts of it
    unmeasured or doubtful ('warn'). With a bundle, also what ties the file to that bundle."""
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
    for name, why in incomplete_sinks(counts).items():
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
    if {r["name"] for r in have if r.get("role", "control") == "control"} != ctrl:
        out.append(("fatal", f"{s}: made with a different controls file (the control regions differ from the bundle's)"))
    known = {n for n, _ in regions}
    unknown = [r["name"] for r in have if r["name"] not in known and r.get("role", "control") != "control"]
    if unknown:
        out.append(("warn", f"{s}: {len(unknown)} non-control regions unknown to the bundle are left out ({', '.join(unknown[:3])}"
                            f"{', ...' if len(unknown) > 3 else ''})"))
    for name, why in panel_mismatch(counts, bundle).items():
        out.append(("warn", f"{s}: {name} is not estimated: {why}"))
    panel, units = _panel(bundle), set(bundle.meta.get("units", {}))
    for c in counts.get("classes", []):
        if c.get("kind") == "positional" and (c["name"] not in panel.classes or c["name"] not in units):
            out.append(("warn", f"{s}: positional class {c['name']} is not estimated: the bundle has no "
                                f"{'panel entry' if c['name'] not in panel.classes else 'unit sequence'} for it"))
    if getattr(bundle, "_panel_hashes", None) is None:
        bundle._panel_hashes = panel_hashes(bundle.panel)
    recorded = counts.get("panel_sha256")
    recorded = [recorded] if isinstance(recorded, str) else (recorded or [])
    if recorded and not bundle._panel_hashes & set(recorded):
        out.append(("warn", f"{s}: none of the panels it was counted with is the bundle's panel file ({bundle.panel.name})"))
    return out
