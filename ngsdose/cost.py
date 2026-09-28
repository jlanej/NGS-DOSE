"""What a fetch reads: the bytes of a CRAM that a set of intervals makes `ngs-dose count -m fetch` decode, from its index.

A CRAM is read slice by slice: a region query decodes every slice whose alignment span overlaps
the region, whole. The engine merges the intervals of a plan where they touch or overlap
(src/count.rs `fetch_plan`: per contig, sorted, a start at or before the running end joins) and
makes one indexed fetch per merged interval (`read_interval`), and htslib decodes a slice again for
every fetch whose interval overlaps it. So a fetch of a plan reads, per merged interval, every slice
overlapping it plus the compression header of each container those slices are in, and a slice under
k merged intervals that do not touch is decoded k times: that is what `CraiIndex.price` charges (the
figure the plan reports), computed with the engine's own merge rule (`merge`). Reading every slice
once (`price_once`) is the floor a reader that sorted the plan's slices would reach; the shipped
engine does not, and the gap is largest where a few slices hold many small intervals (the sparse
chrY slices of a woman's CRAM, decoded once per chrY truth region). A consequence of the rule: an
interval that bridges two intervals of the plan can make it cheaper (one decode of their shared
slices instead of two), so adding an option does not always make a plan dearer.

The .crai index lists every slice (gzip-compressed TSV: reference id, 1-based alignment start,
alignment span, container byte offset, slice byte offset within the container's data, slice size in
bytes). The reference id is the contig's position in the file's @SQ header, so the contig names come
from the header of the CRAMs (or the .fai or .dict of the reference they were aligned to, which
lists the same contigs in the same order). A counts file's `contigs` list will not do: it holds only
the contigs that had reads. A slice that holds reads of several contigs (CRAM's multi-reference
slices, used for the small decoy and unplaced contigs where many sinks are) is listed once per contig
with the same offsets: it is one slice, counted once within one fetch, and decoded again by a fetch
on each of its other contigs. Slices of reads without a coordinate are listed under id -1: the
unmapped bin, which a fetch reads once, whole, only with `ngs-dose count --unmapped`.

The bytes of one fetch are the sizes of the slices it overlaps plus, once per container among them,
the bytes between the container's data start and its first listed slice (its compression header).
Container headers themselves (a few dozen bytes) are not in the index and not counted. A .crai says
nothing of the file it indexes: take it from the same place as the CRAM (an index of another copy of
the same sample gives other numbers), and check a slice decodes if in doubt. BAM (.bai, .csi) is not
supported: the NYGC 1000 Genomes files and the biobank files this is for are CRAMs.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict

UNMAPPED = "*"                                             # the contig name of the unmapped bin in an interval list


def read_contigs(path) -> list[tuple[str, int | None]]:
    """(name, length) in @SQ order from a .fai, a .dict or SAM header (`samtools view -H`), or a
    list of names one per line (length None)."""
    from .io import _open
    sq, fai, names = [], [], []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("@"):
                if line.startswith("@SQ\t"):
                    tags = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
                    if "SN" not in tags:
                        raise ValueError(f"{path}: an @SQ line without SN")
                    sq.append((tags["SN"], int(tags["LN"]) if tags.get("LN", "").isdigit() else None))
                continue
            p = line.split("\t")
            if len(p) >= 2 and p[1].isdigit():
                fai.append((p[0], int(p[1])))
            else:
                names.append((line.split()[0], None))
    got = sq or fai or names
    if not got or sum(map(bool, (sq, fai, names))) > 1:
        raise ValueError(f"{path}: expected a .fai, a .dict or SAM header, or one contig name per line")
    return got


def merge(intervals) -> dict[str, tuple[list[int], list[int]]]:
    """The fetches `ngs-dose count -m fetch` makes for a set of intervals ((contig, 0-based start,
    end), contig UNMAPPED for the unmapped bin): per contig, sorted and merged where they touch or
    overlap, exactly as src/count.rs `fetch_plan` merges them (a start at or before the running end
    joins), as (starts, ends). The unmapped bin is one fetch, ([0], [0])."""
    by: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for c, s, e in intervals:
        by[c].append((s, e))
    out = {}
    for c, v in by.items():
        v.sort()
        starts: list[int] = []
        ends: list[int] = []
        for s, e in v:
            if ends and s <= ends[-1]:
                if e > ends[-1]:
                    ends[-1] = e
            else:
                starts.append(s)
                ends.append(e)
        out[c] = (starts, ends)
    return out


class CraiIndex:
    """The slices of one CRAM, from its .crai and the contig list of its header (`read_contigs`)."""

    def __init__(self, path, contigs: list[tuple[str, int | None]]):
        from .io import _open
        self.path = str(path)
        self.size: dict[tuple[int, int], int] = {}         # (container offset, slice offset) -> slice bytes
        head: dict[int, int] = {}                          # container offset -> offset of its first listed slice
        per: dict[str, list[tuple[int, int, tuple[int, int]]]] = defaultdict(list)
        self.unmapped: set[tuple[int, int]] = set()
        with _open(path) as fh:
            for i, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                try:
                    ref, st, span, coff, soff, size = map(int, line.split("\t"))
                except ValueError:
                    raise ValueError(f"{path}:{i}: not a .crai line (six integers, tab-separated)") from None
                key = (coff, soff)
                if self.size.setdefault(key, size) != size:
                    raise ValueError(f"{path}:{i}: the slice at {coff}+{soff} is listed with two sizes")
                head[coff] = min(head.get(coff, soff), soff)
                if ref == -1:
                    self.unmapped.add(key)
                    continue
                if not 0 <= ref < len(contigs):
                    raise ValueError(f"{path}:{i}: reference id {ref}, but the contig list has {len(contigs)} contigs: "
                                     "the index and the contig list do not belong together")
                name, length = contigs[ref]
                a, b = st - 1, st - 1 + span
                if length is not None and b > length + 1000:
                    raise ValueError(f"{path}:{i}: a slice of {name} ends at {b:,}, beyond the contig's {length:,} bp: the index "
                                     "and the contig list do not belong together (the @SQ order of the CRAM's header is needed)")
                per[name].append((a, b, key))
        self.head = head
        self.total = sum(self.size.values()) + sum(head.values())
        self.slices: dict[str, tuple[list[int], list[int], list[tuple[int, int]], list[int]]] = {}
        for name, v in per.items():
            v.sort()
            top, run = [], 0
            for _, b, _ in v:
                run = max(run, b)
                top.append(run)                            # running maximum of the ends: bisectable
            self.slices[name] = ([a for a, _, _ in v], [b for _, b, _ in v], [k for _, _, k in v], top)
        self.contigs = {n for n, _ in contigs}
        self._fetch: dict[tuple[str, int, int], int] = {}  # fetch(): the same merged interval is priced in plan after plan

    def overlapping(self, contig, s0, e0) -> set[tuple[int, int]]:
        """The slices one interval overlaps (the unmapped bin's for UNMAPPED; none on a contig the
        header lacks or that has no slices)."""
        if contig == UNMAPPED:
            return set(self.unmapped)
        v = self.slices.get(contig)
        if not v:
            return set()
        starts, ends, keys, top = v
        i, j = bisect_right(top, s0), bisect_left(starts, e0)       # the first slice that could end after s0; past the last that starts before e0
        return {keys[x] for x in range(i, j) if ends[x] > s0}

    def keys(self, intervals) -> tuple[set[tuple[int, int]], int]:
        """Every slice the intervals overlap, each once (the floor's slices), and the number of
        intervals on contigs the header does not have (a fetch leaves those out)."""
        out, absent = set(), 0
        for contig, s0, e0 in intervals:
            if contig != UNMAPPED and contig not in self.contigs:
                absent += 1
                continue
            out |= self.overlapping(contig, s0, e0)
        return out, absent

    def bytes(self, keys) -> int:
        """The bytes of decoding a set of slices together: their sizes and the compression header
        of each container among them, once."""
        return sum(self.size[k] for k in keys) + sum(self.head[c] for c in {k[0] for k in keys})

    def fetch(self, contig, s0, e0) -> int:
        """The bytes one indexed fetch of an interval reads (`bytes` of the slices it overlaps)."""
        k = (contig, s0, e0)
        got = self._fetch.get(k)
        if got is None:
            got = self._fetch[k] = self.bytes(self.overlapping(contig, s0, e0))
        return got

    def price(self, intervals) -> int:
        """What `ngs-dose count -m fetch` reads for a set of intervals: merged as the engine merges
        them (`merge`), one fetch per merged interval, a slice under several of them decoded once per
        fetch."""
        return sum(self.fetch(c, s, e) for c, (starts, ends) in merge(intervals).items() for s, e in zip(starts, ends))

    def price_once(self, intervals) -> int:
        """The floor: every slice the intervals overlap decoded once, with each container's
        compression header once. What a reader that sorted the plan's slices would read; the shipped
        engine reads `price`."""
        return self.bytes(self.keys(intervals)[0])

    def extra(self, plan: dict[str, tuple[list[int], list[int]]], contig, s0, e0) -> int:
        """The bytes the fetches of `plan` (from `merge`) grow by when an interval is added: the
        fetch of the merged interval it joins or makes, less the fetches of the plan's intervals it
        joined. 0 when it lies inside one of them (or overlaps no slice), negative when it bridges
        two that decode the same slices."""
        if contig == UNMAPPED:
            return 0 if UNMAPPED in plan else self.fetch(UNMAPPED, 0, 0)
        starts, ends = plan.get(contig, ([], []))
        i, j = bisect_left(ends, s0), bisect_right(starts, e0)      # the plan's intervals ending at or after s0 and starting at or before e0: those it joins
        if i >= j:
            return self.fetch(contig, s0, e0)
        return self.fetch(contig, min(s0, starts[i]), max(e0, ends[j - 1])) - sum(self.fetch(contig, starts[x], ends[x]) for x in range(i, j))


def median(v):
    v = sorted(v)
    n = len(v)
    return (v[n // 2] + v[(n - 1) // 2]) / 2 if n else float("nan")


def component_costs(index: CraiIndex, components: dict[str, list[tuple[str, int, int]]]) -> dict[str, dict]:
    """Per named component, for all of them together ('union') and cumulatively in the order given
    ('cum_*': it and all before it): what the engine reads ('bytes': `price`), the floor with every
    slice once ('bytes_once': `price_once`), the slices overlapped and the number of intervals on
    contigs the header lacks."""
    out: dict[str, dict] = {}
    seen: list[tuple[str, int, int]] = []
    for name, iv in components.items():
        k, absent = index.keys(iv)
        seen = seen + list(iv)
        out[name] = dict(bytes=index.price(iv), bytes_once=index.bytes(k), slices=len(k), absent=absent,
                         cum_bytes=index.price(seen), cum_bytes_once=index.price_once(seen))
    total = dict(bytes=index.price(seen), bytes_once=index.price_once(seen), slices=len(index.keys(seen)[0]),
                 absent=sum(v["absent"] for v in out.values()))
    out["union"] = dict(total, cum_bytes=total["bytes"], cum_bytes_once=total["bytes_once"])
    return out
