"""What a fetch reads: the bytes of a CRAM that a set of intervals makes a reader decode, from its index.

A CRAM is read slice by slice: a region query decodes every slice whose alignment span overlaps
the region, whole, and a slice that overlaps several regions is read once. The .crai index lists
every slice (gzip-compressed TSV: reference id, 1-based alignment start, alignment span, container
byte offset, slice byte offset within the container's data, slice size in bytes). The reference id
is the contig's position in the file's @SQ header, so the contig names come from the header of the
CRAMs (or the .fai or .dict of the reference they were aligned to, which lists the same contigs in
the same order). A counts file's `contigs` list will not do: it holds only the contigs that had
reads. A slice that holds reads of several contigs (CRAM's multi-reference slices, used for the
small decoy and unplaced contigs where many sinks are) is listed once per contig with the same
offsets; its bytes are counted once. Slices of reads without a coordinate are listed under id -1:
the unmapped bin, which a fetch reads only with `ngs-dose count --unmapped`.

The bytes of a set of slices are their sizes plus, once per container read, the bytes between the
container's data start and its first listed slice (its compression header). Container headers
themselves (a few dozen bytes) are not in the index and not counted. A .crai says nothing of the
file it indexes: take it from the same place as the CRAM (an index of another copy of the same
sample gives other numbers), and check a slice decodes if in doubt. BAM (.bai, .csi) is not
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

    def keys(self, intervals) -> tuple[set[tuple[int, int]], int]:
        """The slices that intervals ((contig, 0-based start, end), contig UNMAPPED for the unmapped
        bin) overlap, and the number of intervals on contigs the header does not have (a fetch
        leaves those out)."""
        out, absent = set(), 0
        for contig, s0, e0 in intervals:
            if contig == UNMAPPED:
                out |= self.unmapped
                continue
            if contig not in self.contigs:
                absent += 1
                continue
            v = self.slices.get(contig)
            if not v:
                continue
            starts, ends, keys, top = v
            i, j = bisect_right(top, s0), bisect_left(starts, e0)       # the first slice that could end after s0; past the last that starts before e0
            out.update(keys[x] for x in range(i, j) if ends[x] > s0)
        return out, absent

    def bytes(self, keys) -> int:
        return sum(self.size[k] for k in keys) + sum(self.head[c] for c in {k[0] for k in keys})


def median(v):
    v = sorted(v)
    n = len(v)
    return (v[n // 2] + v[(n - 1) // 2]) / 2 if n else float("nan")


def component_costs(index: CraiIndex, components: dict[str, list[tuple[str, int, int]]]) -> dict[str, dict]:
    """Per named component and for their union ('union', each slice read once), and cumulatively in
    the order given ('cum_bytes' of each component: it and all before it): bytes, slices and the
    number of intervals on contigs the header lacks."""
    out, seen = {}, set()
    for name, iv in components.items():
        k, absent = index.keys(iv)
        seen |= k
        out[name] = dict(bytes=index.bytes(k), slices=len(k), absent=absent, cum_bytes=index.bytes(seen), cum_slices=len(seen))
    out["union"] = dict(bytes=index.bytes(seen), slices=len(seen), absent=sum(v["absent"] for v in out.values()),
                        cum_bytes=index.bytes(seen), cum_slices=len(seen))
    return out
