"""Readers for the engine's counts JSON, k-mer panels and FASTA resources."""
from __future__ import annotations

import gzip
import io as _stdio
import json
import zlib
from dataclasses import dataclass

import numpy as np

COUNTS_FORMAT = "ngs-dose-counts/1"


class CountsError(ValueError):
    """A counts file that cannot be read or is not a counts file; the message starts with its path."""


def _open(path):
    """Text handle on a plain or gzip-compressed file, told apart by the gzip magic bytes, not the name.
    The file is opened once and the magic bytes are peeked at, so a pipe or <(...) works too."""
    fh = open(path, "rb")
    try:
        gz = fh.peek(2)[:2] == b"\x1f\x8b"
    except BaseException:
        fh.close()
        raise
    if not gz:
        return _stdio.TextIOWrapper(fh)
    g = gzip.GzipFile(fileobj=fh, mode="rb")
    g.myfileobj = fh                                       # closing the handle closes the file, as gzip.open(path) does
    return _stdio.TextIOWrapper(g)


def load_counts(path) -> dict:
    """Load one per-sample counts file written by `ngs-dose count`."""
    try:
        with _open(path) as fh:
            d = json.load(fh)
    except (OSError, EOFError, zlib.error, UnicodeDecodeError, json.JSONDecodeError) as e:
        kind = "zlib.error, a damaged gzip stream" if isinstance(e, zlib.error) else type(e).__name__
        raise CountsError(f"{path}: unreadable counts file ({kind}: {e})") from e
    if not isinstance(d, dict) or d.get("format") != COUNTS_FORMAT:
        got = d.get("format") if isinstance(d, dict) else type(d).__name__
        raise CountsError(f"{path}: not an {COUNTS_FORMAT} file (format={got!r})")
    return d


def read_fasta(path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return {k: v.upper() for k, v in seqs.items()}


@dataclass
class PanelClass:
    name: str
    kind: str
    length: int
    circular: bool
    kmer_pos: np.ndarray  # unit start positions of the retained k-mers (positional classes)


@dataclass
class Panel:
    k: int
    classes: dict[str, PanelClass]


def load_panel(path) -> Panel:
    k, defs, pos = 0, {}, {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##k="):
                k = int(line[4:])
            elif line.startswith("##class\t"):
                kv = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:])
                defs[int(kv["id"])] = kv
                pos[int(kv["id"])] = []
            elif not line.startswith("#"):
                p = line.split("\t")
                pos[int(p[1])].append(int(p[2]))
    classes = {}
    for i, kv in defs.items():
        if kv["name"] in classes:
            raise ValueError(f"{path}: class name {kv['name']!r} is defined twice")
        classes[kv["name"]] = PanelClass(kv["name"], kv["kind"], int(kv["length"]), kv["circular"] == "1",
                                         np.array(sorted(pos[i]), dtype=np.int64))
    return Panel(k, classes)


def read_features(path) -> dict[str, list[tuple[str, int, int]]]:
    """Class feature table: BED-like `class start end name` in unit coordinates."""
    out: dict[str, list[tuple[str, int, int]]] = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c, s, e, n = line.rstrip("\n").split("\t")[:4]
            out.setdefault(c, []).append((n, int(s), int(e)))
    return out


class FastaIndex:
    """Random access to an uncompressed, faidx-indexed FASTA (resource-building scripts only)."""

    def __init__(self, path):
        self.path = str(path)
        self.index = {}
        with open(self.path + ".fai") as fh:
            for line in fh:
                n, ln, off, lb, lw = line.split("\t")[:5]
                self.index[n] = (int(ln), int(off), int(lb), int(lw))
        self._fh = open(self.path, "rb")

    def length(self, name) -> int:
        return self.index[name][0]

    def fetch(self, name, start, end) -> str:
        ln, off, lb, lw = self.index[name]
        start, end = max(0, start), min(ln, end)
        if end <= start:
            return ""
        first = off + (start // lb) * lw + start % lb
        last = off + ((end - 1) // lb) * lw + (end - 1) % lb
        self._fh.seek(first)
        raw = self._fh.read(last - first + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()
