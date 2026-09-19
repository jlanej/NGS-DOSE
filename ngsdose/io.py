"""Readers for the engine's counts JSON, k-mer panels and FASTA resources."""
from __future__ import annotations

import gzip
import json
from dataclasses import dataclass

import numpy as np

COUNTS_FORMAT = "ngs-dose-counts/1"


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_counts(path) -> dict:
    """Load one per-sample counts file written by `ngs-dose count`."""
    with _open(path) as fh:
        d = json.load(fh)
    if d.get("format") != COUNTS_FORMAT:
        raise ValueError(f"{path}: not an {COUNTS_FORMAT} file (format={d.get('format')!r})")
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
