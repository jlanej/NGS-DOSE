"""Resource bundles: everything that is specific to one reference build.

A bundle is a directory holding `bundle.json` and the files it names: the k-mer panel, the
control-region FASTA, the class sinks used by fetch mode, the unit sequences and feature tables
of positional classes, and (optionally) a window-efficiency table learned from a cohort.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

from . import estimate, io


def default_bundle() -> Path:
    env = os.environ.get("NGSDOSE_RESOURCES")
    if env:
        return Path(env)
    return Path(__file__).resolve().parent.parent / "resources" / "GRCh38"


class Bundle:
    def __init__(self, path=None):
        self.dir = Path(path) if path else default_bundle()
        meta = self.dir / "bundle.json"
        if not meta.exists():
            raise FileNotFoundError(f"{meta} not found: pass --resources or set NGSDOSE_RESOURCES")
        self.meta = json.loads(meta.read_text())

    def _p(self, key) -> Path:
        return self.dir / self.meta[key]

    @property
    def panel(self) -> Path:
        return self._p("panel")

    @property
    def controls(self) -> Path:
        return self._p("controls")

    @property
    def sinks(self) -> Path:
        return self._p("sinks")

    def units(self) -> dict[str, str]:
        out = {}
        for cls, rel in self.meta.get("units", {}).items():
            seqs = io.read_fasta(self.dir / rel)
            if len(seqs) != 1:
                raise ValueError(f"{rel}: a unit FASTA must hold exactly one record")
            out[cls] = next(iter(seqs.values()))
        return out

    def features(self) -> dict:
        rel = self.meta.get("features")
        return io.read_features(self.dir / rel) if rel else {}

    def anchors(self) -> dict[str, list[tuple[int, int]]]:
        """Per class, the unit intervals that set the absolute level (empty: use the GC rule)."""
        rel = self.meta.get("anchors")
        if not rel or not (self.dir / rel).exists():
            return {}
        tab = json.loads((self.dir / rel).read_text())
        return {cls: [tuple(iv) for iv in v["intervals"]] for cls, v in tab.items()}

    def contig_lengths(self) -> dict[str, int]:
        """Lengths of the build's primary contigs: what tells GRCh38 from a look-alike (hg19 has the same names)."""
        return {k: int(v) for k, v in self.meta.get("contig_lengths", {}).items()}

    def efficiencies(self) -> Path | None:
        rel = self.meta.get("efficiencies")
        return self.dir / rel if rel and (self.dir / rel).exists() else None

    def regions(self) -> list[tuple[str, str]]:
        """(name, role) of every region of the controls file, in file order (the row order of `region_tables`)."""
        if not hasattr(self, "_regions"):
            out = []
            with io._open(self.controls) as fh:
                for line in fh:
                    if line.startswith(">"):
                        f = line[1:].split()
                        out.append((f[0], next((t[5:] for t in f[1:] if t.startswith("role=")), "control")))
            self._regions = out
        return self._regions

    def region_tables(self, L: int) -> np.ndarray:
        """Per-control-region GC tables for window length L, cached beside the user's cache dir."""
        st = self.controls.stat()
        key = hashlib.sha1(f"{self.controls.resolve()}:{st.st_size}:{int(st.st_mtime)}:{L}".encode()).hexdigest()[:16]
        cache = Path(os.environ.get("NGSDOSE_CACHE", Path.home() / ".cache" / "ngsdose"))
        f = cache / f"region_tables.{key}.npy"
        if f.exists():
            return np.load(f)
        tab = estimate.control_region_tables(self.controls, L)
        try:
            cache.mkdir(parents=True, exist_ok=True)
            np.save(f, tab)
        except OSError:
            pass
        return tab
