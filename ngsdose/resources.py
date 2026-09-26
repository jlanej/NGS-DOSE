"""Resource bundles: everything that is specific to one reference build.

A bundle is a directory holding `bundle.json` and the files it names: the k-mer panel, the
control-region FASTA, the class sinks used by fetch mode, and the unit sequences, feature tables
and anchor intervals of positional classes. (A cohort's window-efficiency table is kept outside
the bundle and applied with `ngsdose cohort --efficiencies`.)
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
            if not (self.dir / rel).exists():
                raise FileNotFoundError(f"{self.dir / rel}: the unit of {cls} named in bundle.json is missing")
            seqs = io.read_fasta(self.dir / rel)
            if len(seqs) != 1:
                raise ValueError(f"{rel}: a unit FASTA must hold exactly one record")
            out[cls] = next(iter(seqs.values()))
        return out

    def features(self) -> dict:
        rel = self.meta.get("features")
        return io.read_features(self.dir / rel) if rel else {}

    def check_units(self, panel: io.Panel, units: dict[str, str] | None = None):
        """Every positional class of the panel needs a unit of its length, and every unit a positional class."""
        units = self.units() if units is None else units
        bad = [f"{n}: no unit in bundle.json 'units'" for n, pc in panel.classes.items() if pc.kind == "positional" and n not in units]
        for n, seq in units.items():
            pc = panel.classes.get(n)
            if pc is None or pc.kind != "positional":
                bad.append(f"{n}: a unit in bundle.json but no positional class of the panel")
            elif len(seq) != pc.length:
                bad.append(f"{n}: the unit has {len(seq)} bp, the panel says {pc.length}")
        if bad:
            raise ValueError(f"{self.dir}: the bundle's units and panel disagree - " + "; ".join(bad))

    def anchors(self) -> dict[str, list[tuple[int, int]]]:
        """Per class, the unit intervals that set the absolute level (empty: use the GC rule)."""
        rel = self.meta.get("anchors")
        if not rel:
            return {}
        if not (self.dir / rel).exists():
            raise FileNotFoundError(f"{self.dir / rel} named in bundle.json is missing; restore it, or pass --gc-rule-anchors "
                                    "to set the level by the fragment-GC rule instead")
        tab = json.loads((self.dir / rel).read_text())
        return {cls: [tuple(iv) for iv in v["intervals"]] for cls, v in tab.items()}

    def contig_lengths(self) -> dict[str, int]:
        """Lengths of the build's primary contigs: what tells GRCh38 from a look-alike (hg19 has the same names)."""
        return {k: int(v) for k, v in self.meta.get("contig_lengths", {}).items()}

    def regions(self) -> list[tuple[str, str]]:
        """(name, role) of every region of the controls file, in file order (the row order of `region_tables`)."""
        if not hasattr(self, "_regions"):
            out = []
            with io._open(self.controls) as fh:
                for line in fh:
                    if line.startswith(">"):
                        f = line[1:].split()
                        tags = dict(t.split("=", 1) for t in reversed(f[1:]) if "=" in t)       # the first of a key, as the engine
                        role, label = tags.get("role", "control"), tags.get("label")
                        # as the engine: a misspelt role would otherwise turn a control into a truth region
                        if role not in ("control", "test", "dosage") or (role != "control" and not label):
                            raise ValueError(f"{self.controls}: record {f[0]} has role={role}"
                                             f"{'' if label is None else ' label=' + label}; allowed are control, "
                                             "test with a label= and dosage with a label=")
                        out.append((f[0], role))
            if not any(role == "control" for _, role in out):
                raise ValueError(f"{self.controls}: no region with role=control")
            self._regions = out
        return self._regions

    def region_tables(self, L: int) -> np.ndarray:
        """Per-control-region GC tables for window length L, cached beside the user's cache dir."""
        st = self.controls.stat()
        key = hashlib.sha1(f"{self.controls.resolve()}:{st.st_size}:{int(st.st_mtime)}:{L}".encode()).hexdigest()[:16]
        cache = Path(os.environ.get("NGSDOSE_CACHE", Path.home() / ".cache" / "ngsdose"))
        f = cache / f"region_tables.{key}.npy"
        try:
            tab = np.load(f)
            if tab.shape == (len(self.regions()), estimate.gcmodel.GC_BINS):
                return tab
        except (OSError, EOFError, ValueError):
            pass                                        # absent, or cut short by a killed or concurrent writer: recompute
        tab = estimate.control_region_tables(self.controls, L)
        tmp = cache / f"region_tables.{key}.{os.getpid()}.tmp.npy"
        try:
            cache.mkdir(parents=True, exist_ok=True)
            np.save(tmp, tab)
            os.replace(tmp, f)                          # readers see the old file or the whole new one, never part of it
        except OSError:
            tmp.unlink(missing_ok=True)
        return tab
