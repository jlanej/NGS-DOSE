"""Resource bundles: everything that is specific to one reference build.

A bundle is a directory holding `bundle.json` and the files it names: the k-mer panel, the
control-region FASTA, the class sinks used by fetch mode, and the unit sequences, feature tables
and anchor intervals of positional classes. (A cohort's window-efficiency table is kept outside
the bundle and applied with `ngsdose cohort --efficiencies`.) Beside them it may hold named
subsets of its control regions (`controls.<name>.bed`, Bundle.control_subsets), and units of
positional classes from experimental panels are looked for outside it (ExperimentalUnits).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

from . import estimate, io

EXTRA_UNITS_ENV = "NGSDOSE_EXTRA_UNITS"
UNIT_SUFFIXES = (".fa", ".fa.gz", ".fasta", ".fasta.gz")


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

    def control_subsets(self) -> dict[str, frozenset]:
        """Named subsets of the bundle's control regions, which counts may have been made with instead
        of all of them (a fetch with a lighter controls file: fewer bytes read, a GC curve fitted on
        fewer regions). Each is a BED in the bundle directory named `controls.<name>.bed`, or one that
        bundle.json lists under `control_subsets` ({name: file}), in the bundle BED's format: its
        control lines (name column empty or `control`) are the subset; `test:`/`dosage:` lines may be
        there too and do not count. The engine's controls FASTA for it is built from the BED with
        `ngs-dose controls` (by convention `controls.<name>.fa.gz`). Every region of a subset must be
        a control region of the bundle, under the same coordinates.

        The GRCh38 bundle ships `lite200`: 200 of the 800 control regions (at least three on every
        chromosome that has them, spread over each chromosome's GC range, the draw whose GC
        distribution is closest to all 800's) and all 182 known-truth and dosage regions. On the
        1000 Genomes CRAM of HG00096 its controls file reads 85 MB instead of 213 MB. Against all
        800, on the 1,748 cohort scans (GC tables rebuilt per region from the counts; the rebuild
        was checked against engine runs with subset controls files): the 45S rDNA headline moves by
        a median of +0.2% (SD 0.3%, at most 1.1%), 5S rDNA SD 0.7%, DJ SD 0.3%, truth.auto SD 0.2%,
        the satellite masses SD 0.3-0.8% except HSat1B (+1.5%, SD 2.3%); trio reliabilities and the
        pilot's cross-library agreement are unchanged. The all-window 45S estimate rises 1.5%,
        because a narrower GC support leaves fewer GC-rich windows usable (usable fraction 0.85 ->
        0.79); the chromosome test flags 20 of the 32 samples it flags with all 800 (on the same
        rebuilt tables), and no other."""
        if hasattr(self, "_subsets"):
            return self._subsets
        files = {p.name[len("controls."):-len(".bed")]: p for p in sorted(self.dir.glob("controls.*.bed"))
                 if p.name != Path(self.meta.get("controls_bed", "controls.bed")).name}
        files.update({n: self.dir / rel for n, rel in self.meta.get("control_subsets", {}).items()})
        bundle = {n for n, role in self.regions() if role == "control"}
        out = {}
        for name, f in files.items():
            if not f.exists():
                raise FileNotFoundError(f"{f}: the control subset {name} named in bundle.json is missing")
            names = set()
            with io._open(f) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    if len(p) < 3:
                        raise ValueError(f"{f}: not a BED line: {line.strip()[:60]!r}")
                    if len(p) < 4 or p[3] in ("", "control"):
                        names.add(f"{p[0]}:{int(p[1])}-{int(p[2])}")
            extra = sorted(names - bundle, key=estimate._natural)
            if extra or not names:
                raise ValueError(f"{f}: a control subset must name control regions of the bundle, and "
                                 + (f"{len(extra)} of its {len(names)} are not ({', '.join(extra[:3])})" if extra else "it names none"))
            out[name] = frozenset(names)
        self._subsets = out
        return out

    def experimental(self) -> "ExperimentalUnits":
        """Units of positional classes the bundle does not carry: see ExperimentalUnits. One object
        per bundle, so that the contract check and the estimate share what it has read. Paths are
        reported relative to the directory that holds the bundle's `resources` directory (the
        repository, or the image's install root), so that a host and a container say the same."""
        if not hasattr(self, "_experimental"):
            self._experimental = ExperimentalUnits(experimental_unit_dirs(self.dir), root=self.dir.resolve().parent.parent)
        return self._experimental

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


def experimental_unit_dirs(bundle_dir) -> list[Path]:
    """Where units of experimental positional classes are looked for, in order: the repository's
    `resources/experimental/candidates/units` (beside the bundle's directory), then every directory
    named by NGSDOSE_EXTRA_UNITS (several separated by the path separator, ':')."""
    dirs = [Path(bundle_dir).resolve().parent / "experimental" / "candidates" / "units"]
    dirs += [Path(d) for d in os.environ.get(EXTRA_UNITS_ENV, "").split(os.pathsep) if d]
    return dirs


class ExperimentalUnits:
    """Units and panel entries of positional classes from experimental panels.

    A counts file can carry a positional class the bundle has neither a panel entry nor a unit for
    (a scan made with an extra panel, `ngs-dose count -p`). The estimator needs two things to
    measure it: the class's k-mer positions on its unit (which read starts are callable) and the
    unit's sequence (the expected count of each window). The unit is `<class>.fa` (or .fa.gz,
    .fasta, .fasta.gz), one record, in the first of `dirs` that has it; the panel entry comes from
    a panel file (`*.tsv`, `*.tsv.gz` starting `##ngs-dose-panel`) in the unit's directory or the
    one above it, and when the counts record the sha256 of the panels they were made with, only a
    panel among them is used (a panel rebuilt since would put the k-mers elsewhere). Without a
    unit the class stays skipped, as before, and its reason says where a unit was looked for."""

    def __init__(self, dirs, root=None):
        self.dirs = [Path(d) for d in dirs]
        self.root = Path(root).resolve() if root else None
        self._seqs: dict[Path, dict[str, str]] = {}
        self._headers: dict[Path, dict[str, dict]] = {}
        self._hashes: dict[Path, set[str]] = {}
        self._panels: dict[Path, io.Panel] = {}

    def unit_file(self, name: str) -> Path | None:
        for d in self.dirs:
            for suf in UNIT_SUFFIXES:
                if (d / f"{name}{suf}").is_file():
                    return d / f"{name}{suf}"
        return None

    def show(self, path) -> str:
        """A path as results and messages give it: relative to `root` when it lies under it."""
        if self.root is not None:
            try:
                return str(Path(path).resolve().relative_to(self.root))
            except ValueError:
                pass
        return str(path)

    def missing(self, name: str) -> str:
        """What a skipped class's reason says when no unit was found: which files, where."""
        return f"no experimental unit {name}.fa in {', '.join(self.show(d) for d in self.dirs) or 'no directory'}"

    def _header(self, path: Path) -> dict[str, dict]:
        """name -> class definition of a panel file, from its header lines only; {} if it is not a panel."""
        if path not in self._headers:
            out = {}
            try:
                with io._open(path) as fh:
                    first = fh.readline()
                    if first.startswith("##ngs-dose-panel"):
                        for line in fh:
                            if not line.startswith("##"):
                                break
                            if line.startswith("##class\t"):
                                kv = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:])
                                out[kv["name"]] = kv
            except (OSError, UnicodeDecodeError, EOFError, ValueError):
                out = {}
            self._headers[path] = out
        return self._headers[path]

    def lookup(self, name: str, counts: dict, pc: io.PanelClass | None = None):
        """(PanelClass, unit sequence, unit path as `show` gives it) for a positional class; a string
        saying what is missing when a unit exists but cannot be used; None when there is no unit at
        all (`missing` says where it was looked for). `pc`: the
        bundle panel's entry for the class, when it has one (only the unit is then looked for)."""
        from .contract import panel_hashes
        f = self.unit_file(name)
        if f is None:
            return None
        if f not in self._seqs:
            self._seqs[f] = io.read_fasta(f)
        seqs, shown = self._seqs[f], self.show(f)
        if len(seqs) != 1:
            return f"the experimental unit {shown} holds {len(seqs)} records, not one"
        seq = next(iter(seqs.values()))
        if pc is None:
            files = sorted({p for d in (f.parent, f.parent.parent) for pat in ("*.tsv", "*.tsv.gz") for p in d.glob(pat)})
            cands = [p for p in files if self._header(p).get(name, {}).get("kind") == "positional"]
            if not cands:
                return f"the experimental unit {shown} has no panel beside it that defines {name} as positional"
            recorded = counts.get("panel_sha256")
            recorded = {recorded} if isinstance(recorded, str) else set(recorded or ())
            if recorded:
                for p in cands:
                    if p not in self._hashes:
                        self._hashes[p] = panel_hashes(p)
                cands = [p for p in cands if self._hashes[p] & recorded]
                if not cands:
                    return (f"no panel beside the experimental unit {shown} that defines {name} is one the counts were made with "
                            "(their panel_sha256)")
            p = cands[0]
            if p not in self._panels:
                self._panels[p] = io.load_panel(p)
            panel = self._panels[p]
            if counts.get("k") is not None and panel.k != counts["k"]:
                return f"the experimental panel {p.name} has k={panel.k}, the counts k={counts['k']}"
            pc = panel.classes[name]
        if len(seq) != pc.length:
            return f"the experimental unit {shown} has {len(seq)} bp, the panel says {pc.length}"
        return pc, seq, shown
