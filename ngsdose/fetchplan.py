"""`ngsdose fetchplan`: choose what a fetch reads, trading the bytes read against what is measured.

The menu (resources/fetch_menu.tsv; its header comment gives the format) lists the options: one row
per class, with the panel file that defines it and the sinks BED that says where its reads land,
plus a row for the control regions (always read), one for the unmapped bin, and rows for named
subsets of a class's learned intervals (kind `subset`: one array of a satellite family, such as
DXZ1 in aSatHOR's sinks). Sinks are learned from whole-file scans (`ngsdose sinks`), never assumed
from a class's reference coordinates: a class whose status is `candidate` has none yet, so it can
be selected for scanning only.

Classes are chosen by name, by preset, or by a byte budget; a capture target per class drops the
intervals that yield least while the expected capture (from `ngsdose sinks --stats`) stays at or
above it. Without CRAM indexes the intervals are kept in order of share per read (the statistics'
rank); with them (--crai) in order of share per byte of the interval's own fetch (the CRAM slices
it overlaps, with their containers' compression headers: what dropping it saves), and the expected
capture is then a lower bound (see _trim), unless keeping per read reaches the target with fewer
bytes. Classes share intervals (the chr2:32.91 Mb pile-up bin is in the sinks of ten), and the
engine counts every read of the plan's union for its class: a class keeps any interval the rest of
the plan reads anyway, and with CRAM indexes is trimmed again with the intervals that add no bytes
to the plan's fetches free (_share_between_options). With CRAM indexes the bytes each option and
the plan read are reported as the engine reads them (ngsdose.cost: one fetch per run of touching
or overlapping intervals, a slice under several runs decoded once per run), beside the floor of
each slice once, and what trimming saved of the plan and cost. The outputs are what `ngs-dose count -m fetch` takes:
PREFIX.sinks.bed (--sinks), PREFIX.panels.txt (one -p per line), PREFIX.count_flags.txt (further
flags, one token per line), PREFIX.controls.txt (the controls FASTA to pass as -c: the one matching
the control regions the plan was costed on), plus PREFIX.scan_panels.txt (the panels to load in the
whole-file scans, candidates included) and PREFIX.plan.tsv. The engine fetches the intervals of the
BED it is given.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from . import cost, sinks

COLUMNS = ("class", "group", "kind", "status", "tier", "presets", "panel", "sinks", "measures", "truth", "notes")
KINDS = ("controls", "unmapped", "positional", "compositional", "subset")
CLASS_KINDS = ("positional", "compositional", "subset")    # options that fetch sink intervals
STATUSES = ("shipped", "experimental", "candidate")
TIERS = ("A", "B", "C", "D")
ENGINE_PAD = 600       # `ngs-dose count --pad` default (src/count.rs DEFAULT_PAD)
READLEN_MAX = 400      # the engine refuses a smaller pad (src/count.rs READLEN_MAX)
PLAN_COLUMNS = ("option", "status", "tier", "kind", "intervals", "dropped", "bp", "capture_target", "expected_capture", "full_capture",
                "capture_lost", "mb_median", "mb_saved", "cum_mb_median", "cum_mb_floor", "cum_pct_median", "order", "note")


def default_menu() -> Path:
    """The menu beside the resource bundle ($NGSDOSE_RESOURCES/../fetch_menu.tsv, as in the container, where
    the package is installed away from resources/), else resources/fetch_menu.tsv of a source checkout."""
    env = os.environ.get("NGSDOSE_RESOURCES")
    if env:
        p = Path(env).parent / "fetch_menu.tsv"
        if p.exists():
            return p
        # not the source checkout's: its rows would name the checkout's panels and sinks, and the plan would mix two bundles
        raise ValueError(f"no fetch menu beside the resource bundle ({p}, NGSDOSE_RESOURCES={env}): pass one with --menu, or keep "
                         "resources/fetch_menu.tsv beside the bundle's directory as the repository, the image and the release tarball do")
    return Path(__file__).resolve().parent.parent / "resources" / "fetch_menu.tsv"


def controls_fasta(bed: Path) -> Path:
    """The controls FASTA that goes with a controls BED: controls.NAME.bed -> controls.NAME.fa.gz beside it,
    as the bundle names them (built by `ngs-dose controls -b BED -T REF --flank 1000`)."""
    name = bed.name[:-4] if bed.name.endswith(".bed") else bed.name
    return bed.with_name(name + ".fa.gz")


@dataclass
class Option:
    name: str
    group: str
    kind: str
    status: str
    status_note: str
    tier: str
    presets: tuple[str, ...]
    panel: str
    sinks: str
    measures: str = ""
    truth: str = ""
    notes: str = ""
    stats: str = ""
    order: int = 0
    parent: str = ""                                       # kind subset: the class whose learned intervals it names

    @property
    def counted(self) -> str:
        """The class the engine counts for this option: its own name, or a subset's class."""
        return self.parent or self.name


@dataclass
class Menu:
    path: Path
    options: dict[str, Option]
    presets: dict[str, str]                                # name -> description
    controls: Option
    unmapped: Option | None

    def resolve(self, rel: str) -> Path:
        return Path(rel) if Path(rel).is_absolute() else self.path.parent / rel

    def members(self, preset: str) -> list[str]:
        return [o.name for o in self.options.values() if preset in o.presets]


def _blank(x: str) -> str:
    return "" if x.strip() in ("", "-", ".") else x.strip()


def subset_parent(path) -> str:
    """The one class a subset BED names in its fourth column (every row must name it)."""
    rows = sinks.read_bed(path)
    names = {r[3] for r in rows}
    if not rows or len(names) != 1 or "" in names:
        raise ValueError(f"{path}: a subset's BED names one class in its fourth column on every row "
                         f"({'no rows' if not rows else 'it names ' + (', '.join(sorted(n or '(none)' for n in names)))})")
    return names.pop()


def read_menu(path=None) -> Menu:
    """Parse a fetch menu: '##preset<TAB>name<TAB>description' lines, other '#' lines are comments,
    then a header naming at least COLUMNS (in any order; a 'stats' column is optional) and one row per option."""
    path = Path(path) if path else default_menu()
    if not path.exists():
        raise ValueError(f"no fetch menu at {path}: pass it with --menu (or set NGSDOSE_RESOURCES to the bundle beside fetch_menu.tsv)")
    presets, rows, head = {}, [], None
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if line.startswith("##preset"):
                p = line.split("\t")
                if len(p) < 2 or not p[1].strip():
                    raise ValueError(f"{path}:{i}: a ##preset line needs a name: ##preset<TAB>name<TAB>description")
                presets[p[1].strip()] = p[2].strip() if len(p) > 2 else ""
                continue
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if head is None:
                head = p
                miss = [c for c in COLUMNS if c not in head]
                if miss:
                    raise ValueError(f"{path}:{i}: the header lacks the column(s) {', '.join(miss)}")
                continue
            if len(p) != len(head):
                raise ValueError(f"{path}:{i}: {len(p)} tab-separated fields, the header has {len(head)}")
            rows.append((i, dict(zip(head, p))))
    if head is None:
        raise ValueError(f"{path}: no header line")
    options: dict[str, Option] = {}
    menu = Menu(path, options, presets, None, None)       # resolve() only needs the path
    for i, r in rows:
        where = f"{path}:{i}"
        name = r["class"].strip()
        st, _, note = r["status"].partition(":")
        o = Option(name=name, group=r["group"].strip(), kind=r["kind"].strip(), status=st.strip(), status_note=note.strip(),
                   tier=r["tier"].strip(), presets=tuple(x.strip() for x in _blank(r["presets"]).split(",") if x.strip()),
                   panel=_blank(r["panel"]), sinks=_blank(r["sinks"]), measures=r["measures"].strip(), truth=_blank(r["truth"]),
                   notes=_blank(r["notes"]), stats=_blank(r.get("stats", "")), order=len(options))
        if not name or name in options:
            raise ValueError(f"{where}: {'an empty' if not name else 'a second'} class name {name!r}")
        if name == "union":
            raise ValueError(f"{where}: 'union' names the whole plan in the cost report; call the option something else")
        if o.kind not in KINDS:
            raise ValueError(f"{where}: kind {o.kind!r}; allowed are {', '.join(KINDS)}")
        if o.status not in STATUSES:
            raise ValueError(f"{where}: status {o.status!r}; allowed are {', '.join(STATUSES)} (optionally followed by ': a note')")
        if o.tier not in TIERS:
            raise ValueError(f"{where}: tier {o.tier!r}; allowed are {', '.join(TIERS)}")
        if o.kind == "controls" and (not o.sinks or o.status != "shipped"):
            raise ValueError(f"{where}: the controls row needs its BED of regions in the sinks column and status shipped")
        if o.kind in CLASS_KINDS:
            if not o.panel:
                raise ValueError(f"{where}: {name} has no panel file")
            if o.status != "candidate" and not o.sinks:
                raise ValueError(f"{where}: {name} is {o.status} but names no sinks file: a class without learned sinks is a candidate")
        if o.kind == "subset":
            if o.status == "candidate" or not o.sinks:
                raise ValueError(f"{where}: {name} is a subset: it names a BED of some of a class's learned intervals, so it needs that "
                                 "file in the sinks column and cannot be a candidate")
            if o.stats:
                raise ValueError(f"{where}: {name} is a subset, which is fetched whole: it takes no statistics")
            bed = menu.resolve(o.sinks)
            if not bed.exists():
                raise ValueError(f"{where}: {name}: its subset BED {bed} does not exist")
            o.parent = subset_parent(bed)
            if o.parent == name:
                raise ValueError(f"{where}: {name} is a subset of itself: its BED must name the class it is part of")
        options[name] = o
    special = {k: [o for o in options.values() if o.kind == k] for k in ("controls", "unmapped")}
    if len(special["controls"]) != 1 or len(special["unmapped"]) > 1:
        raise ValueError(f"{path}: needs exactly one row of kind controls and at most one of kind unmapped")
    for o in options.values():
        if o.kind == "subset" and o.parent in options and options[o.parent].kind == "subset":
            raise ValueError(f"{path}: {o.name} is a subset of {o.parent}, which is itself a subset")
        for p in o.presets:
            presets.setdefault(p, "")
    for p in presets:
        if not any(p in o.presets for o in options.values()):
            raise ValueError(f"{path}: preset {p} has no member (list it in the presets column of its rows)")
    menu.controls, menu.unmapped = special["controls"][0], (special["unmapped"] or [None])[0]
    return menu


def panel_classes(path) -> list[str]:
    """The class names a panel file defines (its ##class header lines)."""
    from .io import _open
    out = []
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##class\t"):
                kv = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:] if "=" in f)
                out.append(kv.get("name", ""))
    return out


_OPTION = re.compile(r"^\s*(?:-\w,\s*)?(--[A-Za-z][\w-]*)", re.M)


def probe_engine(engine) -> tuple[set[str] | None, str]:
    """The long options `ENGINE count --help` lists (a line of the help that starts with the option), and
    what was probed: (None, why) when no engine was given, found or could be run. Engines differ: this
    branch's takes --classes and --allow-missing-sinks, main's 7772e32 --allow-missing-sinks (added in
    645ae55) but not --classes, and fae1124, which runs the 1000 Genomes cohort, neither (nor --pad)."""
    if not engine:
        return None, "no engine was probed"
    exe = shutil.which(str(engine)) or (str(engine) if Path(engine).is_file() else None)
    if not exe:
        return None, f"no engine {engine} was found to probe"
    try:
        r = subprocess.run([exe, "count", "--help"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError) as e:
        return None, f"{engine} count --help could not be run ({e})"
    if r.returncode != 0:
        return None, f"{engine} count --help exited with status {r.returncode}"
    return set(_OPTION.findall(r.stdout + r.stderr)), f"{engine} count --help"


@dataclass
class Row:
    option: Option
    intervals: list[tuple[str, int, int]] = field(default_factory=list)
    dropped: int = 0
    target: float | None = None
    expected: float | None = None
    full: float | None = None
    mb: float | None = None                                # the option's own intervals alone, as the engine fetches them (cost.price)
    mb_full: float | None = None                           # all of the option's intervals, before a capture target trimmed them
    mb_saved: float | None = None                          # what trimming the option saved of the plan (the rest of the plan as it is)
    whole: list[tuple[str, int, int]] = field(default_factory=list)   # all of a class's intervals, before trimming
    stats: list[dict] | None = None                        # the statistics rows trimming used
    per_byte: bool = False                                 # the statistics allow the byte order (share_max)
    shared: int = 0                                        # intervals below the target kept because they add no bytes to the plan's fetches
    cum_mb: float | None = None                            # the plan up to this option, as the engine fetches it
    cum_mb_floor: float | None = None                      # the same with every slice decoded once (cost.price_once)
    cum_pct: float | None = None
    order: str = ""                                        # 'read' or 'byte': the order capture trimming kept intervals in
    held_out: bool | None = None                           # what the statistics file says of its scans
    stats_scans: int | None = None
    note: str = ""
    warn: list[str] = field(default_factory=list)          # notes that matter only if the option is in the plan


@dataclass
class Plan:
    rows: list[Row]
    scan_only: list[Option]
    panels: list[str]
    scan_panels: list[str]
    flags: list[str]
    notes: list[str]
    n_index: int = 0
    controls_fasta: str = ""                              # the -c FASTA matching the control regions costed
    classes_flag: str = ""                                # what count_flags.txt says of unselected classes, and why
    capture_note: str = ""                                # what the expected capture is (plan.tsv header)

    def bed(self) -> list[tuple[str, int, int, str]]:
        """The sinks BED: every interval under the class the engine counts (a subset's under its class), once."""
        return sorted({(c, s, e, r.option.counted) for r in self.rows if r.option.kind in CLASS_KINDS for c, s, e in r.intervals})


def _load_bed_rows(paths) -> dict[str, list[tuple[str, int, int]]]:
    by = defaultdict(list)
    for p in paths:
        rows = sinks.read_bed(p)
        blank = sum(1 for r in rows if not r[3])
        if blank:
            raise ValueError(f"{p}: {blank} intervals have no class column; a fetch plan picks intervals by class")
        for c, s, e, n in rows:
            by[n].append((c, s, e))
    return by


def _matches(stats, ivs) -> str:
    """'' when the statistics are for exactly these intervals, else what differs."""
    have = {(r["contig"], r["start"], r["end"]) for r in stats}
    return "" if have == set(ivs) else f"{len(have & set(ivs))} of {len(set(ivs))} intervals match"


def _key(r) -> tuple[str, int, int]:
    return (r["contig"], r["start"], r["end"])


def _bounds(stats, order, col, full) -> list[float]:
    """The expected capture of keeping the first i + 1 rows of `order` (statistics rows in any order), for
    every i: a lower bound that holds scan by scan, the larger of the capture of all the intervals less the
    share_max of each one not kept (when the statistics have share_max), and the statistics' curve up to
    the first of their ranks not kept. For a run of the statistics' own order that is the curve, exact."""
    has_max = all(r.get("share_max") is not None for r in stats)
    after, run = [0.0] * len(order), 0.0
    for i in range(len(order) - 1, -1, -1):
        after[i] = run
        run += order[i]["share_max"] if has_max else 0.0
    held, prefix, out = set(), 0, []
    for i, r in enumerate(order):
        held.add(r["rank"])
        while prefix < len(stats) and stats[prefix]["rank"] in held:
            prefix += 1
        curve = stats[prefix - 1][col] if prefix else 0.0
        out.append(full if len(held) == len(stats) else max(full - after[i], curve) if has_max else curve)
    return out


def _capture_of(stats, kept, stat) -> float | None:
    """The expected capture of the intervals `kept` of a class (as _bounds: exact for a run of the
    statistics' own order, else a lower bound)."""
    if not stats:
        return None
    col = f"cum_capture_{stat}"
    order = [r for r in stats if _key(r) in kept]
    rest = [r for r in stats if _key(r) not in kept]
    return _bounds(stats, order + rest, col, stats[-1][col])[len(order) - 1] if order else 0.0


def _trim(ivs, stats, target, stat, nbytes=None, price_mb=None, free=frozenset()):
    """(kept intervals, dropped, expected capture, capture of all, order) of one class, by its stats rows.

    Without `nbytes` the intervals are kept in the statistics' rank (share per read) and the expected
    capture is the statistics' cumulative capture, exact for those scans. With `nbytes` ((contig,
    start, end) -> bytes of the interval's own fetch: the CRAM slices it overlaps with their
    containers' compression headers, what dropping it saves) they are kept in order of median share
    per byte: an interval of low median share on costly slices (a pile-up bin, a multi-reference decoy
    slice) ranks behind the intervals worth their bytes. The statistics hold the capture curve of their
    own order only, so the expected capture of a set kept per byte is a lower bound that holds scan by
    scan (_bounds). An interval is therefore dropped only while the largest share it held in any one
    scan still leaves the target reached: TEL keeps the chr2:32.91 Mb pile-up bin above about 0.979.
    `price_mb(intervals)` prices a set as the engine fetches it: of the two orders, the one that
    reaches the target with fewer bytes is kept (order 'byte', or 'read, fewer bytes' when the
    statistics' own order does, its capture then exact). `free`: intervals that add no bytes to the
    rest of the plan's fetches; they are kept first, in either order, and cost nothing."""
    if not stats:
        return ivs, 0, None, None, ""
    col = f"cum_capture_{stat}"
    full = stats[-1][col]

    def pick(order, how):
        if target is None:
            return sorted(_key(r) for r in order), 0, full, full, how
        b = _bounds(stats, order, col, full)
        j = next((i + 1 for i in range(len(order)) if b[i] >= target), len(order))
        return sorted(_key(r) for r in order[:j]), len(order) - j, full if j == len(order) else b[j - 1], full, how

    if free:
        by_read = pick([r for r in stats if _key(r) in free] + [r for r in stats if _key(r) not in free], "read")
    else:
        k = len(stats)
        if target is not None:
            k = next((i + 1 for i, r in enumerate(stats) if r[col] >= target), len(stats))
        by_read = (sorted(_key(r) for r in stats[:k]), len(stats) - k, stats[k - 1][col], full, "read")
    if nbytes is None:
        return by_read
    ranked = sorted(stats, key=lambda r: (_key(r) not in free, -r["share_median"] / max(nbytes[_key(r)], 1), -r["share_median"], r["rank"]))
    by_byte = pick(ranked, "byte")
    if price_mb is not None and price_mb(by_read[0]) < price_mb(by_byte[0]):
        return (*by_read[:4], "read, fewer bytes")
    return by_byte


def _inside_all(sub, ivs) -> list[tuple[str, int, int]]:
    """The intervals of `sub` not wholly inside the (merged) intervals `ivs`."""
    index = sinks._index([(c, s, e, "x") for c, s, e in ivs])
    return [(c, s, e) for c, s, e in sub if not sinks._inside(index, "x", c, s, e)]


def make_plan(menu: Menu, classes=(), presets=(), budget_mb=None, fill=False, statuses=("shipped", "experimental"), capture=None, capture_class=None,
              capture_stat="p10", sinks_files=(), stats_files=(), controls=None, pad=600, crais=(), contigs=None, engine=None,
              unmarked_companions=False, log=None) -> Plan:
    """Select options and cost them; see the module docstring and `ngsdose fetchplan --help`. `engine`: the
    `ngs-dose` that will run the fetch, whose `count --help` says which of the plan's count flags it takes
    (probe_engine); a plan it cannot run is refused. `unmarked_companions`: go on, with a warning, when the
    loaded panels define classes the plan does not select and the engine can neither leave them out
    (--classes) nor mark them (--allow-missing-sinks)."""
    log = log or (lambda m: print(m, file=sys.stderr))
    capture_class = dict(capture_class or {})
    if pad < READLEN_MAX:
        raise ValueError(f"--pad {pad} is below {READLEN_MAX} bp: reverse-strand 5' ends at the start of control regions would be missed, "
                         "and `ngs-dose count` refuses it")
    notes: list[str] = []
    chosen: list[str] = []
    for p in presets:
        if p not in menu.presets:
            raise ValueError(f"no preset {p} in {menu.path} (presets: {', '.join(sorted(menu.presets))})")
        chosen += menu.members(p)
    for c in classes:
        if c not in menu.options:
            raise ValueError(f"no option {c} in {menu.path} (options: {', '.join(menu.options)})")
        chosen.append(c)
    chosen = [c for c in dict.fromkeys(chosen) if menu.options[c].kind != "controls"]
    for c in capture_class:
        if c not in menu.options:
            raise ValueError(f"--capture-class names {c}, which is not in {menu.path}")
        if menu.options[c].kind == "subset":
            raise ValueError(f"--capture-class {c}: {c} is a named subset of {menu.options[c].parent}'s intervals, fetched whole; "
                             "a capture target applies to a class")
    if budget_mb is None and not chosen:
        raise ValueError("choose options with --classes, --preset or --budget-mb (--list shows the menu)")
    if budget_mb is not None and not crais:
        raise ValueError("--budget-mb needs the CRAM index(es) to cost the options: pass --crai (and --contigs)")
    if crais and not contigs:
        raise ValueError("--crai needs --contigs: the .fai or .dict of the reference the CRAMs were aligned to, or `samtools view -H` of one "
                         "(the index numbers contigs by their place in the header)")
    pool = [menu.options[c] for c in chosen] if chosen else [o for o in menu.options.values() if o.kind != "controls"]
    if budget_mb is not None:
        bad = [s for s in statuses if s not in STATUSES]
        if bad:
            raise ValueError(f"--status {', '.join(bad)}: allowed are {', '.join(STATUSES)}")
        if not chosen and "candidate" not in statuses and any(o.status == "candidate" for o in pool):
            notes.append(f"{sum(o.status == 'candidate' for o in pool)} candidate option(s) are left out: they have no learned sinks yet")
        pool = [o for o in pool if o.status in statuses or o.status == "candidate" and o.name in chosen]
    scan_only = [o for o in pool if o.status == "candidate"]
    fetch = [o for o in pool if o.status != "candidate"]
    admitted = [o.name for o in scan_only if o.name not in chosen]
    if admitted:
        notes.append(f"{len(admitted)} candidate option(s) admitted by --status have no learned sinks: a fetch cannot measure them, so "
                     "they are not fetched; their panels are in scan_panels.txt, to load in the whole-file scans where their sinks "
                     "are learned")
    for o in scan_only:
        if o.name in chosen:
            notes.append(f"{o.name} is a candidate: no sinks have been learned for it, so a fetch cannot measure it. Its panel is in "
                         "scan_panels.txt: load it in the whole-file scans, learn its sinks there (`ngsdose sinks SCANS"
                         f"{' --classes ' + o.name if o.kind == 'compositional' else ''}`), check their capture on held-out scans "
                         "(`ngsdose sinks HELD_OUT --evaluate BED`), then give it the sinks file in the menu")

    # the CRAM indexes: bytes per option, and per interval for the byte order of capture trimming
    indexes = []
    if crais:
        ctg = cost.read_contigs(contigs)
        indexes = [cost.CraiIndex(p, ctg) for p in crais]

    def nbytes(ivs):
        """(contig, start, end) -> median bytes of the interval's own fetch over the indexes."""
        return {iv: cost.median([ix.fetch(*iv) for ix in indexes]) for iv in ivs}

    def price_mb(ivs):
        """What the engine reads for a set of intervals (median over the indexes), in MB."""
        return cost.median([ix.price(ivs) for ix in indexes]) / 1e6

    # intervals of every fetchable option in the pool
    override = _load_bed_rows(sinks_files) if sinks_files else None
    cache: dict[Path, dict] = {}
    stats: dict[str, list[dict]] = {}
    source: dict[str, Path] = {}                          # class -> the statistics file its rows came from
    shipped = [menu.resolve(o.stats) for o in fetch if o.stats]
    for p in [q for q in dict.fromkeys(shipped) if q.exists()] + [Path(q) for q in stats_files]:
        got = sinks.read_stats(p)
        stats.update(got)                                  # a later file (the command line's) wins
        source.update(dict.fromkeys(got, p))
    held = {p: sinks.stats_held_out(p) for p in set(source.values())}
    rows: list[Row] = []
    menu_ctrl = menu.resolve(menu.controls.sinks)
    ctrl_bed = Path(controls) if controls else menu_ctrl
    ctrl_fa = controls_fasta(ctrl_bed)
    if not ctrl_fa.exists():
        notes.append(f"no controls FASTA {ctrl_fa} beside {ctrl_bed}: build it with `ngs-dose controls -b {ctrl_bed} -T REFERENCE "
                     f"--flank 1000 -o {ctrl_fa}`; the fetch must use the FASTA of exactly these regions")
    if ctrl_bed.resolve() != menu_ctrl.resolve():
        notes.append(f"the plan is costed on the control regions of {ctrl_bed.name}, not the menu's {menu_ctrl.name}: the fetch must pass "
                     f"-c {ctrl_fa.name} (written to controls.txt); a fetch with another controls file reads other regions than costed "
                     "here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names")
    ctrl = Row(menu.controls, [(c, max(0, s - pad), e + pad) for c, s, e, _ in sinks.read_bed(ctrl_bed)])
    ctrl.note = f"{len(ctrl.intervals)} regions of {ctrl_bed.name}, padded by {pad} bp as the engine reads them"
    for o in fetch:
        r = Row(o)
        if o.kind == "unmapped":
            r.intervals = [(cost.UNMAPPED, 0, 0)]
            r.note = "every read without a coordinate (count --unmapped)"
        elif o.kind == "subset":
            p = menu.resolve(o.sinks)
            ivs = sorted({(c, s, e) for c, s, e, _ in sinks.read_bed(p)})
            # a subset names some of its class's learned intervals: it holds only where the sinks the plan takes for the class
            # (--sinks, or the class's row of the menu) hold them too - not with another pipeline's sinks
            parent = menu.options.get(o.parent)
            if override is not None:
                have, src = override.get(o.parent, []), ", ".join(map(str, sinks_files))
            elif parent is not None and parent.sinks and menu.resolve(parent.sinks).exists():
                q = menu.resolve(parent.sinks)
                if q not in cache:
                    cache[q] = _load_bed_rows([q])
                have, src = cache[q].get(o.parent, []), str(q)
            else:
                have, src = None, ""
            if have is not None:
                out = _inside_all(ivs, have)
                if out:
                    raise ValueError(f"{o.name}: {len(out)} of its {len(ivs)} intervals (in {p.name}) are not inside {o.parent}'s intervals of "
                                     f"{src}: a subset is part of its class's learned sinks, so define it from the sinks this plan uses")
            r.intervals = ivs
            r.note = (f"a named subset of {o.parent}'s learned intervals ({p.name}), fetched whole and counted as {o.parent}: "
                      f"`ngsdose estimate` reports {o.name}.mass_Mb")
        else:
            if override is not None:
                ivs = override.get(o.name)
                src = ", ".join(map(str, sinks_files))
            else:
                p = menu.resolve(o.sinks)
                if not p.exists():
                    raise ValueError(f"{o.name}: its sinks file {p} does not exist. Sinks are learned from whole-file scans "
                                     f"(`ngsdose sinks SCANS{' --classes ' + o.name if o.kind == 'compositional' else ''} -o BED`); "
                                     "pass the BED with --sinks")
                if p not in cache:
                    cache[p] = _load_bed_rows([p])
                ivs, src = cache[p].get(o.name), str(p)
            if not ivs:
                raise ValueError(f"{o.name}: no interval of the class in {src}")
            target = capture_class.get(o.name, capture)
            st = stats.get(o.name)
            differ = _matches(st, ivs) if st else ""
            if differ and target is not None:
                raise ValueError(f"the statistics of {o.name} are for other intervals than its sinks ({differ}): compute them for this BED "
                                 "with `ngsdose sinks SCANS --evaluate BED --stats FILE`")
            if differ:
                r.warn.append(f"{o.name}: its statistics are for other intervals than its sinks ({differ}): no expected capture")
                st = None
            if target is not None and not st:
                if o.name in capture_class:
                    raise ValueError(f"--capture-class {o.name}={target}: no statistics for {o.name}; pass --stats "
                                     "(`ngsdose sinks SCANS --evaluate BED --stats FILE`)")
                r.warn.append(f"{o.name}: no statistics, so all its {len(ivs)} intervals are kept")
                target = None
            by_byte = None
            if st and indexes and target is not None:
                if all(x.get("share_max") is not None for x in st):
                    by_byte = nbytes([(x["contig"], x["start"], x["end"]) for x in st])
                else:
                    r.warn.append(f"{o.name}: its statistics have no share_max (written before it was added), so its intervals are kept "
                                  "by share per read, not per byte; write them again with `ngsdose sinks SCANS --evaluate BED --stats FILE`")
            r.whole = sorted(set(ivs))
            r.intervals, r.dropped, r.expected, r.full, r.order = _trim(r.whole, st, target, capture_stat, by_byte, price_mb)
            r.target, r.order, r.stats, r.per_byte = target, (r.order if target is not None else ""), st, by_byte is not None
            if st:
                r.held_out, r.stats_scans = held[source[o.name]], max(x["scans"] for x in st)
            if target is not None and r.full is not None and r.full < target:
                r.warn.append(f"{o.name}: all its intervals together capture {r.full:.4f} ({capture_stat} over the scans of its statistics), "
                             f"below the target {target}: all are kept")
            if r.dropped and indexes:
                r.mb_full = price_mb(r.whole)
        rows.append(r)

    # cost, and the budget: each option alone, then the plan as it grows, priced as the engine fetches it (all its intervals
    # merged together, so an option can even make the plan cheaper by bridging two runs that decode the same slices)
    absent = {}
    for r in [ctrl, *rows]:
        absent[id(r)] = cost.median([ix.keys(r.intervals)[1] for ix in indexes]) if indexes else 0
        r.mb = price_mb(r.intervals) if indexes else None
        if r.mb is not None and r.mb_full is None and r.target is not None:
            r.mb_full = r.mb                                # a capture target that dropped nothing saved nothing
    tier = {t: i for i, t in enumerate(TIERS)}
    rows.sort(key=lambda r: (tier[r.option.tier], r.mb if r.mb is not None else 0.0, r.option.order))
    seen: list[tuple[str, int, int]] = []
    plan = []
    for r in [ctrl, *rows]:
        if indexes:
            now = [ix.price(seen + r.intervals) for ix in indexes]
            mb = cost.median(now) / 1e6
            if budget_mb is not None and mb > budget_mb:
                if r is ctrl:
                    raise ValueError(f"the control regions alone read {mb:.1f} MB, more than the budget of {budget_mb:g} MB")
                if fill:
                    notes.append(f"{r.option.name} (tier {r.option.tier}) is left out: it would take the plan to {mb:.1f} MB, past the budget "
                                 f"of {budget_mb:g} MB")
                    continue
                notes.append(f"the budget of {budget_mb:g} MB stops the plan before {r.option.name} (tier {r.option.tier}), which would take it "
                             f"to {mb:.1f} MB; the options after it are left out (--fill adds those that still fit)")
                break
            seen = seen + r.intervals
            r.cum_mb, r.cum_pct = mb, 100 * cost.median([n / ix.total for n, ix in zip(now, indexes)])
            r.cum_mb_floor = cost.median([ix.price_once(seen) for ix in indexes]) / 1e6
        if absent[id(r)]:
            r.note = (r.note + "; " if r.note else "") + f"{absent[id(r)]:g} intervals on contigs the CRAM header lacks (median): not read"
        plan.append(r)
    if any(r.target is not None for r in plan):
        notes += _share_between_options(plan, indexes, capture_stat)
    for r in plan:
        notes += r.warn
    classes_in = [r.option for r in plan[1:] if r.option.kind in CLASS_KINDS]
    if not classes_in:
        raise ValueError("no class with sinks is selected (or the budget admits none): a fetch reads the sinks of at least one class")
    counted = list(dict.fromkeys(o.counted for o in classes_in))
    whole = {o.name for o in classes_in if o.kind != "subset"}
    parts = defaultdict(list)
    for o in classes_in:
        if o.kind == "subset":
            parts[o.parent].append(o.name)
    for parent, subs in parts.items():
        if parent in whole:
            notes.append(f"{', '.join(subs)} ({'a named subset' if len(subs) == 1 else 'named subsets'} of {parent}'s intervals) and {parent} are both fetched: "
                         "their intervals are read once")
        else:
            notes.append(f"{', '.join(subs)} {'is a named subset' if len(subs) == 1 else 'are named subsets'} of {parent}'s intervals, fetched "
                         f"without the rest of {parent}'s sinks: the fetch counts {parent} only inside them, so `ngsdose estimate` reports "
                         f"{', '.join(s + '.mass_Mb' for s in subs)} and marks {parent} itself as not measured (subset_only)")

    # panels, companions and flags
    panels, scan_panels, companions = [], [], []
    for o in classes_in + scan_only:
        p = menu.resolve(o.panel)
        if p.exists():
            names = panel_classes(p)
            if o.counted not in names:
                raise ValueError(f"{o.name}: its panel {p} defines {', '.join(names) or 'no class'}, not {o.counted}")
            if o in classes_in:
                companions += [n for n in names if n not in counted]
        elif o in classes_in:
            raise ValueError(f"{o.name}: its panel file {p} does not exist")
        else:
            notes.append(f"{o.name}: its panel file {p} does not exist yet")
        if o in classes_in and o.panel not in panels:
            panels.append(o.panel)
        if o.panel not in scan_panels:
            scan_panels.append(o.panel)
    flags = []
    if pad != ENGINE_PAD:
        flags.append(f"--pad={pad}")       # the fetch must read the control regions as costed here (one token per line)
    if any(r.option.kind == "unmapped" for r in plan):
        flags.append("--unmapped")
    companions = list(dict.fromkeys(companions))
    opts, probed = probe_engine(engine) if (flags or companions) else (None, "")
    if opts is not None:
        lack = [f for f in flags if f.split("=", 1)[0] not in opts]
        if lack:
            raise ValueError(f"the plan needs {', '.join(lack)} in count_flags.txt, which the engine does not take ({probed} lists no "
                             f"{', '.join(dict.fromkeys(f.split('=', 1)[0] for f in lack))}; fae1124 takes neither --pad nor --classes): "
                             "use the default --pad 600, or make the plan with --engine naming an engine that takes them")
    classes_flag = ""
    if companions:
        own = ", ".join(companions)
        if opts is None:
            # not probed: the flag engines since 645ae55 take (main's 7772e32 does). fae1124 refuses it and the fetch then fails
            # loudly (nothing is counted), never silently
            flags.append("--allow-missing-sinks")
            classes_flag = f"--allow-missing-sinks ({probed}; engines since 645ae55 take it, fae1124 does not)"
            notes.append(f"the panels loaded also define {own}, not selected. {probed}, so count_flags.txt has --allow-missing-sinks: "
                         "the counts file then lists those classes in sinks_missing_classes and `ngsdose estimate` marks them "
                         "incomplete. Engines since 645ae55 (main's 7772e32) take the flag; fae1124 does not and refuses the fetch; an engine that takes "
                         "--classes would count only the selected classes. Make the plan with --engine naming the engine that runs the fetch")
        elif "--classes" in opts:
            flags.append("--classes=" + ",".join(counted))
            classes_flag = f"--classes ({probed} lists --classes)"
            notes.append(f"the panels loaded also define {own}, not selected: count_flags.txt has --classes={','.join(counted)} "
                         f"({probed} lists it), so the fetch counts and reports only the selected classes and reads only their "
                         "sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming "
                         "the engine that runs the fetch")
        elif "--allow-missing-sinks" in opts:
            flags.append("--allow-missing-sinks")
            classes_flag = f"--allow-missing-sinks ({probed} lists it, not --classes)"
            notes.append(f"the panels loaded also define {own}, not selected: the fetch counts their reads only where they fall inside "
                         "the plan, the counts file lists them in sinks_missing_classes, and `ngsdose estimate` marks them incomplete. "
                         f"The engine refuses such a fetch unless --allow-missing-sinks (in count_flags.txt; {probed} lists it but not "
                         "--classes, which would count only the selected classes)")
        elif unmarked_companions:
            classes_flag = (f"none ({probed} lists neither --classes nor --allow-missing-sinks; --unmarked-companions: "
                            f"{own} are counted and not marked incomplete)")
            notes.append(f"WARNING: the panels loaded also define {own}, not selected, and {probed} lists neither --classes nor "
                         "--allow-missing-sinks (an engine such as fae1124). count_flags.txt has no flag for them: the fetch counts "
                         "their reads only where they fall inside the plan's intervals, an undercount, and the counts file does not say "
                         "so. `ngsdose estimate` marks them incomplete only when given this plan's sinks BED (--fetch-sinks "
                         "PREFIX.sinks.bed); without it their values are undercounts presented as measurements")
        else:
            raise ValueError(f"the panels loaded also define {own}, which the plan does not select, and the engine takes neither "
                             f"--classes nor --allow-missing-sinks ({probed} lists neither; fae1124 is such an engine). It would count "
                             "those classes only where their reads fall inside the plan's intervals and not mark them as incomplete. "
                             "Make the plan with --engine naming an engine that takes --classes, select those classes too, or pass "
                             "--unmarked-companions and give `ngsdose estimate` the plan's sinks BED (--fetch-sinks PREFIX.sinks.bed), "
                             "which then marks them incomplete")
    if any(o.status == "experimental" for o in classes_in):
        notes.append("experimental options selected (" + ", ".join(o.name for o in classes_in if o.status == "experimental") +
                     "): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet")
    capture_note = _capture_note([r for r in plan if r.expected is not None], capture_stat)
    if any(r.target is not None for r in plan):
        notes.append(capture_note)
    for n in notes:
        log(f"[fetchplan] {n}")
    return Plan(plan, scan_only, panels, scan_panels, flags, notes, len(indexes), str(ctrl_fa), classes_flag, capture_note)


def _holders(plan, r, ivs, own) -> list[str]:
    """The options of the plan other than `r` whose intervals (`own(row)`) overlap any of `ivs`."""
    return [x.option.name for x in plan if x is not r and any(c == c2 and s < e2 and s2 < e for c, s, e in ivs for c2, s2, e2 in own(x))]


def _share_between_options(plan: list[Row], indexes, stat) -> list[str]:
    """Trim again with the rest of the plan in view, and say what trimming saved of the plan.

    Classes share sink intervals (the chr2:32.91 Mb pile-up bin is in the sinks of TEL, rDNA45S and
    eight satellite families), and the engine counts every read of the union of the plan's intervals
    for whichever class it belongs to. An interval one class drops while another keeps it saves no
    byte and costs that class no read. The engine fetches the plan's intervals merged where they
    touch or overlap, so an interval's price in a plan is what its fetch adds to the plan's fetches
    (cost.CraiIndex.extra): nothing when it lies inside intervals of other options, its own fetch
    when it stands apart, less when it bridges two runs that decode the same slices. With CRAM
    indexes, each trimmed class is trimmed again with the intervals that add nothing to the rest of
    the plan's fetches (in every index) priced at nothing: those are kept first, and the class may
    then drop others (the set is changed only when the plan reads fewer bytes, or as many with a
    higher expected capture). Then, with or without indexes, every interval a class dropped that adds
    nothing to what the plan reads (with indexes: its fetch adds no bytes to the plan's, in every
    index; without, it lies inside the plan's other intervals) is given back, and the class's expected
    capture is recomputed. Each option's mb_saved is then what trimming it saved of the plan as it is
    (the plan with all the option's intervals, less the plan as planned, both priced as fetched), and
    the rows' own and cumulative bytes are recomputed."""
    notes: list[str] = []
    trimmed = [r for r in plan if r.dropped and r.stats]
    alone = {id(r): set(r.intervals) for r in trimmed}

    def rest_of(r) -> list[tuple[str, int, int]]:
        return [iv for x in plan if x is not r for iv in x.intervals]

    def adds_nothing(fetches, iv) -> bool:
        """The interval's fetch adds no bytes to `fetches` (from cost.merge), in every index."""
        return all(ix.extra(fetches, *iv) <= 0 for ix in indexes)

    if indexes:
        for r in trimmed:
            rest = rest_of(r)
            fetches = cost.merge(rest)
            extra = {iv: [ix.extra(fetches, *iv) for ix in indexes] for iv in r.whole}
            free = frozenset(iv for iv, b in extra.items() if all(x <= 0 for x in b))
            if not free - set(r.intervals):
                continue                                   # nothing the class dropped is read anyway: its own trim stands
            nb = {iv: cost.median(b) for iv, b in extra.items()}   # what each adds to the rest of the plan: at most its own fetch

            def total(ivs, rest=rest):
                return cost.median([ix.price(rest + list(ivs)) for ix in indexes])

            cur = sorted(set(r.intervals) | free)
            cur_exp = _capture_of(r.stats, set(cur), stat)
            got = _trim(r.whole, r.stats, r.target, stat, nb if r.per_byte else None, total if r.per_byte else None, free)
            if (total(got[0]), -got[2]) < (total(cur), -cur_exp):
                r.intervals, r.dropped, r.expected, r.order = got[0], got[1], got[2], got[4]
            else:
                r.intervals, r.dropped, r.expected = cur, len(r.whole) - len(cur), cur_exp
    for r in trimmed:
        kept = set(r.intervals)
        if indexes:
            # against the plan before any is given back: intervals that each add nothing add nothing together either
            fetches = cost.merge([iv for x in plan for iv in x.intervals])
            back = [iv for iv in r.whole if iv not in kept and adds_nothing(fetches, iv)]
        else:
            index = sinks._index([(c, s, e, "x") for x in plan for c, s, e in x.intervals if c != cost.UNMAPPED])
            back = [iv for iv in r.whole if iv not in kept and sinks._inside(index, "x", *iv)]
        if back:
            r.intervals = sorted(kept | set(back))
            r.dropped, r.expected = len(r.whole) - len(r.intervals), _capture_of(r.stats, set(r.intervals), stat)
    for r in trimmed:
        more = sorted(set(r.intervals) - alone[id(r)])
        if indexes and more:
            fetches = cost.merge(rest_of(r))
            more = [iv for iv in more if adds_nothing(fetches, iv)]
        if not more:
            continue

        # name the options that read the intervals of their own accord, before any interval was given back, if any do
        by = (_holders(plan, r, more, lambda x: alone[id(x)] & set(x.intervals) if id(x) in alone else x.intervals)
              or _holders(plan, r, more, lambda x: x.intervals))
        r.shared = len(more)
        r.note = (r.note + "; " if r.note else "") + (
            f"{len(more)} interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway"
            f"{' (for ' + ', '.join(by) + ')' if by else ''}, so they cost nothing and add to its capture")
    if not indexes:
        return notes
    seen: list[tuple[str, int, int]] = []
    for r in plan:
        r.mb = cost.median([ix.price(r.intervals) for ix in indexes]) / 1e6
        seen = seen + r.intervals
        now = [ix.price(seen) for ix in indexes]
        r.cum_mb, r.cum_pct = cost.median(now) / 1e6, 100 * cost.median([b / ix.total for b, ix in zip(now, indexes)])
        r.cum_mb_floor = cost.median([ix.price_once(seen) for ix in indexes]) / 1e6
    for r in plan:
        if r.target is not None:
            r.mb_saved = cost.median([ix.price(seen + r.whole) - ix.price(seen) for ix in indexes]) / 1e6
    whole = cost.median([ix.price([iv for r in plan for iv in (r.whole or r.intervals)]) for ix in indexes]) / 1e6
    planned = plan[-1].cum_mb
    trims = [r.option.name for r in plan if r.target is not None]
    notes.append(f"capture targets saved {whole - planned:.1f} MB of the plan: {whole:.1f} MB with every option's intervals whole, "
                 f"{planned:.1f} MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it "
                 f"is: an interval the plan's other options fetch anyway is no saving, so the options' savings ({', '.join(trims)}) need not "
                 "add up to the plan's")
    return notes


def _capture_note(rows, stat) -> str:
    """What expected_capture is: over which scans (held-out or not, as the statistics files say), and in
    which order it was reached."""
    if not rows:
        return ""
    by = defaultdict(list)
    for r in rows:
        by[r.held_out].append(f"{r.option.name} ({r.stats_scans} scans)")
    parts = []
    if by.get(True):
        parts.append("held-out for " + ", ".join(by[True]) + ": scans not used to learn the sinks")
    if by.get(False):
        parts.append("in-sample for " + ", ".join(by[False]) + ": the scans the sinks were learned from, so check the kept intervals "
                     "on other scans (`ngsdose sinks HELD_OUT --evaluate PREFIX.sinks.bed`)")
    if by.get(None):
        parts.append("of unrecorded provenance for " + ", ".join(by[None]) + ": their statistics files do not say whether the scans "
                     "were held out ('# held-out:', written by `ngsdose sinks --evaluate --held-out --stats`)")
    how = f"the {'10th percentile' if stat == 'p10' else 'median'} over the scans of its statistics of the capture of the intervals kept"
    byte = [r.option.name for r in rows if r.order == "byte" or r.shared]
    order = (f"; for {', '.join(byte)}, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, "
             "it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval "
             "held in any one scan, and the capture of the longest run of the statistics' own order kept whole" if byte else "")
    return f"expected capture is {how}{order}; " + "; ".join(parts)


def _fmt(x, spec):
    return "NA" if x is None else format(x, spec)


def table_lines(plan: Plan) -> list[str]:
    out = ["\t".join(PLAN_COLUMNS)]
    for r in plan.rows:
        o = r.option
        bp = sum(e - s for _, s, e in r.intervals) if o.kind != "unmapped" else None
        lost = r.full - r.expected if r.full is not None and r.expected is not None else None
        saved = r.mb_saved
        out.append("\t".join((o.name, o.status, o.tier, o.kind, str(len(r.intervals)) if o.kind != "unmapped" else "*", str(r.dropped),
                              _fmt(bp, "d"), _fmt(r.target, "g"), _fmt(r.expected, ".5f"), _fmt(r.full, ".5f"), _fmt(lost, ".5f"),
                              _fmt(r.mb, ".1f"), _fmt(saved, ".1f"), _fmt(r.cum_mb, ".1f"), _fmt(r.cum_mb_floor, ".1f"), _fmt(r.cum_pct, ".2f"),
                              {"read": "per read", "byte": "per byte", "read, fewer bytes": "per read (fewer bytes than per byte)"}.get(r.order, "NA"),
                              r.note or ".")))
    for o in plan.scan_only:
        where = " (panel in scan_panels.txt)" if o.panel in plan.scan_panels else ""
        out.append("\t".join((o.name, o.status, o.tier, o.kind, "0", "0", "0", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                              f"scan only: no sinks learned yet{where}")))
    return out


def write(plan: Plan, prefix: str, menu: Menu, panel_root=None, header: str = ""):
    root = Path(panel_root) if panel_root else None

    def where(rel):
        return str(root / rel) if root and not Path(rel).is_absolute() else str(menu.resolve(rel))

    with open(f"{prefix}.sinks.bed", "w") as fh:
        fh.write(f"# ngsdose fetchplan: {header}\n")
        for c, s, e, n in plan.bed():
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    Path(f"{prefix}.panels.txt").write_text("".join(where(p) + "\n" for p in plan.panels))
    Path(f"{prefix}.scan_panels.txt").write_text("".join(where(p) + "\n" for p in plan.scan_panels))
    Path(f"{prefix}.count_flags.txt").write_text("".join(f + "\n" for f in plan.flags))
    # under the menu's directory: a path relative to it, mapped by --panel-root like the panels'. Compared as written first (the
    # menu's controls are named relative to it, and a symlinked bundle resolves elsewhere), then by their real paths
    fa, base = Path(plan.controls_fasta), menu.path.parent
    if not fa.is_relative_to(base):
        fa, base = fa.resolve(), base.resolve()
    rel = str(fa.relative_to(base)) if fa.is_relative_to(base) else str(fa)
    Path(f"{prefix}.controls.txt").write_text(where(rel) + "\n")
    with open(f"{prefix}.plan.tsv", "w") as fh:
        fh.write(f"# ngsdose fetchplan: {header}\n")
        if plan.n_index:
            fh.write(f"# MB = 1e6 bytes of CRAM slices (and their containers' compression headers), median over {plan.n_index} index(es), "
                     "as `ngs-dose count -m fetch` reads them: one indexed fetch per run of touching or overlapping intervals of the plan, "
                     "each decoding every slice that overlaps it, so a slice under several runs is decoded once per run; cum_*: this option "
                     "and all above it; cum_mb_floor: the same with every slice decoded once, the floor a reader that sorted the plan's "
                     "slices would reach (the engine does not); mb_saved: what the option's capture target saved of the plan, the rest of "
                     "the plan as it is (an interval the plan's other options fetch anyway is no saving), capture_lost what it cost of its "
                     "expected capture\n")
        if plan.capture_note:
            fh.write(f"# {plan.capture_note}\n")
        if plan.classes_flag:
            fh.write(f"# count flags for the classes the panels define but the plan does not select: {plan.classes_flag}\n")
        fh.write(f"# controls FASTA for the fetch (-c): {Path(plan.controls_fasta).name}\n")
        for n in plan.notes:
            if n != plan.capture_note:
                fh.write(f"# {n}\n")
        fh.write("\n".join(table_lines(plan)) + "\n")


def list_menu(menu: Menu) -> list[str]:
    out = ["class\tstatus\ttier\tgroup\tkind\tpresets\tmeasures"]
    for o in sorted(menu.options.values(), key=lambda o: (TIERS.index(o.tier), o.order)):
        kind = f"subset of {o.parent}" if o.kind == "subset" else o.kind
        out.append("\t".join((o.name, o.status, o.tier, o.group, kind, ",".join(o.presets) or ".", o.measures)))
    out.append("")
    out.append("preset\tmembers\tdescription")
    for p, d in menu.presets.items():
        out.append(f"{p}\t{' '.join(menu.members(p))}\t{d}")
    return out
