"""Command line for the modelling layer: `ngsdose estimate | cohort | adjust | pcsweep | trios | sinks | selftest`."""
from __future__ import annotations

import argparse
import contextlib
import csv
import json
import re
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from . import __version__, cohort, contract, estimate, io, pcselect, resources, sinks, trios
from .tables import dump as _dump, load_result as _load_result, num as _num, summary_row, write_table


_EST: dict = {}


def _estimate_init(a):
    res = resources.Bundle(a.resources)
    panel, units = io.load_panel(res.panel), res.units()
    res.check_units(panel, units)
    _EST.update(args=a, res=res, panel=panel, units=units, feats=res.features(),
                anchors={} if a.gc_rule_anchors else res.anchors(), tables={}, lengths=res.contig_lengths(),
                regions=None if a.no_control_qc else res.regions())


def _output_name(path) -> str:
    """The estimate file of one input: its name with .json.gz (or .json) replaced by .estimate.json.gz."""
    n = Path(path).name
    for suf in (".json.gz", ".json", ".gz"):
        if n.endswith(suf) and len(n) > len(suf):
            n = n[: -len(suf)]
            break
    return n + ".estimate.json.gz"


def _estimate_one(job):
    """One counts file -> {path, row, resources, classes, warn}, or {path, error}: a file that cannot be
    read or estimated is reported, not raised, so that the other files still get their rows."""
    path, out = job
    a, res = _EST["args"], _EST["res"]
    try:
        counts = contract.load(path)
        found = contract.issues(counts, res)
        if a.no_control_qc:                                 # the bundle's control regions are not used then
            found = [("warn" if lv == "fatal" and "controls file" in m else lv, m) for lv, m in found]
        fatal = [m for lv, m in found if lv == "fatal"]
        if fatal:
            return dict(path=path, error="; ".join(fatal))
        L = estimate.nearest_table(counts, a.L)["l"]
        if not a.no_control_qc and L not in _EST["tables"]:
            _EST["tables"][L] = res.region_tables(L)
        r = estimate.estimate_sample(counts, _EST["panel"], _EST["units"], _EST["feats"], region_tables=_EST["tables"].get(L), L=a.L,
                                     window=a.window, min_kmers=a.min_kmers, anchors=_EST["anchors"],
                                     contig_lengths=_EST["lengths"], regions=_EST["regions"])
        r["resources"] = {k: counts.get(k) for k in ("panel_sha256", "controls_sha256", "sinks_sha256")}
        if out:
            _dump(r, out)
        classes = {c["name"]: (c.get("kind"), c.get("length"), c.get("panel_kmers")) for c in counts["classes"]}
        return dict(path=path, row=summary_row(r), resources=r["resources"], classes=classes, warn=[m for lv, m in found if lv == "warn"])
    except Exception as e:  # noqa: BLE001 - whatever it is, it is this file's failure, not the run's
        msg = str(e)
        return dict(path=path, error=msg[len(str(path)):].lstrip(": ") if msg.startswith(str(path)) else f"{type(e).__name__}: {msg}")


def _mixed_resources(done: list[dict]) -> list[str]:
    """What says that the counts files were not made with one set of resources."""
    out = []
    sets: dict[frozenset, str] = {}
    for o in done:
        h = o["resources"].get("panel_sha256")
        s = frozenset([h] if isinstance(h, str) else (h or ()))
        if s:
            sets.setdefault(s, o["row"]["sample"])
    # a scan made with extra panels holds the bundle panel's hash among its own: only sets of which
    # neither holds the other are different panels
    clash = next(((sets[p], sets[q]) for p in sets for q in sets if not (p <= q or q <= p)), None)
    if clash:
        out.append(f"the counts files were made with different panel files ({clash[0]} and {clash[1]}: neither's panel_sha256 holds the "
                   "other's). A different panel changes what is counted: do not analyse such files as one cohort.")
    fields = ("kind", "length", "panel_kmers")
    seen: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for o in done:
        for name, sig in o["classes"].items():
            for f, v in zip(fields, sig):
                if v is not None:
                    seen[name][f].add(v)
    for name, per in seen.items():
        diff = [f"{f} {' / '.join(map(str, sorted(v, key=str)))}" for f, v in per.items() if len(v) > 1]
        if diff:
            out.append(f"{name} was not counted with the same k-mers in every file ({'; '.join(diff)}): its values are not comparable across them.")
    for k, what in (("controls_sha256", "controls"), ("sinks_sha256", "sinks")):
        if len({o["resources"].get(k) for o in done if o["resources"].get(k)}) > 1:
            out.append(f"the counts files were not all made with the same {what} file. " + (
                "Different controls files are accepted only if their control regions are identical (extra known-truth regions are then "
                "simply missing for some samples)." if what == "controls" else
                "Different sinks change what a fetch reads: do not analyse such files as one cohort."))
    return out


def cmd_estimate(a):
    names = [_output_name(p) for p in a.counts] if a.outdir else [None] * len(a.counts)
    if a.outdir:
        dup = sorted(n for n, k in Counter(names).items() if k > 1)
        if dup:
            raise SystemExit(f"ngsdose estimate: inputs with the same file name would write the same estimate file ({', '.join(dup[:5])}): "
                             "each estimate is named after its input; estimate them into separate directories (scan and fetch counts "
                             "of the same genomes, for one), or give the inputs distinct names")
        Path(a.outdir).mkdir(parents=True, exist_ok=True)
    jobs = [(p, Path(a.outdir) / n if n else None) for p, n in zip(a.counts, names)]
    try:
        _estimate_init(a)                                   # the bundle is checked once, before any worker starts
    except (OSError, ValueError) as e:
        raise SystemExit(f"ngsdose estimate: {e}") from None
    out = []
    with contextlib.ExitStack() as stack:
        if a.jobs > 1 and len(jobs) > 1:
            import multiprocessing as mp
            pool = stack.enter_context(mp.get_context("spawn").Pool(a.jobs, initializer=_estimate_init, initargs=(a,)))
            it = pool.imap(_estimate_one, jobs, chunksize=4)
        else:
            it = map(_estimate_one, jobs)
        for o in it:                                        # in input order, each file's messages as it finishes
            for m in o.get("warn", ()):
                print(f"[estimate] WARNING: {m}", file=sys.stderr)
            if "error" in o:
                print(f"[estimate] FAILED {o['path']}: {o['error']}", file=sys.stderr)
            out.append(o)
    done = [o for o in out if "error" not in o]
    for m in _mixed_resources(done):
        print(f"[estimate] WARNING: {m}", file=sys.stderr)
    by = defaultdict(list)
    for o in done:
        by[o["row"]["sample"]].append(o)
    for s, os_ in by.items():
        if len(os_) > 1:
            which = ", ".join("%s (%s)" % (o["path"], o["row"]["mode"]) for o in os_)
            print(f"[estimate] WARNING: sample {s} is in {len(os_)} inputs ({which}); "
                  f"{'each has its own estimate file, but ' if a.outdir else ''}a cohort takes one estimate per sample", file=sys.stderr)
    if "unknown" in by:
        print("[estimate] WARNING: sample 'unknown' - the counts name no sample (no @RG SM in the input and no -s when counting)", file=sys.stderr)
    na = Counter((k[: -len(".status")], v) for o in done for k, v in o["row"].items() if k.endswith(".status") and v not in (None, "ok"))
    if na:
        print("[estimate] not estimated (NA in the table): " + "; ".join(f"{c} in {n} sample(s): {v}" for (c, v), n in sorted(na.items())),
              file=sys.stderr)
    write_table([o["row"] for o in done], a.table)
    failed = len(out) - len(done)
    if failed:
        raise SystemExit(f"ngsdose estimate: {failed} of {len(out)} counts files failed (listed above); the table holds the other {len(done)}")


def _estimates(paths):
    for p in paths:
        try:
            yield _load_result(p)
        except (OSError, EOFError, ValueError) as e:
            raise SystemExit(f"ngsdose cohort: {p}: unreadable estimate file ({type(e).__name__}: {e})") from None


def cmd_cohort(a):
    # one estimate file at a time: only the window vectors, the control residuals and the summary
    # row of each sample are kept, so a cohort of thousands fits in a few hundred MB
    try:
        anchors = {} if a.gc_rule_anchors else resources.Bundle(a.resources).anchors()
        eff_in = json.loads(Path(a.efficiencies).read_text()) if a.efficiencies else None
    except (OSError, ValueError) as e:
        raise SystemExit(f"ngsdose cohort: {e}") from None
    try:
        rows, eff, _ = cohort.cohort_table(_estimates(a.estimates), anchors, max_window_sd=a.max_window_sd,
                                           n_profile_pcs=a.profile_pcs, n_control_pcs=a.control_pcs, mp_margin=a.mp_margin,
                                           efficiencies=eff_in, log=lambda m: print(m, file=sys.stderr))
    except ValueError as e:
        raise SystemExit(f"ngsdose cohort: {e}") from None
    if a.save_efficiencies:
        Path(a.save_efficiencies).write_text(json.dumps(eff))
    write_table(rows, a.table)


def _read_table(path) -> tuple[list[str], list[dict]]:
    """Header and rows of a per-sample TSV. A sample given twice is refused: which of its rows a
    per-sample analysis would use is arbitrary (a scan and a fetch of one genome, two files named alike)."""
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = list(rd)
        header = list(rd.fieldnames or [])
    if "sample" not in header:
        raise SystemExit(f"{path}: no 'sample' column")
    dup = sorted(s for s, n in Counter(r["sample"] for r in rows).items() if n > 1)
    if dup:
        raise SystemExit(f"{path}: sample {', '.join(dup[:5])}{' and others' if len(dup) > 5 else ''} given more than once; "
                         "give every sample one row (keep scan and fetch estimates of the same genomes in separate tables)")
    return header, rows


def _need_columns(header, columns, path, cmd):
    missing = [c for c in dict.fromkeys(columns) if c not in header]
    if missing:
        raise SystemExit(f"ngsdose {cmd}: no column {', '.join(missing)} in {path}")


def _values(rows, column, path):
    vals = {}
    for r in rows:
        try:
            vals[r["sample"]] = float(r[column])
        except (ValueError, KeyError, TypeError):
            pass
    if not vals:
        raise ValueError(f"no numeric values in column {column!r} of {path}")
    return vals


@contextlib.contextmanager
def _pedigree_warnings(cmd):
    """The pedigree's and the trio analysis's warnings, each printed once, as the command's own."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", trios.PedigreeWarning)
        try:
            yield
        finally:
            said = set()
            for w in caught:
                if issubclass(w.category, trios.PedigreeWarning):
                    if str(w.message) not in said:
                        said.add(str(w.message))
                        print(f"[{cmd}] WARNING: {w.message}", file=sys.stderr)
                else:
                    warnings.showwarning(w.message, w.category, w.filename, w.lineno)


def cmd_trios(a):
    header, rows = _read_table(a.table)
    _need_columns(header, list(a.columns) + ([a.compare_to] if a.compare_to else []), a.table, "trios")
    out, failed = {}, {}
    with _pedigree_warnings("trios"):
        ped, pop = trios.load_pedigree(a.pedigree, population=a.population)
        if not ped:
            raise SystemExit(f"ngsdose trios: {a.pedigree} holds no complete trio (a child with both parents given); check its layout")
        centre = None if a.no_population_centring else pop
        for col in a.columns:
            try:
                out[col] = trios.transmission(_values(rows, col, a.table), ped, centre, n_perm=a.perm)
            except (ValueError, KeyError) as err:           # too few complete trios, or no numbers: the other columns still get their row
                failed[col] = str(err)
        ci = lambda t, k: ("(%.3f,%.3f)" % t[k + "_ci95"]) if k + "_ci95" in t else "NA"
        hdr = ("column", "n_trios", "R_midparent", "CI95", "slope", "se", "R_single", "R_mendel", "spousal_r", "CI95", "error_cv", "perm_null")
        print("\t".join(hdr))
        for col in a.columns:
            if col in failed:
                print("\t".join([col, "0"] + ["NA"] * (len(hdr) - 2)) + f"\t# {failed[col]}")
                continue
            t = out[col]
            print("\t".join([col, str(t["n_trios"]), f"{t['reliability_midparent']:.3f}", ci(t, "reliability_midparent"),
                             f"{t['midparent_slope']:.3f}", f"{t['midparent_slope_se']:.3f}", f"{t['reliability_single_parent']:.3f}",
                             f"{t['reliability_mendel']:.3f}", f"{t['spousal_r']:+.3f}", ci(t, "spousal_r"), f"{t['error_cv']:.4f}",
                             f"{t.get('perm_null_mean', float('nan')):+.3f}"]))
        if a.compare_to:
            try:
                base = _values(rows, a.compare_to, a.table)
            except ValueError as err:
                raise SystemExit(f"ngsdose trios: --compare-to: {err}") from None
            print("\n# paired family bootstrap: R(column) - R(%s)" % a.compare_to)
            print("column\tdelta_R\tCI95\tP(column better)")
            for col in a.columns:
                if col == a.compare_to or col in failed:
                    continue
                try:
                    c = trios.compare(_values(rows, col, a.table), base, ped, centre)
                except (ValueError, KeyError) as err:
                    print(f"{col}\tNA\tNA\tNA\t# {err}")
                    continue
                out[col]["vs_" + a.compare_to] = c
                print(f"{col}\t{c['delta']:+.4f}\t({c['ci95'][0]:+.4f},{c['ci95'][1]:+.4f})\t{c['p_a_better']:.3f}")
    if a.json:
        _dump({**out, **{col: dict(error=msg) for col, msg in failed.items()}}, a.json)
    if failed and len(failed) == len(set(a.columns)):
        raise SystemExit(f"ngsdose trios: none of the {len(failed)} columns could be analysed (reasons above)")


def _float(x) -> float:
    try:
        return float(x)
    except ValueError:
        return float("nan")


def load_pcs(path, n_pc=None, strip=(".by1000.", ".")):
    """NGS-PCA svd.pcs.txt -> {sample: PCs}. NGS-PCA names samples after the mosdepth file, which
    leaves a suffix such as `.by1000.`; it is removed so that ids match the alignment's SM tag.
    A value that is not a number is NaN."""
    pcs = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            name = p[0]
            for suf in strip:
                if name.endswith(suf):
                    name = name[: -len(suf)]
                    break
            pcs[name] = [_float(x) for x in (p[1:] if n_pc is None else p[1:1 + n_pc])]
    return pcs


def _clamp(asked: int, available: int, what: str):
    """A number of PCs asked for outright is used as far as it can be: a pipeline should not die at its
    last step because one of two PC sets is smaller than the other, but it has to say so."""
    if asked > available:
        print(f"WARNING: --n-pc {asked}, but there are {available} {what}: using {available}", file=sys.stderr)
        return available, f"--n-pc {asked}, clamped to the {available} available"
    return asked, f"--n-pc {asked}"


def _pcs_and_choice(a, rows):
    """{sample: every available PC (NaN where the table has NA)}, the number to use, and how it was
    chosen. `--n-pc mp` (the default) counts the components above the Marchenko-Pastur edge of the
    noise bulk: for NGS-PCA's PCs from the singular values and bin list it writes beside svd.pcs.txt,
    for the table's own control-region PCs from the count `ngsdose cohort` recorded (ctrlPC_mp)."""
    auto = str(a.n_pc).lower() == "mp"
    if a.pcs:
        pcs = load_pcs(a.pcs)
        if not pcs:
            raise SystemExit(f"{a.pcs} holds no sample")
        n_avail = min(len(v) for v in pcs.values())
        if not auto:
            return pcs, *_clamp(int(a.n_pc), n_avail, "PCs in " + str(a.pcs))
        d = Path(a.pcs).parent
        sv_path = Path(a.singular_values) if a.singular_values else d / "svd.singularvalues.txt"
        n_feat = a.n_features
        if n_feat is None and (d / "svd.bins.txt").exists():
            with open(d / "svd.bins.txt") as fh:
                n_feat = sum(1 for line in fh if line.strip())
        if not sv_path.exists() or not n_feat:
            raise SystemExit(f"--n-pc mp needs the singular values ({sv_path}) and the number of bins (svd.bins.txt beside the PCs, or "
                             "--n-features) of the SVD the PCs came from; or give a number: --n-pc 20")
        sv = []
        with open(sv_path) as fh:
            for line in fh:
                try:
                    sv.append(float(line.split()[-1]))
                except (ValueError, IndexError):
                    continue
        sel = pcselect.mp_select(sv, len(pcs), n_feat, margin=a.mp_margin)
        return pcs, min(sel.n_pc, n_avail), sel.describe()
    cols = [c for c in rows[0] if re.fullmatch(r"ctrlPC\d+", c)] if rows else []
    cols.sort(key=lambda c: int(c[6:]))
    if not cols:
        raise SystemExit("no --pcs given and the table has no ctrlPC columns; run `ngsdose cohort` with --control-pcs")
    pcs = {r["sample"]: [_num(r, c) for c in cols] for r in rows}
    if not auto:
        return pcs, *_clamp(int(a.n_pc), len(cols), "control-region PCs the table carries (`ngsdose cohort --control-pcs N` writes more)")
    mp = next((r["ctrlPC_mp"] for r in rows if r.get("ctrlPC_mp", "NA") not in ("NA", "")), None)
    if mp is None:
        raise SystemExit("--n-pc mp: the table does not record a Marchenko-Pastur count (ctrlPC_mp); re-run `ngsdose cohort`, or give --n-pc N")
    return pcs, min(int(float(mp)), len(cols)), "the Marchenko-Pastur count recorded by `ngsdose cohort` (ctrlPC_mp)"


def _with_pcs(a, rows, pcs, k, need, cmd):
    """The rows whose sample has its first k PCs; the others are named. Fewer than `need` such rows
    is an error."""
    keep = [r for r in rows if r["sample"] in pcs and all(np.isfinite(pcs[r["sample"]][:k]))]
    miss = [r["sample"] for r in rows if not (r["sample"] in pcs and all(np.isfinite(pcs[r["sample"]][:k])))]
    what = f"PCs in {a.pcs}" if a.pcs else f"finite ctrlPC1..ctrlPC{k} (`ngsdose cohort` writes NA for a sample without control residuals)"
    hint = ""
    if a.pcs and len(miss) > len(rows) / 2:
        hint = (f". Most ids do not match: the table has, for example, {rows[0]['sample']!r} and {a.pcs} {next(iter(pcs))!r} "
                "(suffixes .by1000. and . are removed from the PC file's names)")
    if miss:
        print(f"[{cmd}] {len(miss)} of {len(rows)} samples have no {what}: {', '.join(miss[:10])}{', ...' if len(miss) > 10 else ''}"
              f"{'; their .adj values are NA' if cmd == 'adjust' else '; they are left out'}{hint}", file=sys.stderr)
    if len(keep) < need:
        raise SystemExit(f"only {len(keep)} samples have coverage PCs; need more than n_pc + 10 = {need}{hint}")
    return keep


def cmd_adjust(a):
    header, rows = _read_table(a.table)
    _need_columns(header, a.columns, a.table, "adjust")
    pcs, k, how = _pcs_and_choice(a, rows)
    keep = _with_pcs(a, rows, pcs, k, k + 10, "adjust")
    log = not a.no_log
    print(f"[adjust] regressing out {k} PCs ({how})", file=sys.stderr)
    if k > len(keep) / 10:
        print(f"[adjust] WARNING: {k} PCs for {len(keep)} samples is more than one regressor per ten samples - a count chosen on a larger "
              "cohort's SVD does not carry over to a subset; consider --n-pc N (and see `ngsdose pcsweep`)", file=sys.stderr)
    X = np.array([pcs[r["sample"]][:k] for r in keep]).reshape(len(keep), k)
    for col in a.columns:
        for r in rows:                                      # samples without PCs keep their row, with NA
            r[f"{col}.adj"] = "NA"
        y = np.array([_num(r, col) for r in keep])
        if k:
            n_ok = int(cohort.adjust_mask(y, X, log).sum())
            if n_ok < k + 10:
                print(f"[adjust] WARNING: {col}: {n_ok} usable values ({'positive' if log else 'finite'}) among the {len(keep)} samples with PCs; "
                      f"{k} PCs need at least {k + 10}: {col}.adj is NA", file=sys.stderr)
                continue
            if k > n_ok / 10 and n_ok < len(keep):
                print(f"[adjust] WARNING: {col}: {k} PCs for its {n_ok} values is more than one regressor per ten", file=sys.stderr)
            adj, r2 = cohort.adjust_for_covariates(y, X, log=log)
        else:
            adj, r2 = y.copy(), 0.0
        for r, v in zip(keep, adj):
            r[f"{col}.adj"] = "NA" if not np.isfinite(v) else round(float(v), 4)
        n = int(np.isfinite(adj).sum())
        # k regressors explain k/(n-1) of pure noise: report what is left after that
        r2_adj = 1 - (1 - r2) * (n - 1) / max(n - k - 1, 1)
        print(f"[adjust] {col}: {k} PCs explain {100 * r2:.1f}% of the variance of {'log ' if log else ''}{col} "
              f"(n={n}; {100 * k / max(n - 1, 1):.1f}% expected by chance; adjusted R2 {100 * r2_adj:.1f}%)", file=sys.stderr)
    write_table(rows, a.out)


def _sex_from(rows):
    """'M' / 'F' per sample from the known-truth columns themselves: chrY if it was measured, else chrX."""
    out = {}
    for r in rows:
        y, x = r.get("truth.chrY", "NA"), r.get("truth.chrX", "NA")
        if y not in ("NA", ""):
            out[r["sample"]] = "M" if float(y) > 0.5 else "F"
        elif x not in ("NA", ""):
            out[r["sample"]] = "M" if float(x) < 1.5 else "F"
    return out


def _expected_copies(a) -> dict[str, float]:
    """Copies per genome that a class has in everyone (bundle.json `expected_copies`: DJ = 10 in GRCh38)."""
    try:
        return dict(resources.Bundle(a.resources).meta.get("expected_copies", {}))
    except FileNotFoundError as e:
        if a.resources:
            raise SystemExit(f"ngsdose pcsweep: {e}") from None
        print(f"[pcsweep] no bundle ({e}): no class is scored against expected copies", file=sys.stderr)
        return {}


def cmd_pcsweep(a):
    header, rows = _read_table(a.table)
    extra = {}
    for spec in a.truth or []:
        col, sep, val = spec.rpartition("=")
        try:
            v = float(val)
        except ValueError:
            v = float("nan")
        if not sep or not col or not (np.isfinite(v) and v > 0):
            raise SystemExit(f"ngsdose pcsweep: --truth {spec!r}: expected COLUMN=VALUE with a positive number")
        extra[col] = v
    _need_columns(header, list(a.columns) + list(extra), a.table, "pcsweep")
    if a.population and not a.pedigree:
        raise SystemExit("ngsdose pcsweep: --population needs --pedigree")
    pcs, k_mp, how = _pcs_and_choice(a, rows)
    n_avail = min(len(v) for v in pcs.values())
    want = min(n_avail, a.max_pc if a.max_pc else max(2 * k_mp, 20))
    keep = _with_pcs(a, rows, pcs, want, k_mp + 10, "pcsweep")
    max_pc = min(want, max(len(keep) // 5, 1))
    P = np.array([pcs[r["sample"]][:max_pc] for r in keep]).reshape(len(keep), max_pc)
    samples = [r["sample"] for r in keep]
    num = lambda col: np.array([_num(r, col) for r in keep])
    sex = _sex_from(keep)
    male = np.array([sex.get(s) == "M" for s in samples]); female = np.array([sex.get(s) == "F" for s in samples])
    # the known truths the bundle builds in, where the table has them; --truth COLUMN=VALUE adds others
    truths = {}
    if "truth.auto" in header:
        truths["truth.auto"] = np.full(len(keep), 2.0)
    if "truth.chrX" in header:
        truths["truth.chrX"] = np.where(male, 1.0, np.where(female, 2.0, np.nan))
    if "truth.chrY" in header:
        truths["truth.chrY"] = np.where(male, 1.0, np.nan)               # a female's truth is zero, which has no log
    for cls, n in sorted(_expected_copies(a).items()):
        for cand in (f"{cls}.cn", f"{cls}.cn_single"):
            if cand in header:
                truths[cand] = np.full(len(keep), float(n))
                break
    for col, v in extra.items():
        truths[col] = np.full(len(keep), v)
    table = {c: num(c) for c in list(truths) + list(a.columns)}
    trio_list = population = None
    with _pedigree_warnings("pcsweep"):
        if a.pedigree:
            trio_list, population = trios.load_pedigree(a.pedigree, population=a.population)
            if a.no_population_centring:
                population = None
        res = pcselect.sweep(table, P, max_pc, truths, list(a.columns), samples, trio_list, population, folds=a.folds, n_boot=a.boot)
    write_table([{k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()} for r in res], a.out)
    print(f"[pcsweep] {len(keep)} samples, PCs 0..{max_pc}, {a.folds}-fold cross-validation; default choice: {k_mp} PCs ({how})", file=sys.stderr)
    rec = pcselect.recommend(res)
    print("column\tkind\tmeasure\tat_0_PCs\tat_default\tbest_n_pc\tat_best\tfewest_within_1se\tthere", file=sys.stderr)
    for col, r in rec.items():
        at_def = next((x for x in res if x["column"] == col and x["n_pc"] == min(k_mp, max_pc)), {})
        key = "sd_log_robust" if r["kind"] == "truth" else "R_midparent"
        print(f"{col}\t{r['kind']}\t{'robust SD of log(estimate/truth), cross-validated' if r['kind'] == 'truth' else 'transmission reliability'}"
              f"\t{r['at_0']:.4f}\t{at_def.get(key, float('nan')):.4f}\t{r['best']}\t{r['at_best']:.4f}\t{r['pick']}\t{r['at_pick']:.4f}", file=sys.stderr)


def cmd_sinks(a):
    if a.evaluate:
        try:
            bed = sinks.read_bed(a.evaluate)
        except (OSError, ValueError) as e:
            sys.exit(f"ngsdose sinks: {e}")
        print("sample\tclass\tscan_reads\tcaptured\tfraction\tunmapped")
        bad = []
        for path in a.counts:                              # one file at a time: a cohort of scans fits
            try:
                c = io.load_counts(path)
                sinks.check_scan(c, a.allow_cut)
            except ValueError as e:
                print(f"[sinks] not evaluated: {path}: {e}", file=sys.stderr)
                bad.append(path)
                continue
            if sinks.looks_cut(c):
                print(f"# {c['sample']}: {100 * sinks.cut_share(c):.0f}% of its primary reads lie in the control regions - a cut along a fetch plan, "
                      "not a whole-file scan; its capture says nothing", file=sys.stderr)
            unmapped = Counter()
            for p in c["placements"]:
                if p["contig"] == "*":
                    unmapped[p["class"]] += p["reads"]
            for cls, (tot, inside) in sinks.capture(c, bed).items():
                print(f"{c['sample']}\t{cls}\t{tot}\t{inside}\t{inside / max(tot, 1):.5f}\t{unmapped[cls]}")
        if bad:
            sys.exit(f"ngsdose sinks: {len(bad)} of {len(a.counts)} counts files were not evaluated (listed above)")
        return
    try:
        rows, stats = sinks.learn(a.counts, a.min_frac, a.pad, a.min_reads, classes=a.classes, allow_cut=a.allow_cut,
                                  allow_mixed_panels=a.allow_mixed_panels, log=lambda m: print(m, file=sys.stderr))
    except ValueError as e:
        sys.exit(f"ngsdose sinks: {e}")
    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    try:
        for contig, s0, e0, cls in rows:
            fh.write(f"{contig}\t{s0}\t{e0}\t{cls}\n")
    finally:
        if fh is not sys.stdout:
            fh.close()
    for cls, per in sorted(stats.items()):
        v = sorted(per.values())
        n = sum(1 for r in rows if r[3] == cls)
        print(f"[sinks] {cls}: {n} intervals, {sum(r[2] - r[1] for r in rows if r[3] == cls):,} bp; "
              f"capture in these scans min {v[0]:.5f} median {v[len(v) // 2]:.5f} over {len(v)} scan(s)"
              f"{'' if n else ' - WARNING: no bin reached the rule, so the BED has no sinks for it'}", file=sys.stderr)


def cmd_selftest(a):
    from . import selftest
    sys.exit(0 if selftest.run(verbose=True) else 1)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ngsdose", description=__doc__)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("estimate", help="counts JSON -> per-sample estimates")
    e.add_argument("counts", nargs="+")
    e.add_argument("-r", "--resources", default=None, help="resource bundle directory (default: packaged GRCh38)")
    e.add_argument("-o", "--outdir", help="write one estimate per input, named after the input file: X.json.gz -> X.estimate.json.gz")
    e.add_argument("-t", "--table", default="-", help="summary TSV (default stdout); a file that fails is listed on stderr, "
                   "the table holds the others, and the exit status is 1")
    e.add_argument("-L", type=int, default=None, help="fragment-GC window (default: nearest tabulated to the insert median)")
    e.add_argument("--window", type=int, default=250)
    e.add_argument("--min-kmers", type=int, default=20)
    e.add_argument("-j", "--jobs", type=int, default=1, help="parallel processes")
    e.add_argument("--no-control-qc", action="store_true")
    e.add_argument("--gc-rule-anchors", action="store_true", help="ignore the bundle's anchor intervals; anchor on fragment GC 40-60%%")
    e.set_defaults(fn=cmd_estimate)

    c = sub.add_parser("cohort", help="calibrate window efficiencies across samples")
    c.add_argument("estimates", nargs="+")
    c.add_argument("-t", "--table", default="-")
    c.add_argument("--efficiencies", help="apply a saved efficiency table instead of fitting one")
    c.add_argument("--save-efficiencies")
    c.add_argument("-r", "--resources", default=None)
    c.add_argument("--gc-rule-anchors", action="store_true")
    c.add_argument("--max-window-sd", type=float, default=None)
    c.add_argument("--profile-pcs", type=int, default=3)
    c.add_argument("--mp-margin", type=float, default=0.01,
                   help="how far (relative) a component has to clear the fitted noise edge to count towards ctrlPC_mp (default 0.01)")
    c.add_argument("--control-pcs", default="mp",
                   help="components of the control regions' residual depth to write as ctrlPC columns: a number, or 'mp' (default): "
                        "count the components above the Marchenko-Pastur edge of the noise bulk, record the count (ctrlPC_mp) "
                        "and write twice as many (at least 20, but never more than one per five samples; nothing, not even "
                        "ctrlPC_mp, below 10 samples with control residuals) so that `pcsweep` can look beyond it")
    c.set_defaults(fn=cmd_cohort)

    t = sub.add_parser("trios", help="transmission reliability from a pedigree")
    t.add_argument("table")
    t.add_argument("-p", "--pedigree", required=True, help="the 1000 Genomes layout (FamilyID SampleID FatherID MotherID Sex Population ...), a PLINK "
                   "PED/FAM, or a child father mother [population] table; whitespace-separated (or tab-separated under a header). "
                   "A population is read only from a column a header names (or from --population FILE); without a header "
                   "the file gives trios only")
    t.add_argument("-c", "--columns", nargs="+", required=True)
    t.add_argument("--perm", type=int, default=1000)
    t.add_argument("--compare-to", help="baseline column: paired bootstrap of the reliability difference of every other column against it")
    t.add_argument("--population", metavar="FILE", help="sample -> population (or ancestry cluster) to centre within: two columns, or a header "
                   "naming a sample and a population column; its labels replace the pedigree's")
    t.add_argument("--no-population-centring", action="store_true",
                   help="do not centre within populations (a cohort of one population)")
    t.add_argument("--json")
    t.set_defaults(fn=cmd_trios)

    def pc_source(sp):
        sp.add_argument("--pcs", help="NGS-PCA svd.pcs.txt (default: the ctrlPC columns written by `ngsdose cohort`)")
        sp.add_argument("--n-pc", default="mp",
                        help="how many PCs to regress out: a number, or 'mp' (default): the components above the Marchenko-Pastur "
                             "edge of the noise bulk")
        sp.add_argument("--mp-margin", type=float, default=0.01,
                        help="for --n-pc mp: how far (relative) a component has to clear the fitted noise edge to count (default 0.01)")
        sp.add_argument("--singular-values", help="singular values of the SVD behind --pcs (default: svd.singularvalues.txt beside it)")
        sp.add_argument("--n-features", type=int, help="number of bins in that SVD (default: the lines of svd.bins.txt beside --pcs)")

    j = sub.add_parser("adjust", help="residualise estimates on coverage PCs",
                       description="Every row of the table is written; a sample without PCs, or a column with fewer than n_pc + 10 "
                                   "usable values, gets NA in its .adj column.")
    j.add_argument("table")
    pc_source(j)
    j.add_argument("-c", "--columns", nargs="+", required=True)
    j.add_argument("--no-log", action="store_true")
    j.add_argument("-o", "--out", default="-")
    j.set_defaults(fn=cmd_adjust)

    w = sub.add_parser("pcsweep", help="what each further PC does: cross-validated error of the known truths, transmission of the classes",
                       description="Regress out 0, 1, 2, ... PCs and report, for every number: the cross-validated error of the "
                                   "known-truth columns (held-out autosomal = 2, chrX and chrY by sex, a class's expected copies from "
                                   "the bundle - DJ = 10; --truth adds others) and, with a pedigree, the transmission reliability of "
                                   "the class columns and its paired difference from no adjustment. The evidence on which the default "
                                   "number of PCs can be overruled.")
    w.add_argument("table")
    pc_source(w)
    w.add_argument("-c", "--columns", nargs="+", default=[], help="class columns (estimates without a known truth)")
    w.add_argument("--truth", nargs="+", metavar="COLUMN=VALUE", help="further columns with a known value")
    w.add_argument("-r", "--resources", default=None, help="bundle whose expected_copies give known truths (default: packaged GRCh38)")
    w.add_argument("-p", "--pedigree")
    w.add_argument("--population", metavar="FILE", help="as for `ngsdose trios`")
    w.add_argument("--no-population-centring", action="store_true")
    w.add_argument("--max-pc", type=int, default=0, help="sweep up to this many PCs (default: twice the default choice, at least 20; "
                   "always capped at the PCs available and at one per five samples)")
    w.add_argument("--folds", type=int, default=10)
    w.add_argument("--boot", type=int, default=300)
    w.add_argument("-o", "--out", default="-")
    w.set_defaults(fn=cmd_pcsweep)

    k = sub.add_parser("sinks", help="learn fetch-mode sink intervals from scan-mode counts")
    k.add_argument("counts", nargs="+")
    k.add_argument("-o", "--out", default="-")
    k.add_argument("--min-frac", type=float, default=1e-5,
                   help="keep the placement bins that hold class reads in any aligned 10-kb window which, in at least one scan, holds "
                        ">= this fraction of the class's reads (unmapped included) and >= --min-reads reads (default 1e-5)")
    k.add_argument("--evaluate", metavar="BED", help="instead of learning, report how much of each class an existing sinks BED captures "
                   "(a placement bin counts when all of it lies inside the class's intervals; unmapped reads count as missed)")
    k.add_argument("--min-reads", type=int, default=25, help="minimum class reads in such a 10-kb window (default 25)")
    k.add_argument("--pad", type=int, default=1000, help="merge a class's kept bins when the gap between them is at most this many bp, then "
                   "extend each merged stretch by this many bp on both sides, clipped to the contig (default 1000)")
    k.add_argument("--classes", nargs="+", metavar="CLASS",
                   help="also learn sinks for these compositional classes (TEL and the satellite families: in NYGC bwa-mem scans against "
                        "the 1000 Genomes analysis set, the aligner concentrates each family's reads in a stable set of intervals, 0.2-60 Mb, "
                        "that held >= 99.8%% of the class in held-out scans (99.85%% in two of three random draws); HSat1B is the exception at >= 96.6%%, part of it unmapped. "
                        "No satellite sinks ship with the bundle. Check capture on held-out scans with --evaluate, and re-learn for "
                        "every aligner and reference)")
    k.add_argument("--allow-cut", action="store_true", help="accept counts that look like a cut along a fetch plan")
    k.add_argument("--allow-mixed-panels", action="store_true", help="pool scans counted with different panel files")
    k.set_defaults(fn=cmd_sinks)

    s = sub.add_parser("selftest", help="simulation-based check of the estimator; needs no data")
    s.set_defaults(fn=cmd_selftest)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
