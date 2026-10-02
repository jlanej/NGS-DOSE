"""Command line for the modelling layer: `ngsdose estimate | cohort | adjust | pcsweep | trios | sinks | fetchplan | control-sets | selftest`."""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import json
import re
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from . import __version__, cohort, contract, estimate, fetchplan, io, karyotype, pcselect, resources, sinks, trios
from .tables import dump as _dump, load_result as _load_result, num as _num, summary_row, write_table


_EST: dict = {}


def _estimate_init(a):
    res = resources.Bundle(a.resources)
    panel, units = io.load_panel(res.panel), res.units()
    res.check_units(panel, units)
    sub_options, known_sinks = estimate.fetch_context(res.dir, res.sinks, getattr(a, "fetch_sinks", None) or ())
    ex_dir, looked = estimate.experimental_layout(res.dir)
    _EST.update(args=a, res=res, panel=panel, units=units, feats=res.features(),
                # a bare copy of the bundle's directory has no experimental resources: say so, and let estimate_sample
                # treat fetches through other sinks BEDs as unverified rather than measured
                experimental_missing=None if ex_dir is not None else str(looked),
                bundle_sinks_sha256=hashlib.sha256(Path(res.sinks).read_bytes()).hexdigest(),
                sinks_pipeline=(res.meta.get("sinks_learned_from") or {}).get("pipeline"),
                anchors={} if a.gc_rule_anchors else res.anchors(), tables={}, lengths=res.contig_lengths(),
                regions=None if a.no_control_qc else res.regions(),
                # a fetch made with a lighter controls file the bundle names, and units of positional
                # classes from experimental panels (resources/experimental/candidates/units, NGSDOSE_EXTRA_UNITS)
                subsets=res.control_subsets(), experimental=res.experimental(),
                # named subsets of classes' sinks (resources/experimental/subsets), and the sinks BEDs a fetch
                # may have been made with, by hash: the bundle's, the experimental ones, the subsets', --fetch-sinks
                sub_options=sub_options, known_sinks=known_sinks)


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
        found = contract.issues(counts, res, _EST["known_sinks"])
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
                                     contig_lengths=_EST["lengths"], regions=_EST["regions"],
                                     control_subsets=_EST["subsets"], experimental=_EST["experimental"],
                                     sub_options=_EST["sub_options"], known_sinks=_EST["known_sinks"],
                                     experimental_missing=_EST["experimental_missing"], bundle_sinks_sha256=_EST["bundle_sinks_sha256"],
                                     sinks_pipeline=_EST["sinks_pipeline"])
        r["resources"] = {k: counts.get(k) for k in ("panel_sha256", "controls_sha256", "sinks_sha256")}
        unverified = {k: v["parent"] for k, v in (r.get("sub_options") or {}).items() if v["status"] == "unverified"}
        parents = sorted(n for n, c in r["classes"].items() if c.get("status") == "unverified")
        if unverified or parents:
            # the fetch's intervals are not known: whether it read these classes' sinks whole, or only some of them, cannot be told
            found.append(("warn", f"{counts.get('sample')}: a fetch made with a sinks BED not known here ({counts.get('sinks')}, sha256 "
                                  f"{str(counts.get('sinks_sha256'))[:12]}...): {', '.join(sorted(set(parents) | set(unverified.values())))} "
                                  "may have been read whole or only at the intervals of the sub-options "
                                  f"{', '.join(unverified) or 'of the fetch menu'}, so these classes get no value and the sub-options no "
                                  "mass (status unverified); pass the BED with --fetch-sinks"))
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
    sets_used = Counter(o["row"].get("controls_subset") or "all" for o in done if o["row"].get("controls_used") is not None)
    if len(sets_used) > 1:
        out.append("the counts files were made with different sets of control regions (" + ", ".join(f"{n}: {k} file(s)" for n, k in
                   sorted(sets_used.items())) + "). A lighter set moves levels a little (for lite200 against all the controls: the 45S "
                   "rDNA by +0.2%, HSat1B by +1.5% with an SD of 2.3%): say so, or analyse them apart.")
    for k, what in (("controls_sha256", "controls"), ("sinks_sha256", "sinks")):
        if k == "controls_sha256" and len(sets_used) > 1:
            continue                                        # said above
        if len({o["resources"].get(k) for o in done if o["resources"].get(k)}) > 1:
            out.append(f"the counts files were not all made with the same {what} file. " + (
                "Different controls files are accepted only if their control regions are identical, or one of the bundle's named subsets "
                "(extra known-truth regions are then simply missing for some samples)." if what == "controls" else
                "Different sinks change what a fetch reads: do not analyse such files as one cohort."))
    return out


def cmd_estimate(a):
    for p in getattr(a, "fetch_sinks", None) or ():           # a path that is not there would be passed over without a word
        if not Path(p).is_file():
            raise SystemExit(f"ngsdose estimate: --fetch-sinks {p}: no such file")
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
    if _EST["experimental_missing"]:
        print(f"[estimate] WARNING: no experimental resources beside the bundle (looked for {_EST['experimental_missing']}): the sub-options, "
              "the experimental sinks and the candidate units are unknown here, so a fetch made with a sinks BED other than the bundle's "
              "gets its satellite families reported unverified, and positional candidates are skipped. Keep resources/experimental beside "
              "the bundle's directory, or set NGSDOSE_EXPERIMENTAL", file=sys.stderr)
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
    # a sub-option whose class the file never counted was never asked for (an ordinary fetch loads no satellite panel): the
    # table says not_counted, the log does not list it
    subs = set(_EST.get("sub_options") or ())
    na = Counter((k[: -len(".status")], v) for o in done for k, v in o["row"].items() if k.endswith(".status") and v not in (None, "ok", "experimental")
                 and not (v == "not_counted" and k[: -len(".status")] in subs))
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
        rules = {} if a.no_class_rules else resources.Bundle(a.resources).calibration()
        eff_in = json.loads(Path(a.efficiencies).read_text()) if a.efficiencies else None
    except (OSError, ValueError) as e:
        raise SystemExit(f"ngsdose cohort: {e}") from None
    profiles = {} if a.segments else None
    try:
        kar = None if a.no_karyotype else resources.Bundle(a.resources).karyotype()
        if a.karyotype_model:
            if kar is None:
                raise ValueError("--karyotype-model: this bundle says nothing about chromosomes (bundle.json `karyotype`)")
            with io._open(a.karyotype_model) as fh:
                kar.update(model=karyotype.Model.from_json(json.load(fh)), fit_own=False)
    except (OSError, ValueError) as e:
        raise SystemExit(f"ngsdose cohort: {e}") from None
    kout: dict = {}
    try:
        rows, eff, _ = cohort.cohort_table(_estimates(a.estimates), anchors, max_window_sd=a.max_window_sd,
                                           n_profile_pcs=a.profile_pcs, n_control_pcs=a.control_pcs, mp_margin=a.mp_margin,
                                           efficiencies=eff_in, log=lambda m: print(m, file=sys.stderr), rules=rules, profiles=profiles,
                                           karyotype=kar, karyotypes=kout)
    except ValueError as e:
        raise SystemExit(f"ngsdose cohort: {e}") from None
    if a.save_efficiencies:
        Path(a.save_efficiencies).write_text(json.dumps(eff))
    if a.save_karyotype_model:
        if kout.get("model") is None:
            print("[cohort] --save-karyotype-model: no model to save (chromosomes were not read)", file=sys.stderr)
        else:
            _dump(kout["model"].to_json(), a.save_karyotype_model)
    if a.karyotype_table:
        long = karyotype.long_table(kout.get("samples") or [], kout.get("readings") or [])
        if long:
            write_table([{k: r.get(k) for k in karyotype.LONG_COLUMNS} for r in long], a.karyotype_table)
        else:
            Path(a.karyotype_table).write_text("\t".join(karyotype.LONG_COLUMNS) + "\n")
    if a.segments:
        cols = ["sample", "class", "kind", "start", "end", "windows", "state", "mean", "se", "raw", "off_integer", "z", "scale_f", "complete_copies", "call"]
        segs = []
        for cls, p in profiles.items():
            for s in p["samples"]:
                c = p["calls"].get(s)
                if c is None:
                    continue
                same = dict(scale_f=round(c.scale, 3), complete_copies=c.copies, call=c.status)
                # the whole numbers called along the unit, then the stretches that read a fraction of a copy off them
                segs += [dict(sample=s, **{"class": cls}, kind="segment", **g.as_dict(), off_integer=g.off_integer, z=None, **same) for g in c.segments]
                segs += [dict(sample=s, **{"class": cls}, kind="fraction", start=f.start, end=f.end, windows=f.windows, state=f.state, mean=round(f.state + f.offset, 3),
                              se=None, raw=None, off_integer=True, z=f.z, **same) for f in c.fractions]
        segs = [{k: r.get(k) for k in cols} for r in segs]
        if segs:
            write_table(segs, a.segments)
        else:
            print("[cohort] --segments: no class of this bundle has segment calls (calibration.json `segments`)", file=sys.stderr)
            Path(a.segments).write_text("\t".join(cols) + "\n")
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


def _need_columns(header, columns, path, cmd, strict=True) -> list[str]:
    """The requested columns the table lacks. Strict (adjust): the run stops, since its output would lack
    what was asked for. Otherwise (trios, pcsweep) a column that is not there is that column's failure, not
    the run's: it is said once, the caller gives it an NA row with the reason or leaves it out of the sweep,
    and stops only when no column is left."""
    missing = [c for c in dict.fromkeys(columns) if c not in header]
    if missing and strict:
        raise SystemExit(f"ngsdose {cmd}: no column {', '.join(missing)} in {path}")
    if missing:
        print(f"[{cmd}] WARNING: no column {', '.join(missing)} in {path}", file=sys.stderr)
    return missing


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
    """The pedigree's and the trio analysis's warnings, each printed once, as the command's own. Any other
    warning raised inside (numpy's, a library's) is shown as usual - after the recording context has been
    left: shown inside it, each would be recorded again into the very list being walked, without end."""
    caught: list = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", trios.PedigreeWarning)
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
    if a.compare_to and a.compare_to not in header:
        raise SystemExit(f"ngsdose trios: --compare-to: no column {a.compare_to} in {a.table}")
    absent = _need_columns(header, a.columns, a.table, "trios", strict=False)
    if len(absent) == len(set(a.columns)):
        raise SystemExit(f"ngsdose trios: none of the columns ({', '.join(dict.fromkeys(a.columns))}) is in {a.table}")
    out, failed = {}, {c: f"no column {c} in {a.table}" for c in absent}
    with _pedigree_warnings("trios"):
        ped, pop = trios.load_pedigree(a.pedigree, population=a.population)
        if not ped:
            raise SystemExit(f"ngsdose trios: {a.pedigree} holds no complete trio (a child with both parents given); check its layout")
        centre = None if a.no_population_centring else pop
        for col in a.columns:
            if col in failed:
                continue
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
    absent = _need_columns(header, list(a.columns) + list(extra), a.table, "pcsweep", strict=False)
    if any(c in absent for c in extra):
        raise SystemExit(f"ngsdose pcsweep: --truth names a column {a.table} lacks: {', '.join(c for c in extra if c in absent)}")
    a.columns = [c for c in a.columns if c not in absent]     # the others are swept; a column that is not there is left out
    if not a.columns:
        raise SystemExit(f"ngsdose pcsweep: none of the columns is in {a.table}")
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


def _write_stats(a, table, source, held_out):
    try:
        with open(a.stats, "w") as fh:
            sinks.write_stats(fh, table, "\n".join([f"ngsdose sinks --stats: {source}", *(a.stats_note or ())]), held_out)
    except OSError as e:
        sys.exit(f"ngsdose sinks: {e}")
    print(f"[sinks] statistics of {len(table)} intervals ({len({r['class'] for r in table})} classes) -> {a.stats}", file=sys.stderr)


def cmd_sinks(a):
    if a.stats and a.stats == a.out and a.out != "-":
        sys.exit("ngsdose sinks: --stats and -o name the same file")
    if a.held_out and not a.evaluate:
        sys.exit("ngsdose sinks: --held-out describes the scans of --evaluate; learned sinks are measured in-sample")
    if (a.held_out or a.stats_note) and not a.stats:
        sys.exit("ngsdose sinks: --held-out and --stats-note go into the --stats file; pass --stats")
    if a.evaluate:
        try:
            bed = sinks.read_bed(a.evaluate)
            tally = sinks.IntervalTally(bed) if a.stats else None
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
            if tally:
                tally.add(c)
        if tally:
            _write_stats(a, tally.table(), f"{a.evaluate} evaluated on {len(a.counts) - len(bad)} scan(s)", True if a.held_out else None)
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
    if a.stats:                                            # a third pass over the scans, one file at a time
        try:
            table = sinks.interval_stats(a.counts, rows, allow_cut=a.allow_cut)
        except ValueError as e:
            sys.exit(f"ngsdose sinks: {e}")
        _write_stats(a, table, f"the learned intervals, in the {len(a.counts)} training scan(s) (in-sample)", False)


def _capture_class(items) -> dict[str, float]:
    out = {}
    for x in items or ():
        name, eq, v = x.partition("=")
        try:
            out[name] = float(v)
        except ValueError:
            sys.exit(f"ngsdose fetchplan: --capture-class {x}: expected CLASS=FRACTION, e.g. TEL=0.99")
        if not eq or not name or not 0 < out[name] <= 1:
            sys.exit(f"ngsdose fetchplan: --capture-class {x}: expected CLASS=FRACTION with 0 < FRACTION <= 1")
    return out


def _plan_header(a, menu, plan) -> str:
    """Every argument that shapes the plan, as a command line, for the head of PREFIX.sinks.bed and PREFIX.plan.tsv: two plans
    that differ in their intervals, flags or costing differ here too. Defaults are written out; files by path."""
    x = [f"--menu {menu.path}"]
    x += [f"--classes {' '.join(a.classes)}"] if a.classes else []
    x += [f"--preset {' '.join(a.preset)}"] if a.preset else []
    if a.budget_mb is not None:
        x += [f"--budget-mb {a.budget_mb:g}", f"--status {a.status}"] + (["--fill"] if a.fill else [])
    x += [f"--capture {a.capture:g}"] if a.capture is not None else []
    x += [f"--capture-class {' '.join(a.capture_class)}"] if a.capture_class else []
    x += [f"--capture-stat {a.capture_stat}"]
    x += [f"--stats {' '.join(a.stats)}"] if a.stats else []
    x += [f"--sinks {' '.join(a.sinks)}"] if a.sinks else []
    x += [f"--controls {a.controls}"] if a.controls else []
    x += [f"--pad {a.pad}"]
    if a.crai:
        x += [f"--crai ({len(a.crai)} index(es): {' '.join(Path(p).name for p in a.crai)})", f"--contigs {a.contigs}"]
    x += [f"--engine {a.engine}"] if a.engine else []
    x += ["--unmarked-companions"] if a.unmarked_companions else []
    x += [f"--panel-root {a.panel_root}"] if a.panel_root else []
    flags = f"; count flags: {plan.classes_flag}" if plan.classes_flag else ""
    return " ".join(x) + flags


def cmd_fetchplan(a):
    try:
        menu = fetchplan.read_menu(a.menu)
        if a.list:
            print("\n".join(fetchplan.list_menu(menu)))
            return
        if a.capture is not None and not 0 < a.capture <= 1:
            sys.exit(f"ngsdose fetchplan: --capture {a.capture}: expected a fraction in (0, 1]")
        plan = fetchplan.make_plan(menu, classes=a.classes or (), presets=a.preset or (), budget_mb=a.budget_mb, fill=a.fill,
                                   statuses=tuple(x.strip() for x in a.status.split(",") if x.strip()), capture=a.capture,
                                   capture_class=_capture_class(a.capture_class), capture_stat=a.capture_stat, sinks_files=a.sinks or (),
                                   stats_files=a.stats or (), controls=a.controls, pad=a.pad, crais=a.crai or (), contigs=a.contigs,
                                   engine=a.engine, unmarked_companions=a.unmarked_companions)
    except (OSError, ValueError) as e:
        sys.exit(f"ngsdose fetchplan: {e}")
    print("\n".join(fetchplan.table_lines(plan)))
    if a.out:
        try:
            fetchplan.write(plan, a.out, menu, panel_root=a.panel_root, header=_plan_header(a, menu, plan))
        except OSError as e:
            sys.exit(f"ngsdose fetchplan: {e}")
        print(f"[fetchplan] wrote {a.out}.sinks.bed ({len(plan.bed())} intervals), {a.out}.panels.txt, {a.out}.count_flags.txt, "
              f"{a.out}.controls.txt (-c {Path(plan.controls_fasta).name}), {a.out}.scan_panels.txt, {a.out}.plan.tsv", file=sys.stderr)


def cmd_controlsets(a):
    """The bundle's control sets: which of its single-copy regions a fetch reads."""
    try:
        res = resources.Bundle(a.resources)
        sets = fetchplan.control_sets(res.dir / "controls.bed")
        if not sets:
            sys.exit(f"ngsdose control-sets: no controls.bed in {res.dir}")
        print("set\tregions\tMb\tcontrols\ttruth_regions\tkaryotype_pieces\tfasta")
        for name in ["all"] + sorted(n for n in sets if n != "all"):
            rows = sinks.read_bed(sets[name])
            fa = res.controls if name == "all" else fetchplan.controls_fasta(sets[name])
            if a.write and not fa.exists():
                n = fetchplan.subset_fasta(res.controls, [f"{c}:{s0}-{e0}" for c, s0, e0, _ in rows], fa)
                print(f"[control-sets] wrote {fa} ({n} regions, cut from {res.controls.name})", file=sys.stderr)
            role = lambda r: r[3] or "control"
            print("\t".join(map(str, (name, len(rows), f"{sum(e0 - s0 for _, s0, e0, _ in rows) / 1e6:.2f}", sum(1 for r in rows if role(r) == "control"),
                                      sum(1 for r in rows if role(r).startswith("test:") and not role(r).startswith("test:karyotype")),
                                      sum(1 for r in rows if role(r).startswith("test:karyotype")),
                                      fa.name if fa.exists() else "not written: `ngsdose control-sets --write`, or `ngsdose fetchplan --controls NAME` cuts it for a plan"))))
    except (OSError, ValueError) as e:
        sys.exit(f"ngsdose control-sets: {e}")


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
    e.add_argument("--fetch-sinks", nargs="+", metavar="BED", help="sinks BEDs fetches were made with (e.g. a fetchplan PREFIX.sinks.bed), "
                   "known by the sha256 the counts record: the sub-options (resources/experimental/subsets: DXZ1, DYZ3, DYZ1, DYZ2) of a "
                   "fetch are reported only when its BED is known and holds their intervals, and their classes (aSatHOR, HSat3, "
                   "HSat1B) only when it is known (else status unverified: the fetch may have read them at the sub-options only). "
                   "The bundle's sinks.bed and the BEDs of resources/experimental are known without it")
    e.set_defaults(fn=cmd_estimate)

    c = sub.add_parser("cohort", help="calibrate window efficiencies across samples")
    c.add_argument("estimates", nargs="+")
    c.add_argument("-t", "--table", default="-")
    c.add_argument("--efficiencies", help="apply a saved efficiency table instead of fitting one")
    c.add_argument("--save-efficiencies")
    c.add_argument("-r", "--resources", default=None)
    c.add_argument("--gc-rule-anchors", action="store_true")
    c.add_argument("--no-class-rules", action="store_true",
                   help="ignore the bundle's calibration.json: every class is levelled on all its windows, scaled by its anchors, and no copy "
                        "states are called along the unit")
    c.add_argument("--segments", metavar="TSV",
                   help="write the integer copy states called along the unit, one row per genome and segment (kind `segment`; `raw` is the segment's mean as "
                        "the reads give it, before the genome's scale and lean), and one per stretch that reads a fraction of a copy off them (kind `fraction`, "
                        "with its z against the cohort), for the classes whose rules ask for them")
    c.add_argument("--no-karyotype", action="store_true", help="do not read chromosomes in copies (the karyotype columns)")
    c.add_argument("--karyotype-model", metavar="JSON",
                   help="read the chromosomes against this saved model of the single-copy regions instead of learning one on this cohort "
                        "(a cohort of fewer than 50 genomes takes the bundle's model without it): a genome's reading then does not depend on "
                        "which other genomes are given")
    c.add_argument("--save-karyotype-model", metavar="JSON", help="save the model the chromosomes were read against (.gz for gzip)")
    c.add_argument("--karyotype-table", metavar="TSV",
                   help="write every genome's chromosomes, one row each (copies, SE, whole number, distance from it, z, status), and one row per "
                        "stretch where a chromosome holds more than one level (an arm, a segment)")
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
    k.add_argument("--stats", metavar="FILE", help="also write, per interval (of the learned sinks, or of the --evaluate BED): its share of "
                   "the class's reads (median and 10th percentile over the scans), the reads of any class in it (the 'all' counts of the "
                   "class's placement bins there, median) and, per class in order of decreasing share per read, the cumulative capture "
                   "(median and 10th percentile): the capture curve `ngsdose fetchplan --capture` trims by. For learned sinks it is "
                   "in-sample; with --evaluate on held-out scans it is not (say so with --held-out). Also the largest share an interval held "
                   "in any one scan (share_max), which bounds the capture of intervals kept in another order (fetchplan --crai)")
    k.add_argument("--held-out", action="store_true", help="with --evaluate and --stats: none of these scans was used to learn the BED. "
                   "The statistics file says so ('# held-out: yes'), and `ngsdose fetchplan` then reports its expected capture as held-out")
    k.add_argument("--stats-note", action="append", metavar="TEXT", help="a line of provenance for the head of the --stats file (repeatable)")
    k.set_defaults(fn=cmd_sinks)

    f = sub.add_parser("fetchplan", help="choose what a fetch reads: classes by name, preset or byte budget, with capture targets",
                       description="Select options from the fetch menu (resources/fetch_menu.tsv) and write what `ngs-dose count -m fetch` "
                                   "takes: PREFIX.sinks.bed (--sinks), PREFIX.panels.txt (one -p each), PREFIX.count_flags.txt, PREFIX.controls.txt "
                                   "(the -c FASTA matching the control regions costed); plus "
                                   "PREFIX.scan_panels.txt (panels for the whole-file scans, candidates included) and PREFIX.plan.tsv. "
                                   "The control regions are always read. Sinks come from scans (`ngsdose sinks`): a candidate class, "
                                   "which has none yet, is selected for scanning only.")
    f.add_argument("--menu", help="the fetch menu (default: fetch_menu.tsv beside $NGSDOSE_RESOURCES, else resources/fetch_menu.tsv)")
    f.add_argument("--list", action="store_true", help="print the menu's options and presets and stop")
    f.add_argument("--classes", nargs="+", metavar="OPTION", help="options by name (a class, a named subset of a class's intervals such as DXZ1, or 'unmapped')")
    f.add_argument("--preset", nargs="+", metavar="PRESET", help="options by preset (see --list)")
    f.add_argument("--budget-mb", type=float, help="add options in tier order (A to D), cheapest first within a tier, until the next would take "
                   "the plan (controls included) past this many MB (1e6 bytes, median over the --crai indexes); from the --classes and "
                   "--preset options if given, otherwise from the whole menu. Needs --crai")
    f.add_argument("--fill", action="store_true", help="for --budget-mb: skip an option that does not fit and go on adding later ones that "
                   "do (default: stop at the first that does not fit, so that no lower tier displaces a higher one)")
    f.add_argument("--status", default="shipped,experimental", help="for --budget-mb: the statuses to draw from (default shipped,experimental)")
    f.add_argument("--capture", type=float, help="per class, keep the intervals of highest share per read until the expected capture reaches "
                   "this fraction, drop the rest (needs statistics: the menu's stats column or --stats; a class without them keeps all its "
                   "intervals, and a named subset is fetched whole). With --crai, per byte of the interval's own fetch instead (the CRAM "
                   "slices it overlaps, with their containers' compression headers: what dropping it saves): an interval of low median "
                   "share on costly slices (a pile-up bin, a multi-reference decoy slice) ranks behind the intervals worth their bytes. "
                   "The expected capture is then a lower bound (the capture of all the class's intervals "
                   "less the largest share each dropped interval held in any one scan), so an interval is dropped only while its largest "
                   "share in any scan still leaves the target reached (TEL keeps the chr2:32.91 Mb pile-up bin above about 0.979), and "
                   "when keeping per read reaches the target with fewer bytes, that is kept. Classes share intervals: one that another "
                   "selected class keeps is read anyway, so it is kept for every class that has it (with --crai: any interval that adds "
                   "no bytes to the plan's fetches, and each class is trimmed again with those intervals free). The table reports per "
                   "option the MB its trimming saved of the plan (mb_saved) and the capture lost")
    f.add_argument("--capture-class", nargs="+", metavar="CLASS=FRACTION", help="a capture target for one class, e.g. TEL=0.99")
    f.add_argument("--capture-stat", choices=("p10", "median"), default="p10",
                   help="which capture of the statistics' scans the target applies to: the 10th percentile (default) or the median")
    f.add_argument("--stats", nargs="+", metavar="FILE", help="interval statistics from `ngsdose sinks --stats` (later files win per class)")
    f.add_argument("--sinks", nargs="+", metavar="BED", help="take every class's intervals from these sinks BEDs (learned for this pipeline) "
                   "instead of the files the menu names")
    f.add_argument("--controls", metavar="SET", help="the regions a fetch reads for the denominator, the known truths and the chromosomes: the name of one of "
                   "the bundle's control sets (controls.NAME.bed beside its controls: `all`, `base`, `lite200`, `karyotype`, `screen` in the GRCh38 bundle; "
                   "README, 'Choosing what a fetch reads') or a BED (default: the menu's controls row). The fetch must use the controls FASTA of exactly "
                   "these regions, which controls.txt names: controls.NAME.fa.gz where the bundle ships it, else PREFIX.controls.fa.gz, cut from the "
                   "bundle's controls.fa.gz when the plan is written")
    f.add_argument("--pad", type=int, default=600, help="padding of the control regions, as `ngs-dose count --pad` (default 600, at least 400); another value goes to "
                   "count_flags.txt so the fetch reads the regions costed here")
    f.add_argument("--crai", nargs="+", metavar="CRAI", help="CRAM indexes to cost the plan on (median over them); each from the same place "
                   "as its CRAM. The plan is priced as `ngs-dose count -m fetch` reads it: one indexed fetch per run of touching or "
                   "overlapping intervals, each decoding every CRAM slice that overlaps it (with its container's compression header), "
                   "so a slice under several runs is decoded once per run; cum_mb_floor is the plan with every slice decoded once. Over "
                   "HTTP(S) an engine from 0.4.0 on (`count --transport ranges`, its default) moves about cum_mb_floor: it asks for the "
                   "containers its queries read, each once; an older one, or --transport htslib, moves several times the MB")
    f.add_argument("--contigs", metavar="FILE", help="the contigs of the CRAMs in header order: the .fai or .dict of their reference, or "
                   "`samtools view -H` output")
    f.add_argument("--engine", default="ngs-dose", help="the ngs-dose that will run the fetch (default: ngs-dose on PATH). When the plan "
                   "needs count flags, fetchplan runs `ENGINE count --help` and refuses a plan the engine cannot run (--pad other than "
                   "600 on fae1124). When the loaded panels define classes the plan does not select: an engine that takes --classes gets "
                   "'--classes=...' in count_flags.txt (only the selected classes are counted and reported); one that takes only "
                   "--allow-missing-sinks (7772e32) gets that flag (the others are listed as incomplete); one that takes neither "
                   "(fae1124) is refused, see --unmarked-companions. When no engine is found, '--allow-missing-sinks' is written, "
                   "which fae1124 refuses at the fetch. The notes and plan.tsv say which was written")
    f.add_argument("--unmarked-companions", action="store_true", help="with an engine that takes neither --classes nor "
                   "--allow-missing-sinks (fae1124): write the plan anyway. The classes the loaded panels define but the plan does not "
                   "select are then counted only inside the plan's intervals and not marked incomplete in the counts; `ngsdose "
                   "estimate --fetch-sinks PREFIX.sinks.bed` marks them")
    f.add_argument("--panel-root", metavar="DIR", help="write panel paths under DIR (e.g. the resources directory inside a container) instead "
                   "of the menu's directory")
    f.add_argument("-o", "--out", metavar="PREFIX", help="write PREFIX.sinks.bed, .panels.txt, .count_flags.txt, .controls.txt, .scan_panels.txt, .plan.tsv "
                   "(without it, only the table is printed)")
    f.set_defaults(fn=cmd_fetchplan)

    k = sub.add_parser("control-sets", help="list the bundle's control sets (which single-copy regions a fetch reads), and write their FASTAs",
                       description="The bundle's controls.bed holds every single-copy region it has; controls.NAME.bed beside it is a set of them that a "
                                   "fetch can read instead (`ngsdose fetchplan --controls NAME`). The engine takes a set as a FASTA of exactly its "
                                   "regions (-c). --write cuts each set's FASTA that is not there yet from the bundle's controls.fa.gz, into the "
                                   "bundle's directory (an image or an install does this once; no reference genome is needed).")
    k.add_argument("-r", "--resources", default=None, help="resource bundle directory (default: packaged GRCh38)")
    k.add_argument("--write", action="store_true", help="write the FASTAs of the sets that lack one")
    k.set_defaults(fn=cmd_controlsets)

    s = sub.add_parser("selftest", help="simulation-based check of the estimator; needs no data")
    s.set_defaults(fn=cmd_selftest)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
