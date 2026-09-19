"""Command line for the modelling layer: `ngsdose estimate | cohort | adjust | trios | sinks | selftest`."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

from . import __version__, cohort, estimate, io, resources, sinks, trios


class _Enc(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def _dump(obj, path):
    data = json.dumps(obj, cls=_Enc, allow_nan=True)
    if str(path) == "-":
        sys.stdout.write(data + "\n")
    elif str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(data)
    else:
        Path(path).write_text(data)


def _load_result(path) -> dict:
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        return json.load(fh)


def summary_row(r: dict) -> dict:
    row = dict(sample=r["sample"], mode=r["mode"], depth=round(r["depth_equiv"], 3), read_length=r["read_length"],
               insert_median=r["insert_median"], gc_L=r["gc_L"], ctrl_dup_frac=round(r["ctrl_dup_frac"], 4),
               gc_rel_35=r["gc_rel"].get("35"), gc_rel_65=r["gc_rel"].get("65"), gc_curve_max_se=round(r["gc_curve_max_se"], 4),
               ctrl_region_sd=None if not r["control_qc"] else round(r["control_qc"]["region_log_mad_sd"], 4),
               flagged_chromosomes=None if not r["control_qc"] else ",".join(r["control_qc"]["flagged_chromosomes"]))
    for label, t in r.get("truth_regions", {}).items():
        # known-truth sets are scored against their answer; dosage sets (chrM, chrEBV) are copies per cell
        row[f"{label}.copies" if t.get("role") == "dosage" else f"truth.{label}"] = round(t["cn"], 4)
    if r.get("eof_marker"):
        row["eof_marker"] = r["eof_marker"]
    for name, c in r["classes"].items():
        if c["kind"] == "positional":
            row[f"{name}.cn_single"] = round(c["cn"], 3)
            row[f"{name}.cn_anchor"] = round(c["cn_anchor"], 2)
            row[f"{name}.cn_all"] = round(c["cn_all"], 2)
            row[f"{name}.cn_median"] = round(c["cn_median"], 2)
            row[f"{name}.window_log_sd"] = round(c["window_log_sd"], 4)
            row[f"{name}.cn_all_flat"] = round(c["cn_all_flat"], 2)
            for fn, fv in c["features"].items():
                row[f"{name}.{fn}"] = round(fv["cn"], 2)
                row[f"{name}.{fn}.flat"] = round(fv["cn_flat"], 2)
        else:
            row[f"{name}.mass_Mb"] = round(c["mass_Mb"], 4)
    return row


def write_table(rows: list[dict], path):
    cols: list[str] = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    fh = sys.stdout if str(path) == "-" else open(path, "w", newline="")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="NA")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("NA" if v is None or (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in r.items()})
    if fh is not sys.stdout:
        fh.close()


_EST: dict = {}


def _estimate_init(a):
    res = resources.Bundle(a.resources)
    _EST.update(args=a, res=res, panel=io.load_panel(res.panel), units=res.units(), feats=res.features(),
                anchors={} if a.gc_rule_anchors else res.anchors(), tables={}, lengths=res.contig_lengths())


def _estimate_one(path):
    a, res = _EST["args"], _EST["res"]
    counts = io.load_counts(path)
    L = estimate.nearest_table(counts, a.L)["l"]
    if not a.no_control_qc and L not in _EST["tables"]:
        _EST["tables"][L] = res.region_tables(L)
    r = estimate.estimate_sample(counts, _EST["panel"], _EST["units"], _EST["feats"], region_tables=_EST["tables"].get(L), L=a.L,
                                 window=a.window, min_kmers=a.min_kmers, anchors=_EST["anchors"],
                                 contig_lengths=_EST["lengths"], regions=None if a.no_control_qc else res.regions())
    if r["eof_marker"] == "absent":
        print(f"[estimate] WARNING: {r['sample']} was counted from a file without an end-of-file marker (truncated input)", file=sys.stderr)
    r["resources"] = {k: counts.get(k) for k in ("panel_sha256", "controls_sha256", "sinks_sha256")}
    if a.outdir:
        _dump(r, Path(a.outdir) / f"{r['sample']}.estimate.json.gz")
    return summary_row(r), r["resources"]


def cmd_estimate(a):
    if a.outdir:
        Path(a.outdir).mkdir(parents=True, exist_ok=True)
    if a.jobs > 1 and len(a.counts) > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.jobs, initializer=_estimate_init, initargs=(a,)) as pool:
            out = pool.map(_estimate_one, a.counts, chunksize=4)
    else:
        _estimate_init(a)
        out = [_estimate_one(p) for p in a.counts]
    # a cohort assembled from runs against different bundles is not one cohort
    differ = [k for k in ("panel_sha256", "controls_sha256", "sinks_sha256")
              if len({json.dumps(h.get(k)) for _, h in out if h.get(k)}) > 1]
    if differ:
        print(f"[estimate] WARNING: the counts files were not all made with the same {', '.join(k[:-7] for k in differ)} file(s). "
              "A different panel or different sinks change what is counted: do not analyse such files as one cohort. Different "
              "controls files are accepted only if their control regions are identical (extra known-truth regions are then "
              "simply missing for some samples).", file=sys.stderr)
    write_table([row for row, _ in out], a.table)


def cmd_cohort(a):
    # one estimate file at a time: only the window vectors, the control residuals and the summary
    # row of each sample are kept, so a cohort of thousands fits in a few hundred MB
    anchors = {} if a.gc_rule_anchors else resources.Bundle(a.resources).anchors()
    rows, order, win, layout, ctrl = {}, [], {}, {}, []
    for path in a.estimates:
        r = _load_result(path)
        rows[r["sample"]] = summary_row(r)
        order.append(r["sample"])
        for cls, v in r["classes"].items():
            if v["kind"] != "positional":
                continue
            y, gc, starts, ends = cohort.sample_windows(r, cls)
            if cls in layout and len(layout[cls][0]) != len(starts):
                raise SystemExit(f"{path}: {cls} was estimated with a different window layout")
            layout.setdefault(cls, (starts, ends))
            win.setdefault(cls, []).append((y, gc))
        ctrl.append((r.get("control_qc") or {}).get("region_log_ratio"))
    eff = {}
    for cls, per in win.items():
        if len(per) != len(order):
            continue
        a_fixed = None
        if a.efficiencies:
            tab = json.loads(Path(a.efficiencies).read_text())
            a_fixed = np.array([np.nan if v is None else v for v in tab[cls]["a"]], float)
        Y = np.array([p[0] for p in per])
        gc = np.nanmedian(np.array([p[1] for p in per]), axis=0)
        cal = cohort.calibrate_matrix(order, Y, layout[cls][0], layout[cls][1], gc, a_fixed=a_fixed,
                                      max_window_sd=a.max_window_sd, anchors=anchors.get(cls))
        scores, var = cohort.profile_pcs(cal, a.profile_pcs) if len(order) > a.profile_pcs + 2 else (None, None)
        for i, s in enumerate(cal.samples):
            rows[s][f"{cls}.cn"] = round(float(np.exp(cal.c[i])), 2)
            rows[s][f"{cls}.cn_se_rel"] = round(float(cal.c_se[i]), 5)
            rows[s][f"{cls}.profile_sd"] = round(float(cal.resid_sd[i]), 4)
            if scores is not None:
                for k in range(scores.shape[1]):
                    rows[s][f"{cls}.profilePC{k + 1}"] = round(float(scores[i, k]), 5)
        eff[cls] = dict(start=cal.window_start.tolist(), gc=[round(float(g), 4) for g in cal.window_gc],
                        anchor=cal.anchor.tolist(), a=[None if np.isnan(v) else round(float(v), 5) for v in cal.a],
                        window_sd=[None if np.isnan(v) else round(float(v), 5) for v in cal.window_sd],
                        n_samples=len(cal.samples))
    ok = len(order) >= 5 * max(a.control_pcs, 1) and all(x is not None for x in ctrl) and len({len(x) for x in ctrl}) == 1
    cp = cohort.control_pcs(np.array(ctrl, float), a.control_pcs) if ok else None
    if cp is not None:
        scores, var = cp
        for i, s in enumerate(order):
            for k in range(scores.shape[1]):
                rows[s][f"ctrlPC{k + 1}"] = round(float(scores[i, k]), 5)
        print("[cohort] control-region PCs, variance explained: " + " ".join(f"{v:.3f}" for v in var), file=sys.stderr)
    if a.save_efficiencies:
        Path(a.save_efficiencies).write_text(json.dumps(eff))
    write_table([rows[s] for s in order], a.table)


def _read_values(path, column):
    vals = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                vals[r["sample"]] = float(r[column])
            except (ValueError, KeyError):
                pass
    if not vals:
        raise SystemExit(f"no numeric values in column {column!r} of {path}")
    return vals


def cmd_trios(a):
    ped, pop = trios.load_pedigree(a.pedigree)
    centre = None if a.no_population_centring else pop
    out = {}
    for col in a.columns:
        out[col] = trios.transmission(_read_values(a.table, col), ped, centre, n_perm=a.perm)
    ci = lambda t, k: ("(%.3f,%.3f)" % t[k + "_ci95"]) if k + "_ci95" in t else "NA"
    hdr = ("column", "n_trios", "R_midparent", "CI95", "slope", "se", "R_single", "R_mendel", "spousal_r", "CI95", "error_cv", "perm_null")
    print("\t".join(hdr))
    for col, t in out.items():
        print("\t".join([col, str(t["n_trios"]), f"{t['reliability_midparent']:.3f}", ci(t, "reliability_midparent"),
                         f"{t['midparent_slope']:.3f}", f"{t['midparent_slope_se']:.3f}", f"{t['reliability_single_parent']:.3f}",
                         f"{t['reliability_mendel']:.3f}", f"{t['spousal_r']:+.3f}", ci(t, "spousal_r"), f"{t['error_cv']:.4f}",
                         f"{t.get('perm_null_mean', float('nan')):+.3f}"]))
    if a.compare_to:
        base = _read_values(a.table, a.compare_to)
        print("\n# paired family bootstrap: R(column) - R(%s)" % a.compare_to)
        print("column\tdelta_R\tCI95\tP(column better)")
        for col in a.columns:
            if col == a.compare_to:
                continue
            try:
                c = trios.compare(_read_values(a.table, col), base, ped, centre)
            except ValueError as err:
                print(f"{col}\tNA\tNA\tNA\t# {err}")
                continue
            out[col]["vs_" + a.compare_to] = c
            print(f"{col}\t{c['delta']:+.4f}\t({c['ci95'][0]:+.4f},{c['ci95'][1]:+.4f})\t{c['p_a_better']:.3f}")
    if a.json:
        _dump(out, a.json)


def load_pcs(path, n_pc, strip=(".by1000.", ".")):
    """NGS-PCA svd.pcs.txt -> {sample: PCs}. NGS-PCA names samples after the mosdepth file, which
    leaves a suffix such as `.by1000.`; it is removed so that ids match the alignment's SM tag."""
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
            pcs[name] = [float(x) for x in p[1:1 + n_pc]]
    return pcs


def cmd_adjust(a):
    with open(a.table) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if a.pcs:
        pcs = load_pcs(a.pcs, a.n_pc)
    else:                                                  # the cohort table's own control-region PCs
        cols = [f"ctrlPC{k + 1}" for k in range(a.n_pc)]
        if not all(c in rows[0] for c in cols):
            raise SystemExit(f"no --pcs given and the table lacks {cols[-1]}; run `ngsdose cohort --control-pcs {a.n_pc}`")
        pcs = {r["sample"]: [float(r[c]) for c in cols] for r in rows if all(r[c] not in ("NA", "") for c in cols)}
    keep = [r for r in rows if r["sample"] in pcs]
    if len(keep) < a.n_pc + 10:
        raise SystemExit(f"only {len(keep)} samples have coverage PCs; need more than n_pc + 10")
    X = np.array([pcs[r["sample"]] for r in keep])
    for col in a.columns:
        y = np.array([float(r[col]) if r.get(col, "NA") not in ("NA", "") else np.nan for r in keep])
        adj, r2 = cohort.adjust_for_covariates(y, X, log=not a.no_log)
        for r, v in zip(keep, adj):
            r[f"{col}.adj"] = "NA" if not np.isfinite(v) else round(float(v), 4)
        n, k = int(np.isfinite(adj).sum()), X.shape[1]
        # k regressors explain k/(n-1) of pure noise: report what is left after that
        r2_adj = 1 - (1 - r2) * (n - 1) / max(n - k - 1, 1)
        print(f"[adjust] {col}: {k} PCs explain {100 * r2:.1f}% of the variance of {'log ' if not a.no_log else ''}{col} "
              f"(n={n}; {100 * k / max(n - 1, 1):.1f}% expected by chance; adjusted R2 {100 * r2_adj:.1f}%)", file=sys.stderr)
    write_table(keep, a.out)


def cmd_sinks(a):
    if a.evaluate:
        bed = [(p[0], int(p[1]), int(p[2]), p[3].strip()) for p in (l.split("\t") for l in open(a.evaluate) if l.strip())]
        print("sample\tclass\tscan_reads\tcaptured\tfraction")
        for path in a.counts:                              # one file at a time: a cohort of scans fits
            c = io.load_counts(path)
            for cls, (tot, inside) in sinks.capture(c, bed).items():
                print(f"{c['sample']}\t{cls}\t{tot}\t{inside}\t{inside / max(tot, 1):.5f}")
        return
    rows, stats = sinks.learn(a.counts, a.min_frac, a.pad, a.min_reads)
    with (sys.stdout if a.out == "-" else open(a.out, "w")) as fh:
        for contig, s0, e0, cls in rows:
            fh.write(f"{contig}\t{s0}\t{e0}\t{cls}\n")
    for cls, per in stats.items():
        v = sorted(per.values())
        print(f"[sinks] {cls}: {sum(1 for r in rows if r[3] == cls)} intervals, {sum(r[2] - r[1] for r in rows if r[3] == cls):,} bp; "
              f"capture min {v[0]:.5f} median {v[len(v) // 2]:.5f} over {len(v)} scan(s)", file=sys.stderr)


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
    e.add_argument("-o", "--outdir", help="write one <sample>.estimate.json.gz per input")
    e.add_argument("-t", "--table", default="-", help="summary TSV (default stdout)")
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
    c.add_argument("--control-pcs", type=int, default=10,
                   help="components of the control regions' residual depth to report (needs >= 5 samples per component)")
    c.set_defaults(fn=cmd_cohort)

    t = sub.add_parser("trios", help="transmission reliability from a pedigree")
    t.add_argument("table")
    t.add_argument("-p", "--pedigree", required=True)
    t.add_argument("-c", "--columns", nargs="+", required=True)
    t.add_argument("--perm", type=int, default=1000)
    t.add_argument("--compare-to", help="baseline column: paired bootstrap of the reliability difference of every other column against it")
    t.add_argument("--no-population-centring", action="store_true")
    t.add_argument("--json")
    t.set_defaults(fn=cmd_trios)

    j = sub.add_parser("adjust", help="residualise estimates on NGS-PCA coverage PCs")
    j.add_argument("table")
    j.add_argument("--pcs", help="NGS-PCA svd.pcs.txt (default: the ctrlPC columns written by `ngsdose cohort`)")
    j.add_argument("--n-pc", type=int, default=20)
    j.add_argument("-c", "--columns", nargs="+", required=True)
    j.add_argument("--no-log", action="store_true")
    j.add_argument("-o", "--out", default="-")
    j.set_defaults(fn=cmd_adjust)

    k = sub.add_parser("sinks", help="learn fetch-mode sink intervals from scan-mode counts")
    k.add_argument("counts", nargs="+")
    k.add_argument("-o", "--out", default="-")
    k.add_argument("--min-frac", type=float, default=1e-5, help="keep placement bins holding at least this fraction of a class's reads")
    k.add_argument("--evaluate", metavar="BED", help="instead of learning, report how much of each class an existing sinks BED captures")
    k.add_argument("--min-reads", type=int, default=25)
    k.add_argument("--pad", type=int, default=1000)
    k.set_defaults(fn=cmd_sinks)

    s = sub.add_parser("selftest", help="simulation-based check of the estimator; needs no data")
    s.set_defaults(fn=cmd_selftest)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
