"""Command line for the modelling layer: `ngsdose estimate | cohort | adjust | trios | sinks | selftest`."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

from . import __version__, cohort, estimate, io, pcselect, resources, sinks, trios


from .tables import Encoder as _Enc, dump as _dump, load_result as _load_result, summary_row, write_table  # noqa: E402,F401


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
    eff_in = json.loads(Path(a.efficiencies).read_text()) if a.efficiencies else None
    rows, eff, _ = cohort.cohort_table((_load_result(p) for p in a.estimates), anchors, max_window_sd=a.max_window_sd,
                                       n_profile_pcs=a.profile_pcs, n_control_pcs=a.control_pcs, mp_margin=a.mp_margin,
                                       efficiencies=eff_in, log=lambda m: print(m, file=sys.stderr))
    if a.save_efficiencies:
        Path(a.save_efficiencies).write_text(json.dumps(eff))
    write_table(rows, a.table)


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
    out, failed = {}, {}
    for col in a.columns:
        try:
            out[col] = trios.transmission(_read_values(a.table, col), ped, centre, n_perm=a.perm)
        except (ValueError, KeyError, SystemExit) as err:  # too few complete trios, or no such column: the other columns still get their row
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
        base = _read_values(a.table, a.compare_to)
        print("\n# paired family bootstrap: R(column) - R(%s)" % a.compare_to)
        print("column\tdelta_R\tCI95\tP(column better)")
        for col in a.columns:
            if col == a.compare_to or col in failed:
                continue
            try:
                c = trios.compare(_read_values(a.table, col), base, ped, centre)
            except (ValueError, KeyError, SystemExit) as err:
                print(f"{col}\tNA\tNA\tNA\t# {err}")
                continue
            out[col]["vs_" + a.compare_to] = c
            print(f"{col}\t{c['delta']:+.4f}\t({c['ci95'][0]:+.4f},{c['ci95'][1]:+.4f})\t{c['p_a_better']:.3f}")
    if a.json:
        _dump({**out, **{col: dict(error=msg) for col, msg in failed.items()}}, a.json)


def load_pcs(path, n_pc=None, strip=(".by1000.", ".")):
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
            pcs[name] = [float(x) for x in (p[1:] if n_pc is None else p[1:1 + n_pc])]
    return pcs


def _clamp(asked: int, available: int, what: str):
    """A number of PCs asked for outright is used as far as it can be: a pipeline should not die at its
    last step because one of two PC sets is smaller than the other, but it has to say so."""
    if asked > available:
        print(f"WARNING: --n-pc {asked}, but there are {available} {what}: using {available}", file=sys.stderr)
        return available, f"--n-pc {asked}, clamped to the {available} available"
    return asked, f"--n-pc {asked}"


def _pcs_and_choice(a, rows):
    """{sample: every available PC}, the number to use, and how it was chosen. `--n-pc mp` (the
    default) counts the components above the Marchenko-Pastur edge of the noise bulk: for NGS-PCA's
    PCs from the singular values and bin list it writes beside svd.pcs.txt, for the table's own
    control-region PCs from the count `ngsdose cohort` recorded (ctrlPC_mp)."""
    auto = str(a.n_pc).lower() == "mp"
    if a.pcs:
        pcs = load_pcs(a.pcs)
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
    cols = [c for c in rows[0] if re.fullmatch(r"ctrlPC\d+", c)]
    cols.sort(key=lambda c: int(c[6:]))
    if not cols:
        raise SystemExit("no --pcs given and the table has no ctrlPC columns; run `ngsdose cohort` with --control-pcs")
    pcs = {r["sample"]: [float(r[c]) for c in cols] for r in rows if all(r.get(c, "NA") not in ("NA", "") for c in cols)}
    if not auto:
        return pcs, *_clamp(int(a.n_pc), len(cols), "control-region PCs the table carries (`ngsdose cohort --control-pcs N` writes more)")
    mp = next((r["ctrlPC_mp"] for r in rows if r.get("ctrlPC_mp", "NA") not in ("NA", "")), None)
    if mp is None:
        raise SystemExit("--n-pc mp: the table does not record a Marchenko-Pastur count (ctrlPC_mp); re-run `ngsdose cohort`, or give --n-pc N")
    return pcs, min(int(float(mp)), len(cols)), "the Marchenko-Pastur count recorded by `ngsdose cohort` (ctrlPC_mp)"


def cmd_adjust(a):
    with open(a.table) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    pcs, k, how = _pcs_and_choice(a, rows)
    keep = [r for r in rows if r["sample"] in pcs]
    if len(keep) < k + 10:
        raise SystemExit(f"only {len(keep)} samples have coverage PCs; need more than n_pc + 10 = {k + 10}")
    print(f"[adjust] regressing out {k} PCs ({how})", file=sys.stderr)
    if k > len(keep) / 10:
        print(f"[adjust] WARNING: {k} PCs for {len(keep)} samples is more than one regressor per ten samples - a count chosen on a larger "
              "cohort's SVD does not carry over to a subset; consider --n-pc N (and see `ngsdose pcsweep`)", file=sys.stderr)
    X = np.array([pcs[r["sample"]][:k] for r in keep]).reshape(len(keep), k)
    for col in a.columns:
        y = np.array([float(r[col]) if r.get(col, "NA") not in ("NA", "") else np.nan for r in keep])
        adj, r2 = cohort.adjust_for_covariates(y, X, log=not a.no_log) if k else (y.copy(), 0.0)
        for r, v in zip(keep, adj):
            r[f"{col}.adj"] = "NA" if not np.isfinite(v) else round(float(v), 4)
        n = int(np.isfinite(adj).sum())
        # k regressors explain k/(n-1) of pure noise: report what is left after that
        r2_adj = 1 - (1 - r2) * (n - 1) / max(n - k - 1, 1)
        print(f"[adjust] {col}: {k} PCs explain {100 * r2:.1f}% of the variance of {'log ' if not a.no_log else ''}{col} "
              f"(n={n}; {100 * k / max(n - 1, 1):.1f}% expected by chance; adjusted R2 {100 * r2_adj:.1f}%)", file=sys.stderr)
    write_table(keep, a.out)


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


def cmd_pcsweep(a):
    with open(a.table) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    pcs, k_mp, how = _pcs_and_choice(a, rows)
    keep = [r for r in rows if r["sample"] in pcs]
    n_avail = min(len(pcs[r["sample"]]) for r in keep)
    max_pc = min(n_avail, a.max_pc if a.max_pc else max(2 * k_mp, 20), max(len(keep) // 5, 1))
    P = np.array([pcs[r["sample"]][:max_pc] for r in keep]).reshape(len(keep), max_pc)
    samples = [r["sample"] for r in keep]
    num = lambda col: np.array([float(r[col]) if r.get(col, "NA") not in ("NA", "") else np.nan for r in keep])
    sex = _sex_from(keep)
    male = np.array([sex.get(s) == "M" for s in samples]); female = np.array([sex.get(s) == "F" for s in samples])
    # the known truths the bundle builds in, where the table has them; --truth COLUMN=VALUE adds others
    truths = {}
    if "truth.auto" in keep[0]:
        truths["truth.auto"] = np.full(len(keep), 2.0)
    if "truth.chrX" in keep[0]:
        truths["truth.chrX"] = np.where(male, 1.0, np.where(female, 2.0, np.nan))
    if "truth.chrY" in keep[0]:
        truths["truth.chrY"] = np.where(male, 1.0, np.nan)               # a female's truth is zero, which has no log
    for cand in ("DJ.cn", "DJ.cn_single"):
        if cand in keep[0]:
            truths[cand] = np.full(len(keep), 10.0)
            break
    for spec in a.truth or []:
        col, val = spec.rsplit("=", 1)
        truths[col] = np.full(len(keep), float(val))
    table = {c: num(c) for c in list(truths) + list(a.columns)}
    trio_list = population = None
    if a.pedigree:
        trio_list, population = trios.load_pedigree(a.pedigree)
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
        bed = [(p[0], int(p[1]), int(p[2]), p[3].strip()) for p in (l.split("\t") for l in open(a.evaluate) if l.strip())]
        print("sample\tclass\tscan_reads\tcaptured\tfraction")
        for path in a.counts:                              # one file at a time: a cohort of scans fits
            c = io.load_counts(path)
            if sinks.looks_cut(c):
                print(f"# {c['sample']}: {100 * sinks.cut_share(c):.0f}% of its primary reads lie in the control regions - a cut along a fetch plan, "
                      "not a whole-file scan; its capture says nothing", file=sys.stderr)
            for cls, (tot, inside) in sinks.capture(c, bed).items():
                print(f"{c['sample']}\t{cls}\t{tot}\t{inside}\t{inside / max(tot, 1):.5f}")
        return
    try:
        rows, stats = sinks.learn(a.counts, a.min_frac, a.pad, a.min_reads, classes=a.classes)
    except ValueError as e:
        sys.exit(f"ngsdose sinks: {e}")
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
    c.add_argument("--mp-margin", type=float, default=0.01,
                   help="how far (relative) a component has to clear the fitted noise edge to count towards ctrlPC_mp (default 0.01)")
    c.add_argument("--control-pcs", default="mp",
                   help="components of the control regions' residual depth to write as ctrlPC columns: a number, or 'mp' (default): "
                        "count the components above the Marchenko-Pastur edge of the noise bulk, record the count (ctrlPC_mp) "
                        "and write twice as many (at least 20) so that `pcsweep` can look beyond it")
    c.set_defaults(fn=cmd_cohort)

    t = sub.add_parser("trios", help="transmission reliability from a pedigree")
    t.add_argument("table")
    t.add_argument("-p", "--pedigree", required=True, help="the 1000 Genomes layout (FamilyID SampleID FatherID MotherID Sex Population ...), a PLINK "
                   "PED/FAM, or a child father mother [population] table; whitespace-separated, header optional")
    t.add_argument("-c", "--columns", nargs="+", required=True)
    t.add_argument("--perm", type=int, default=1000)
    t.add_argument("--compare-to", help="baseline column: paired bootstrap of the reliability difference of every other column against it")
    t.add_argument("--no-population-centring", action="store_true")
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

    j = sub.add_parser("adjust", help="residualise estimates on coverage PCs")
    j.add_argument("table")
    pc_source(j)
    j.add_argument("-c", "--columns", nargs="+", required=True)
    j.add_argument("--no-log", action="store_true")
    j.add_argument("-o", "--out", default="-")
    j.set_defaults(fn=cmd_adjust)

    w = sub.add_parser("pcsweep", help="what each further PC does: cross-validated error of the known truths, transmission of the classes",
                       description="Regress out 0, 1, 2, ... PCs and report, for every number: the cross-validated error of the "
                                   "known-truth columns (held-out autosomal = 2, chrX and chrY by sex, DJ = 10; --truth adds others) "
                                   "and, with a pedigree, the transmission reliability of the class columns and its paired "
                                   "difference from no adjustment. The evidence on which the default number of PCs can be overruled.")
    w.add_argument("table")
    pc_source(w)
    w.add_argument("-c", "--columns", nargs="+", default=[], help="class columns (estimates without a known truth)")
    w.add_argument("--truth", nargs="+", metavar="COLUMN=VALUE", help="further columns with a known value")
    w.add_argument("-p", "--pedigree")
    w.add_argument("--max-pc", type=int, default=0, help="sweep up to this many PCs (default: twice the default choice, at least 20)")
    w.add_argument("--folds", type=int, default=10)
    w.add_argument("--boot", type=int, default=300)
    w.add_argument("-o", "--out", default="-")
    w.set_defaults(fn=cmd_pcsweep)

    k = sub.add_parser("sinks", help="learn fetch-mode sink intervals from scan-mode counts")
    k.add_argument("counts", nargs="+")
    k.add_argument("-o", "--out", default="-")
    k.add_argument("--min-frac", type=float, default=1e-5, help="keep placement bins holding at least this fraction of a class's reads")
    k.add_argument("--evaluate", metavar="BED", help="instead of learning, report how much of each class an existing sinks BED captures")
    k.add_argument("--min-reads", type=int, default=25)
    k.add_argument("--pad", type=int, default=1000)
    k.add_argument("--classes", nargs="+", metavar="CLASS", help="also learn sinks for these compositional classes (TEL: the aligner concentrates telomeric "
                   "reads at the chromosome ends, so they are fetchable; the satellite families are dispersed and are not)")
    k.set_defaults(fn=cmd_sinks)

    s = sub.add_parser("selftest", help="simulation-based check of the estimator; needs no data")
    s.set_defaults(fn=cmd_selftest)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
