"""Transmission reliability: what fraction of an estimate's variance is real?

Array dosage is a physical DNA quantity inherited additively: a child's diploid total T_c is one
transmitted haploid complement from each parent, so E[T_c | parents] = (T_f + T_m) / 2 exactly,
with a segregation term s of mean zero. With X = T + e (measurement error of variance s_e^2,
independent between individuals) and rho the observed spousal correlation,

    midparent slope  b   = R (1 + rho_T) / (1 + R rho_T)      =>   R = b - rho (1 - b)
    single-parent slope  = (R + rho) / 2                       =>   R = 2 b_p - rho
    Mendelian variance   Var(X_c - X_mid) = s_T^2/2 + 1.5 s_e^2 =>   R = 1.5 - D / V

where R = s_T^2 / (s_T^2 + s_e^2) is the reliability. (The last line assumes random mating and
independent assortment; de novo or culture-induced change inflates D and makes it conservative.)

All three are inflated by error *shared within a family* (a trio libraried and sequenced
together). The spousal correlation is the direct test: spouses share no dosage by descent, so
after removing population means rho measures shared error (plus assortment). Technical
replicates are the clean estimate; trios give an upper bound.

Inheritance is additive on the natural scale, so the regressions are run on copies or Mb, not
on logs. Values are centred within population first, because a pooled regression across
populations with different means measures ancestry, not transmission.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Trio:
    child: str
    father: str
    mother: str
    population: str = ""


def load_pedigree(path) -> tuple[list[Trio], dict[str, str]]:
    """1000 Genomes style pedigree: FamilyID SampleID FatherID MotherID Sex Population ...
    (whitespace separated, header optional). Returns complete trios and sample -> population."""
    trios, pop = [], {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 4 or p[1] in ("SampleID", "IID", "sampleID"):
                continue
            pop[p[1]] = p[5] if len(p) > 5 else ""
            if p[2] not in ("0", "-9", "NA") and p[3] not in ("0", "-9", "NA"):
                trios.append(Trio(p[1], p[2], p[3], p[5] if len(p) > 5 else ""))
    return trios, pop


def _ols(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    xm, ym = x - x.mean(), y - y.mean()
    b = float((xm * ym).sum() / (xm ** 2).sum())
    res = ym - b * xm
    # heteroscedasticity-robust (HC1) standard error
    se = float(np.sqrt(((xm * res) ** 2).sum()) / (xm ** 2).sum() * np.sqrt(len(x) / max(len(x) - 2, 1)))
    return b, se


def centre_within(values: dict[str, float], group: dict[str, str], min_n: int = 5) -> dict[str, float]:
    """Subtract group (population) means estimated from all samples with a value."""
    by = {}
    for s, v in values.items():
        if np.isfinite(v):
            by.setdefault(group.get(s, ""), []).append(v)
    grand = np.mean([v for vs in by.values() for v in vs])
    mean = {g: (np.mean(vs) if len(vs) >= min_n else grand) for g, vs in by.items()}
    return {s: v - mean[group.get(s, "")] + grand for s, v in values.items() if np.isfinite(v)}


def _estimators(c, f, m) -> dict:
    mid = (f + m) / 2
    rho = float(np.corrcoef(f, m)[0, 1])
    b, se = _ols(mid, c)
    bf, sef = _ols(f, c)
    bm, sem = _ols(m, c)
    V = float(np.var(np.r_[f, m], ddof=1))
    D = float(np.var(c - mid, ddof=1))
    r = lambda x, y: float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else float("nan")
    return dict(spousal_r=rho, midparent_slope=b, midparent_slope_se=se, reliability_midparent=b - rho * (1 - b),
                r_midparent=r(mid, c), r_father=r(f, c), r_mother=r(m, c),
                father_slope=bf, father_slope_se=sef, mother_slope=bm, mother_slope_se=sem,
                reliability_single_parent=(bf + bm) - rho, mendel_D_over_V=D / V, reliability_mendel=1.5 - D / V,
                parent_sd=float(np.sqrt(V)), child_minus_midparent_sd=float(np.sqrt(D)))


def transmission(values: dict[str, float], trios: list[Trio], population: dict[str, str] | None = None,
                 n_perm: int = 1000, seed: int = 1, n_boot: int = 1000) -> dict:
    """Reliability estimates from complete trios. `values` are on the natural scale.

    Confidence intervals are percentile intervals from resampling families, which carries the
    uncertainty of the spousal correlation into the corrected reliabilities (a slope's own
    standard error does not)."""
    if population:
        values = centre_within(values, population)
    t = [x for x in trios if all(s in values and np.isfinite(values[s]) for s in (x.child, x.father, x.mother))]
    n = len(t)
    if n < 3:
        raise ValueError(f"need at least 3 complete trios, got {n}")
    c = np.array([values[x.child] for x in t])
    f = np.array([values[x.father] for x in t])
    m = np.array([values[x.mother] for x in t])
    mid = (f + m) / 2
    out = dict(n_trios=n, **_estimators(c, f, m), parent_mean=float(np.mean(np.r_[f, m])),
               child_minus_midparent_mean=float(np.mean(c - mid)))
    V = out["parent_sd"] ** 2
    out["error_cv"] = float(np.sqrt(max(0.0, 1 - out["reliability_midparent"]) * V) / out["parent_mean"])
    rng = np.random.default_rng(seed)
    if n_perm and n >= 10:
        null = np.array([_ols(mid, c[rng.permutation(n)])[0] for _ in range(n_perm)])
        out["perm_null_mean"], out["perm_null_sd"] = float(null.mean()), float(null.std())
        # one-sided: how often a slope at least as large arises when children are shuffled among the families
        out["perm_p"] = float((1 + np.sum(null >= out["midparent_slope"])) / (n_perm + 1))
    if n_boot and n >= 20:
        keys = ("reliability_midparent", "reliability_single_parent", "reliability_mendel", "spousal_r", "r_midparent")
        draws = {k: [] for k in keys}
        for _ in range(n_boot):
            i = rng.integers(0, n, n)
            e = _estimators(c[i], f[i], m[i])
            for k in keys:
                draws[k].append(e[k])
        for k in keys:
            out[k + "_ci95"] = (float(np.percentile(draws[k], 2.5)), float(np.percentile(draws[k], 97.5)))
    return out


def compare(values_a: dict[str, float], values_b: dict[str, float], trios: list[Trio],
            population: dict[str, str] | None = None, n_boot: int = 2000, seed: int = 1) -> dict:
    """Is estimator A more reliable than estimator B? Paired family bootstrap of R_A - R_B.

    Two estimators of the same quantity are strongly correlated, so their reliabilities can be
    told apart far more finely than two separate confidence intervals suggest."""
    if population:
        values_a, values_b = centre_within(values_a, population), centre_within(values_b, population)
    t = [x for x in trios if all(s in values_a and s in values_b for s in (x.child, x.father, x.mother))]
    n = len(t)
    if n < 20:
        raise ValueError(f"need at least 20 complete trios for a paired bootstrap, got {n}")
    arr = lambda v: tuple(np.array([v[getattr(x, who)] for x in t]) for who in ("child", "father", "mother"))
    (ca, fa, ma), (cb, fb, mb) = arr(values_a), arr(values_b)
    point = _estimators(ca, fa, ma)["reliability_midparent"] - _estimators(cb, fb, mb)["reliability_midparent"]
    rng = np.random.default_rng(seed)
    d = np.empty(n_boot)
    for j in range(n_boot):
        i = rng.integers(0, n, n)
        d[j] = _estimators(ca[i], fa[i], ma[i])["reliability_midparent"] - _estimators(cb[i], fb[i], mb[i])["reliability_midparent"]
    return dict(n_trios=n, delta=float(point), ci95=(float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))),
                p_a_better=float(np.mean(d > 0)))


def mendelian_z(values: dict[str, float], se: dict[str, float], trios: list[Trio]) -> list[dict]:
    """Per-trio departure of the child from the midparent, for outlier review. The segregation
    variance is unknown per family, so this is a screen, not a test."""
    rows = []
    for x in trios:
        if all(s in values for s in (x.child, x.father, x.mother)):
            d = values[x.child] - (values[x.father] + values[x.mother]) / 2
            e = np.sqrt(se.get(x.child, 0) ** 2 + (se.get(x.father, 0) ** 2 + se.get(x.mother, 0) ** 2) / 4)
            rows.append(dict(child=x.child, father=x.father, mother=x.mother, delta=float(d), meas_se=float(e),
                             within_parental_range=bool(d + (values[x.father] + values[x.mother]) / 2 <= values[x.father] + values[x.mother])))
    return rows
