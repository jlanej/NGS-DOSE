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

The last line also assumes the children are measured on their parents' scale. In general, with s the
children's standard deviation over the parents' and equal parental variances, it is exactly

    R_mendel = R + 1 + rho/2 - s^2

so it adds to the midparent slope only whether the children vary more or less than their parents:
more when variation arises new in them, and either way when a batch that measured only the children
reads on another scale. `sd_ratio` reports s itself. A batch that multiplies the children's values by
k multiplies the slope, and R with it, by k; `reliability_rescaled` is R with the children first
rescaled to their parents' spread (b / s for b), which such a factor does not move and new variation in
the children lowers. The two bracket the reliability when generation and batch go together, as in
the 1000 Genomes 30x release, which sequenced the children after their parents.

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


MISSING_PARENT = ("0", "-9", "NA", ".", "")
_CHILD_NAMES = ("sampleid", "iid", "sample", "child", "kid", "proband")
_FATHER_NAMES = ("fatherid", "pat", "father", "dad", "paternal_id")
_MOTHER_NAMES = ("motherid", "mat", "mother", "mom", "maternal_id")
_POP_NAMES = ("population", "pop")


def pedigree_layout(first: list[str]) -> tuple[int, int, int, int | None, bool]:
    """(child, father, mother, population column, header present) for the first line of a pedigree file.

    Three layouts are read: the 1000 Genomes one (FamilyID SampleID FatherID MotherID Sex Population ...),
    PLINK PED/FAM (FID IID PAT MAT SEX PHENOTYPE: the sixth column is a phenotype code, not a population),
    and a trios table (child father mother [population]). A header names the columns; without one the
    layout is told from the number of columns and whether the sixth is a phenotype code."""
    low = [x.lower() for x in first]
    if any(x in _CHILD_NAMES for x in low):
        col = lambda names: next((i for i, x in enumerate(low) if x in names), None)
        c, f, m = col(_CHILD_NAMES), col(_FATHER_NAMES), col(_MOTHER_NAMES)
        if f is None or m is None:
            raise ValueError(f"pedigree header names no father/mother column: {' '.join(first)}")
        return c, f, m, col(_POP_NAMES), True
    if len(first) == 3:
        return 0, 1, 2, None, False
    if len(first) == 4 and first[3] not in ("1", "2", "0"):
        return 0, 1, 2, 3, False
    if len(first) >= 6 and first[5] in ("-9", "0", "1", "2"):
        return 1, 2, 3, None, False                        # PLINK: phenotype, not population
    return 1, 2, 3, (5 if len(first) > 5 else None), False


def load_pedigree(path) -> tuple[list[Trio], dict[str, str]]:
    """Complete trios and sample -> population from a pedigree file in any layout `pedigree_layout` reads
    (whitespace-separated). A parent given as 0, -9, NA or . is absent; a population column is optional."""
    trios, pop = [], {}
    with open(path) as fh:
        lines = [line.split() for line in fh if line.strip() and not line.startswith("#")]
    if not lines:
        return trios, pop
    c, f, m, g, header = pedigree_layout(lines[0])
    for p in lines[1:] if header else lines:
        if len(p) <= max(c, f, m):
            continue
        pop[p[c]] = p[g] if g is not None and len(p) > g else ""
        if p[f] not in MISSING_PARENT and p[m] not in MISSING_PARENT:
            trios.append(Trio(p[c], p[f], p[m], pop[p[c]]))
    for t in trios:                                        # a trios table names the parents only on the child's row
        for parent in (t.father, t.mother):
            if t.population and not pop.get(parent):
                pop[parent] = t.population
    return trios, pop


def _ols(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    xm, ym = x - x.mean(), y - y.mean()
    sxx = float((xm ** 2).sum())
    if not sxx > 0:                                        # a constant regressor: no slope to estimate
        return float("nan"), float("nan")
    b = float((xm * ym).sum() / sxx)
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
    if not by:                                             # nothing to centre (a column absent from these rows)
        return {}
    grand = np.mean([v for vs in by.values() for v in vs])
    mean = {g: (np.mean(vs) if len(vs) >= min_n else grand) for g, vs in by.items()}
    return {s: v - mean[group.get(s, "")] + grand for s, v in values.items() if np.isfinite(v)}


def centre_within_sex(values: dict[str, float], population: dict[str, str], sex: dict[str, str], min_n: int = 5) -> dict[str, float]:
    """Subtract the mean of each sample's population and sex (of its sex alone where that group has fewer than
    `min_n`), so that a difference between men and women is not read as transmission."""
    by, by_sex = {}, {}
    for s, v in values.items():
        if np.isfinite(v):
            by.setdefault((population.get(s, ""), sex.get(s, "")), []).append(v)
            by_sex.setdefault(sex.get(s, ""), []).append(v)
    mean = {g: float(np.mean(vs)) for g, vs in by.items() if len(vs) >= min_n}
    sex_mean = {g: float(np.mean(vs)) for g, vs in by_sex.items()}
    return {s: v - mean.get((population.get(s, ""), sex.get(s, "")), sex_mean[sex.get(s, "")]) for s, v in values.items() if np.isfinite(v)}


# parent -> child pairings by sex: (key, the parent, the child's sex)
PAIRS = (("father_son", "father", "M"), ("father_daughter", "father", "F"), ("mother_son", "mother", "M"), ("mother_daughter", "mother", "F"))


def by_sex(values: dict[str, float], trios: list[Trio], population: dict[str, str] | None, sex: dict[str, str], min_pairs: int = 10,
           n_perm: int = 1000, seed: int = 1) -> dict:
    """Parent-to-child transmission split by the sex of parent and child.

    An autosomal quantity passes half of each parent's deviation to every child (single-parent slope
    (R + rho) / 2, as above); a Y-linked one all of the father's to his sons and none to his daughters; an
    X-linked one all of the father's to his daughters, none to his sons, and half of the mother's to every
    child. Values are centred within population and sex. Per pairing: pairs, the slope of child on parent
    with its robust standard error, and Pearson r with a Fisher 95% interval. The correlation is the one to
    compare between pairings: a slope also carries the ratio of the child's spread to the parent's, which
    differs between the sexes for a sex-linked quantity. `father_contrast` compares the father-son with the
    father-daughter correlation (independent groups, Fisher z): near 0 for an autosomal quantity, large and
    positive for a Y-linked one, large and negative for an X-linked one; `father_slope_contrast` does the
    same with the slopes. `heterogeneity` asks whether the four pairings' correlations differ at all, as an
    autosomal quantity's should not: Cochran's Q on Fisher's z (each pairing weighted by n - 3) over the trios
    with all three values, referred to its distribution when the children's sexes are shuffled among the
    families and each family's parents swap roles at random (`p`). `p_normal` refers Q to chi-square with 3
    degrees of freedom instead, which assumes normal values and is far too small for skewed ones."""
    v = centre_within_sex(values, population or {}, sex)
    out = {}
    for key, who, child_sex in PAIRS:
        xy = [(v[getattr(t, who)], v[t.child]) for t in trios if sex.get(t.child) == child_sex and getattr(t, who) in v and t.child in v]
        if len(xy) < min_pairs:
            continue
        x, y = np.array([a for a, _ in xy]), np.array([b for _, b in xy])
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        b, se = _ols(x, y)
        r = float(np.corrcoef(x, y)[0, 1])
        z, h = float(np.arctanh(np.clip(r, -0.999999, 0.999999))), 1.96 / np.sqrt(max(len(xy) - 3, 1))
        out[key] = dict(n=len(xy), slope=b, slope_se=se, r=r, r_lo=float(np.tanh(z - h)), r_hi=float(np.tanh(z + h)))
    if "father_son" in out and "father_daughter" in out:
        a, d = out["father_son"], out["father_daughter"]
        fz = lambda r: float(np.arctanh(np.clip(r, -0.999999, 0.999999)))
        se_z = float(np.sqrt(1 / max(a["n"] - 3, 1) + 1 / max(d["n"] - 3, 1)))
        out["father_contrast"] = dict(diff=a["r"] - d["r"], z=(fz(a["r"]) - fz(d["r"])) / se_z)
        diff, se = a["slope"] - d["slope"], float(np.hypot(a["slope_se"], d["slope_se"]))
        out["father_slope_contrast"] = dict(diff=diff, se=se, z=diff / se if se > 0 else float("nan"))
    full = [t for t in trios if sex.get(t.child) in ("M", "F") and all(s in v for s in (t.child, t.father, t.mother))]
    if n_perm and all(k in out for k, _, _ in PAIRS) and min(sum(sex[t.child] == "M" for t in full), sum(sex[t.child] == "F" for t in full)) >= min_pairs:
        f, m, c = (np.array([v[getattr(t, who)] for t in full]) for who in ("father", "mother", "child"))
        son = np.array([sex[t.child] == "M" for t in full])
        q = _pairing_q(f, m, c, son)
        rng = np.random.default_rng(seed)
        null = np.empty(n_perm)
        for j in range(n_perm):
            swap = rng.random(len(full)) < 0.5
            null[j] = _pairing_q(np.where(swap, m, f), np.where(swap, f, m), c, rng.permutation(son))
        out["heterogeneity"] = dict(q=q, p=float((1 + np.sum(null >= q)) / (n_perm + 1)), p_normal=_chi2_sf3(q), n_trios=len(full), n_perm=n_perm)
    return out


def _pairing_q(f, m, c, son) -> float:
    """Cochran's Q over the four parent-child pairings' Fisher z, each weighted by its pairs less 3."""
    zw = []
    for parent in (f, m):
        for sel in (son, ~son):
            x, y = parent[sel], c[sel]
            r = np.corrcoef(x, y)[0, 1] if np.std(x) > 0 and np.std(y) > 0 else 0.0
            zw.append((np.arctanh(np.clip(r, -0.999999, 0.999999)), len(x) - 3))
    zbar = sum(z * w for z, w in zw) / sum(w for _, w in zw)
    return float(sum(w * (z - zbar) ** 2 for z, w in zw))


def _chi2_sf3(x: float) -> float:
    """Upper tail of chi-square with 3 degrees of freedom, in closed form."""
    from math import erfc, exp, pi, sqrt
    return float(erfc(sqrt(x / 2)) + sqrt(2 * x / pi) * exp(-x / 2))


def _estimators(c, f, m) -> dict:
    mid = (f + m) / 2
    rho = float(np.corrcoef(f, m)[0, 1]) if np.std(f) > 0 and np.std(m) > 0 else float("nan")
    b, se = _ols(mid, c)
    bf, sef = _ols(f, c)
    bm, sem = _ols(m, c)
    V = float(np.var(np.r_[f, m], ddof=1))
    D = float(np.var(c - mid, ddof=1))
    if not V > 0:                                          # constant among the parents (chrEBV in blood): nothing to inherit or to err in
        V = float("nan")
    r = lambda x, y: float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else float("nan")
    s = float(np.std(c, ddof=1) / np.sqrt(V)) if V > 0 else float("nan")
    pm = float(np.mean(np.r_[f, m]))
    return dict(spousal_r=rho, midparent_slope=b, midparent_slope_se=se, reliability_midparent=b - rho * (1 - b),
                r_midparent=r(mid, c), r_father=r(f, c), r_mother=r(m, c),
                father_slope=bf, father_slope_se=sef, mother_slope=bm, mother_slope_se=sem,
                reliability_single_parent=(bf + bm) - rho, mendel_D_over_V=D / V, reliability_mendel=1.5 - D / V,
                parent_sd=float(np.sqrt(V)), child_minus_midparent_sd=float(np.sqrt(D)),
                sd_ratio=s, reliability_rescaled=(b / s) * (1 + rho) - rho if s > 0 else float("nan"),
                mean_ratio=float(np.mean(c)) / pm if pm > 0 else float("nan"))


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
    out["error_cv"] = (float(np.sqrt(max(0.0, 1 - out["reliability_midparent"]) * V) / out["parent_mean"])
                       if np.isfinite(V) and np.isfinite(out["reliability_midparent"]) and out["parent_mean"] > 0 else float("nan"))
    rng = np.random.default_rng(seed)
    if n_perm and n >= 10:
        null = np.array([_ols(mid, c[rng.permutation(n)])[0] for _ in range(n_perm)])
        out["perm_null_mean"], out["perm_null_sd"] = float(null.mean()), float(null.std())
        # one-sided: how often a slope at least as large arises when children are shuffled among the families
        out["perm_p"] = float((1 + np.sum(null >= out["midparent_slope"])) / (n_perm + 1))
    if n_boot and n >= 20:
        keys = ("reliability_midparent", "reliability_single_parent", "reliability_mendel", "spousal_r", "r_midparent",
                "reliability_rescaled", "sd_ratio", "mean_ratio")
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
