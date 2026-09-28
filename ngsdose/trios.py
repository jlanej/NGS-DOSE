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

import re
import warnings
from dataclasses import dataclass, replace

import numpy as np


class PedigreeWarning(UserWarning):
    """A pedigree row that could not be read, or trios whose values could not be centred within population."""


@dataclass
class Trio:
    child: str
    father: str
    mother: str
    population: str = ""


MISSING_PARENT = ("0", "-9", "NA", ".", "")
MISSING_POPULATION = ("", ".", "NA", "na", "-9")
_CHILD_NAMES = ("sampleid", "iid", "sample", "child", "kid", "proband", "individualid", "samplename")
_FATHER_NAMES = ("fatherid", "pat", "father", "dad", "paternalid")
_MOTHER_NAMES = ("motherid", "mat", "mother", "mom", "maternalid")
_POP_NAMES = ("population", "pop")
_norm = lambda x: re.sub(r"[\s_.-]", "", x.lower())      # 'Individual ID', 'paternal_id' -> 'individualid', 'paternalid'


def pedigree_layout(first: list[str]) -> tuple[int, int, int, int | None, bool]:
    """(child, father, mother, population column, header present) for the first line of a pedigree file.

    Three layouts are read: the 1000 Genomes one (FamilyID SampleID FatherID MotherID Sex Population ...),
    PLINK PED/FAM (FID IID PAT MAT SEX PHENOTYPE: the sixth column is a phenotype code, not a population),
    and a trios table (child father mother [population]). A header names the columns (case, spaces and
    underscores aside); without one the layout is told from the number of columns, and a population is read
    only from the fourth column of a four-column table: the sixth column of a headerless file is never taken
    for one, since in PLINK it holds a phenotype (-9, 0, 1, 2, NA or a quantitative value)."""
    low = [_norm(x) for x in first]
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
    return 1, 2, 3, None, False                            # PLINK, or the 1000 Genomes layout without its header


def _split(line: str, by_tab: bool) -> list[str]:
    return [x.strip() for x in line.rstrip("\r\n").split("\t")] if by_tab else line.split()


def _header(line: str, groups) -> tuple[list[str], bool] | None:
    """(fields, split on tabs) when `line` names a column of each of `groups`: split on whitespace, or on tabs
    when the names hold spaces (the 1000 Genomes 20130606_g1k.ped: 'Family ID<TAB>Individual ID<TAB>...')."""
    for by_tab in (False, True):
        if by_tab and "\t" not in line:
            break
        fields = _split(line, by_tab)
        low = [_norm(x) for x in fields]
        if all(any(x in names for x in low) for names in groups):
            return fields, by_tab
    return None


def load_population(path) -> dict[str, str]:
    """sample -> population, or any grouping to centre within (clusters on ancestry PCs for a biobank), from a
    file of sample and group. A header naming a sample column and a population column picks them from a wider
    table (a 1000 Genomes pedigree, a .psam), and a header naming a sample column alone is an error; a file
    without a header gives its first two columns. Other lines starting
    with '#' are comments, and a label of NA, . or -9 is missing."""
    with open(path) as fh:
        raw = [line for line in fh if line.strip()]
    names = _POP_NAMES + ("group", "cluster", "ancestry", "superpop", "superpopulation")    # in order of preference
    h = _header(raw[0].lstrip("#"), (_CHILD_NAMES, names)) if raw else None
    if raw and h is None and _header(raw[0].lstrip("#"), (_CHILD_NAMES,)) is not None:
        # a header with a sample column but no group column (a .psam's '#IID SEX ...'): its second column is not a group
        raise ValueError(f"{path}: the header names a sample column but no population/group column ({', '.join(names)})")
    lines = [line for line in raw[1 if h else 0:] if not line.startswith("#")]
    if h:
        low = [_norm(x) for x in h[0]]
        s, g = next(i for i, x in enumerate(low) if x in _CHILD_NAMES), next(low.index(n) for n in names if n in low)
        by_tab = h[1]
    else:
        s, g, by_tab = 0, 1, False
    out = {}
    for line in lines:
        p = _split(line, by_tab)
        if len(p) > max(s, g) and p[g] not in MISSING_POPULATION:
            out[p[s]] = p[g]
    return out


def load_pedigree(path, population=None) -> tuple[list[Trio], dict[str, str]]:
    """Complete trios and sample -> population from a pedigree file in any layout `pedigree_layout` reads
    (whitespace-separated, or tab-separated under a header whose names hold spaces). A parent given as 0, -9,
    NA or . is absent; a population column is optional, and the dict holds only the samples with a label, so
    it is empty for a file without one. `population` (a dict, or a file for `load_population`) gives labels
    that replace the pedigree's for the samples it names. Rows too short for the layout are skipped with a
    PedigreeWarning."""
    trios, pop = [], {}
    with open(path) as fh:
        raw = [line for line in fh if line.strip()]
    body = [line for line in raw if not line.startswith("#")]
    head, by_tab = None, False
    groups = (_CHILD_NAMES, _FATHER_NAMES, _MOTHER_NAMES)
    # a first line that starts with '#' is the header when it names the child, father and mother columns and
    # has the data's column count (PLINK's "#FID IID PAT MAT ...", a "#kid dad mom ..." table); every other
    # line starting with '#' is a comment
    if raw and raw[0].startswith("#") and body:
        h = _header(raw[0].lstrip("#"), groups)
        if h and len(h[0]) == len(_split(body[0], h[1])):
            head, by_tab = h
    if head is None and body:
        h = _header(body[0], groups)
        if h:
            (head, by_tab), body = h, body[1:]
    rows = [_split(line, by_tab) for line in body]
    if head is None and not rows:
        return trios, pop
    if head is not None:
        c, f, m, g, _ = pedigree_layout(head)
    else:
        pedigree_layout(rows[0])                            # a header that names a child but no parents is an error
        # without a header the layout is guessed from the most common column count, not from the first row alone
        counts = {}
        for p in rows:
            counts[len(p)] = counts.get(len(p), 0) + 1
        n = max(counts, key=lambda k: (counts[k], -k))
        c, f, m, g, _ = pedigree_layout(next(p for p in rows if len(p) == n))
        if len(counts) > 1:
            warnings.warn(f"{path}: no header, so the layout was taken from the {counts[n]} rows with {n} columns; "
                          f"{len(rows) - counts[n]} rows have another column count", PedigreeWarning, stacklevel=2)
    short = 0
    for p in rows:
        if len(p) <= max(c, f, m):
            short += 1
            continue
        label = p[g] if g is not None and len(p) > g and p[g] not in MISSING_POPULATION else ""
        if label:
            pop[p[c]] = label
        if p[f] not in MISSING_PARENT and p[m] not in MISSING_PARENT:
            trios.append(Trio(p[c], p[f], p[m], label))
    if short:
        warnings.warn(f"{path}: {short} rows with {max(c, f, m)} columns or fewer skipped", PedigreeWarning, stacklevel=2)
    for t in trios:                                        # a trios table names the parents only on the child's row
        for parent in (t.father, t.mother):
            if t.population and not pop.get(parent):
                pop[parent] = t.population
    if population is not None:
        extra = population if isinstance(population, dict) else load_population(population)
        pop.update({s: g for s, g in extra.items() if g not in MISSING_POPULATION})
        trios = [replace(t, population=pop.get(t.child, "")) for t in trios]
    return trios, pop


def families(trios: list[Trio]) -> list[list[int]]:
    """Indices of `trios` grouped into families, in order of first appearance: trios that share anyone
    (siblings, half-siblings, a child who is a parent in another trio) are one family."""
    up = list(range(len(trios)))

    def root(i):
        while up[i] != i:
            up[i] = up[up[i]]
            i = up[i]
        return i
    first = {}
    for i, t in enumerate(trios):
        for s in (t.child, t.father, t.mother):
            a, b = root(i), root(first.setdefault(s, i))
            if a != b:
                up[max(a, b)] = min(a, b)
    out = {}
    for i in range(len(trios)):
        out.setdefault(root(i), []).append(i)
    return list(out.values())


def _resample(fam: list[np.ndarray], rng) -> np.ndarray:
    """Trio indices of one bootstrap draw of whole families (with one trio per family, a draw of trios)."""
    k = rng.integers(0, len(fam), len(fam))
    return np.concatenate([fam[j] for j in k])


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
    `min_n`, and of everyone where the sex has fewer too), so that a difference between men and women is not
    read as transmission."""
    by, by_sex = {}, {}
    for s, v in values.items():
        if np.isfinite(v):
            by.setdefault((population.get(s, ""), sex.get(s, "")), []).append(v)
            by_sex.setdefault(sex.get(s, ""), []).append(v)
    if not by_sex:
        return {}
    mean = {g: float(np.mean(vs)) for g, vs in by.items() if len(vs) >= min_n}
    sex_mean = {g: float(np.mean(vs)) for g, vs in by_sex.items() if len(vs) >= min_n}
    grand = float(np.mean([v for vs in by_sex.values() for v in vs]))
    return {s: v - mean.get((population.get(s, ""), sex.get(s, "")), sex_mean.get(sex.get(s, ""), grand)) for s, v in values.items() if np.isfinite(v)}


def sex_code(x) -> str:
    """'M', 'F' or '' (unknown) from M/F, male/female or PLINK's 1/2."""
    return {"m": "M", "male": "M", "1": "M", "f": "F", "female": "F", "2": "F"}.get(str(x).strip().lower(), "")


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
    degrees of freedom instead, which assumes normal values and is far too small for skewed ones.

    Sex is read as M/F, male/female or PLINK's 1/2; a parent without one takes its role's. A ValueError is
    raised when trios are given but no child has a sex."""
    sex = {s: sex_code(x) for s, x in sex.items()}
    for t in trios:                                        # the role fixes a parent's sex where the pedigree gives none
        for s, code in ((t.father, "M"), (t.mother, "F")):
            if not sex.get(s):
                sex[s] = code
    if trios and not any(sex.get(t.child) for t in trios):
        raise ValueError("no child has a sex (M/F or PLINK 1/2): transmission by sex needs the children's")
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


_UNCENTRED = ("values not centred within population: the trios' samples are not labelled with at least two populations "
              "(of 5 or more samples), so a difference between populations can pass for transmission. Give the labels "
              "(a population column, or a population file), or turn centring off for a cohort of one population")


def _populations(values: dict[str, float], trios: list[Trio], population: dict[str, str] | None, min_n: int = 5) -> tuple[int, bool]:
    """(populations that centring separates: labels with `min_n` or more samples with a value, whether the lack
    of them deserves a warning). No warning when centring was turned off (None) or when every trio's samples
    carry one and the same label: one population has nothing to centre."""
    if population is None:
        return 0, False
    n = {}
    for s, v in values.items():
        if np.isfinite(v) and population.get(s):
            n[population[s]] = n.get(population[s], 0) + 1
    k = sum(x >= min_n for x in n.values())
    labels = {population.get(s, "") for x in trios for s in (x.child, x.father, x.mother) if s in values}
    return k, k < 2 and not (len(labels) == 1 and "" not in labels)


def transmission(values: dict[str, float], trios: list[Trio], population: dict[str, str] | None = None,
                 n_perm: int = 1000, seed: int = 1, n_boot: int = 1000) -> dict:
    """Reliability estimates from complete trios. `values` are on the natural scale.

    Confidence intervals are percentile intervals from resampling families: trios that share anyone
    (`families`) are drawn together, so siblings do not count as independent. This carries the
    uncertainty of the spousal correlation into the corrected reliabilities (a slope's own
    standard error does not). `population_centred` says whether values were centred within two or
    more populations. When they were not although `population` was given (a pedigree without
    labels), a PedigreeWarning is issued; `warning` in the result notes that, and a spousal
    correlation above 0.2 without centring."""
    n_pop, warn = _populations(values, trios, population)
    if population:
        values = centre_within(values, population)
    t = [x for x in trios if all(s in values and np.isfinite(values[s]) for s in (x.child, x.father, x.mother))]
    n = len(t)
    if n < 3:
        raise ValueError(f"need at least 3 complete trios, got {n}")
    fam = [np.array(g) for g in families(t)]
    c = np.array([values[x.child] for x in t])
    f = np.array([values[x.father] for x in t])
    m = np.array([values[x.mother] for x in t])
    mid = (f + m) / 2
    out = dict(n_trios=n, n_families=len(fam), population_centred=n_pop >= 2, n_populations=n_pop, **_estimators(c, f, m),
               parent_mean=float(np.mean(np.r_[f, m])), child_minus_midparent_mean=float(np.mean(c - mid)))
    notes = [_UNCENTRED] if warn else []
    if warn:
        warnings.warn(_UNCENTRED, PedigreeWarning, stacklevel=2)
    if n_pop < 2 and out["spousal_r"] > 0.2:
        notes.append("spousal correlation above 0.2 without population centring: population structure, or error shared "
                     "within families, inflates R")
    if notes:
        out["warning"] = "; ".join(notes)
    V = out["parent_sd"] ** 2
    out["error_cv"] = (float(np.sqrt(max(0.0, 1 - out["reliability_midparent"]) * V) / out["parent_mean"])
                       if np.isfinite(V) and np.isfinite(out["reliability_midparent"]) and out["parent_mean"] > 0 else float("nan"))
    rng = np.random.default_rng(seed)
    if n_perm and n >= 10:
        if np.isfinite(out["midparent_slope"]):
            null = np.array([_ols(mid, c[rng.permutation(n)])[0] for _ in range(n_perm)])
            out["perm_null_mean"], out["perm_null_sd"] = float(null.mean()), float(null.std())
            # one-sided: how often a slope at least as large arises when children are shuffled among the families
            out["perm_p"] = float((1 + np.sum(null >= out["midparent_slope"])) / (n_perm + 1))
        else:                                              # no slope (a constant column): no test, not a small p
            out["perm_null_mean"] = out["perm_null_sd"] = out["perm_p"] = float("nan")
    if n_boot and n >= 20:
        keys = ("reliability_midparent", "reliability_single_parent", "reliability_mendel", "spousal_r", "r_midparent",
                "reliability_rescaled", "sd_ratio", "mean_ratio")
        draws = {k: [] for k in keys}
        for _ in range(n_boot):
            i = _resample(fam, rng)
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
    told apart far more finely than two separate confidence intervals suggest. Only trios with finite
    values under both are used, and whole families are resampled, as in `transmission`.
    `p_a_better` is NaN when the difference or any bootstrap draw of it is (a constant column)."""
    n_pop, warn = _populations(values_a, trios, population)
    if warn:
        warnings.warn(_UNCENTRED, PedigreeWarning, stacklevel=2)
    if population:
        values_a, values_b = centre_within(values_a, population), centre_within(values_b, population)
    t = [x for x in trios if all(s in v and np.isfinite(v[s]) for v in (values_a, values_b) for s in (x.child, x.father, x.mother))]
    n = len(t)
    if n < 20:
        raise ValueError(f"need at least 20 complete trios for a paired bootstrap, got {n}")
    fam = [np.array(g) for g in families(t)]
    arr = lambda v: tuple(np.array([v[getattr(x, who)] for x in t]) for who in ("child", "father", "mother"))
    (ca, fa, ma), (cb, fb, mb) = arr(values_a), arr(values_b)
    point = _estimators(ca, fa, ma)["reliability_midparent"] - _estimators(cb, fb, mb)["reliability_midparent"]
    rng = np.random.default_rng(seed)
    d = np.empty(n_boot)
    for j in range(n_boot):
        i = _resample(fam, rng)
        d[j] = _estimators(ca[i], fa[i], ma[i])["reliability_midparent"] - _estimators(cb[i], fb[i], mb[i])["reliability_midparent"]
    p = float(np.mean(d > 0)) if np.isfinite(point) and np.isfinite(d).all() else float("nan")
    return dict(n_trios=n, n_families=len(fam), population_centred=n_pop >= 2, delta=float(point),
                ci95=(float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))), p_a_better=p)
