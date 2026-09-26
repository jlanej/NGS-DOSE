"""Trio analysis on the pedigrees and families that biobanks have: a pedigree without population labels says so
instead of pooling populations silently, and labels can be supplied from a file; the old 1000 Genomes .ped
header and PLINK 2 phenotype codes are read correctly; siblings and three-generation families are resampled
as families; and `compare` does not report a probability next to a NaN difference."""
import warnings

import numpy as np
import pytest

from ngsdose.trios import PedigreeWarning, Trio, compare, families, load_pedigree, load_population, transmission


def _cohort(n_pop=5, per_pop=60, seed=5):
    """Populations with different mean dosage, true reliability 0.6 within each."""
    rng = np.random.default_rng(seed)
    rows, vals = [], {}
    for p in range(n_pop):
        mu = 250 + 40 * p
        for i in range(per_pop):
            c, f, m = f"C{p}_{i}", f"F{p}_{i}", f"M{p}_{i}"
            hf, hm = rng.normal(mu / 2, 20, 2), rng.normal(mu / 2, 20, 2)
            true = {f: hf.sum(), m: hm.sum(), c: rng.choice(hf) + rng.choice(hm)}
            for s, t in true.items():
                vals[s] = t + rng.normal(0, 23)
            rows.append((c, f, m, f"POP{p}"))
    return rows, vals


def _write_1000g(path, rows):
    path.write_text("FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
                    + "".join(f"FAM{i} {c} {f} {m} 1 {p} SUP\nFAM{i} {f} 0 0 1 {p} SUP\nFAM{i} {m} 0 0 2 {p} SUP\n" for i, (c, f, m, p) in enumerate(rows)))


def test_a_pedigree_without_populations_warns_and_a_population_file_restores_the_centred_estimate(tmp_path):
    rows, vals = _cohort()
    g1k = tmp_path / "g1k.txt"
    _write_1000g(g1k, rows)
    t_ok, pop_ok = load_pedigree(g1k)
    with warnings.catch_warnings():
        warnings.simplefilter("error", PedigreeWarning)
        centred = transmission(vals, t_ok, pop_ok, n_perm=0, n_boot=0)
    assert centred["population_centred"] and centred["n_populations"] == 5 and "warning" not in centred
    # the 1000 Genomes 30x release's pedigree_info layout, and a PLINK .fam: no population column
    info = tmp_path / "pedigree_info.txt"
    info.write_text("sampleID fatherID motherID sex\n" + "".join(f"{c} {f} {m} 1\n{f} 0 0 1\n{m} 0 0 2\n" for c, f, m, _ in rows))
    fam = tmp_path / "cohort.fam"
    fam.write_text("".join(f"F{i} {c} {f} {m} 1 -9\nF{i} {f} 0 0 1 -9\nF{i} {m} 0 0 2 -9\n" for i, (c, f, m, _) in enumerate(rows)))
    for path in (info, fam):
        t, pop = load_pedigree(path)
        assert pop == {} and len(t) == len(rows)
        with pytest.warns(PedigreeWarning, match="not centred within population"):
            pooled = transmission(vals, t, pop, n_perm=0, n_boot=0)
        assert not pooled["population_centred"] and "not centred" in pooled["warning"] and "spousal" in pooled["warning"]
        assert pooled["reliability_midparent"] > centred["reliability_midparent"] + 0.1        # ancestry read as transmission
        # centring turned off on purpose: no warning, but the result says so
        with warnings.catch_warnings():
            warnings.simplefilter("error", PedigreeWarning)
            off = transmission(vals, t, None, n_perm=0, n_boot=0)
        assert not off["population_centred"] and off["reliability_midparent"] == pooled["reliability_midparent"]
    # labels from a file: two columns, or the columns a header names in a wider table (here the 1000G pedigree itself)
    labels = tmp_path / "clusters.tsv"
    labels.write_text("".join(f"{s}\t{p}\n" for c, f, m, p in rows for s in (c, f, m)))
    for source in (labels, g1k, {s: p for c, f, m, p in rows for s in (c, f, m)}):
        t, pop = load_pedigree(info, population=source)
        assert pop == pop_ok and [x.population for x in t] == [p for *_, p in rows]
        again = transmission(vals, t, pop, n_perm=0, n_boot=0)
        assert again["reliability_midparent"] == centred["reliability_midparent"] and again["population_centred"]
    assert load_population(g1k) == pop_ok
    psam = tmp_path / "cohort.psam"
    psam.write_text("#IID\tSEX\tCluster\n# a comment\n" + "".join(f"{s}\t1\t{p}\n" for c, f, m, p in rows for s in (c, f, m)) + "X1\t1\tNA\n")
    assert load_population(psam) == pop_ok


def test_one_labelled_population_needs_no_centring_and_gives_no_warning():
    rows, vals = _cohort(n_pop=1, per_pop=80)
    t = [Trio(c, f, m, p) for c, f, m, p in rows]
    pop = {s: p for c, f, m, p in rows for s in (c, f, m)}
    with warnings.catch_warnings():
        warnings.simplefilter("error", PedigreeWarning)
        r = transmission(vals, t, pop, n_perm=0, n_boot=0)
    assert not r["population_centred"] and r["n_populations"] == 1


def test_the_old_1000g_ped_header_and_plink2_phenotypes(tmp_path):
    rows, _ = _cohort(n_pop=2, per_pop=5)
    old = tmp_path / "20130606_g1k.ped"                     # tab-separated, names with spaces
    old.write_text("Family ID\tIndividual ID\tPaternal ID\tMaternal ID\tGender\tPhenotype\tPopulation\tRelationship\tSiblings\tSecond Order\tThird Order\tOther Comments\n"
                   + "".join(f"FAM{i}\t{c}\t{f}\t{m}\t1\t0\t{p}\tchild\t0\t0\t0\t0\n{'FAM%d' % i}\t{f}\t0\t0\t1\t0\t{p}\tfather\t0\t0\t0\t0\n"
                             f"FAM{i}\t{m}\t0\t0\t2\t0\t{p}\tmother\t0\t0\t0\t0\n" for i, (c, f, m, p) in enumerate(rows)))
    t, pop = load_pedigree(old)
    assert t == [Trio(c, f, m, p) for c, f, m, p in rows]                   # no Trio('ID', 'Individual', 'ID') from the header
    assert set(pop.values()) == {"POP0", "POP1"}                           # the Population column, not Phenotype's 0
    # PLINK 2 writes NA for a missing phenotype; a quantitative phenotype is no population either
    for pheno in ("NA", "na", "1.37", "2"):
        fam = tmp_path / f"plink2_{pheno}.fam"
        fam.write_text("".join(f"F{i} {c} {f} {m} 1 {pheno}\nF{i} {f} 0 0 1 1\nF{i} {m} 0 0 2 2\n" for i, (c, f, m, _) in enumerate(rows)))
        t, pop = load_pedigree(fam)
        assert [(x.child, x.father, x.mother) for x in t] == [(c, f, m) for c, f, m, _ in rows] and pop == {}


def test_the_layout_follows_the_most_common_row_and_short_rows_are_reported(tmp_path):
    rows, _ = _cohort(n_pop=1, per_pop=6)
    fam = tmp_path / "short_first.fam"                      # one stray three-column row first, then a PLINK .fam
    fam.write_text("X Y Z\n" + "".join(f"F{i} {c} {f} {m} 1 -9\nF{i} {f} 0 0 1 -9\nF{i} {m} 0 0 2 -9\n" for i, (c, f, m, _) in enumerate(rows)))
    with pytest.warns(PedigreeWarning) as rec:
        t, _ = load_pedigree(fam)
    assert [(x.child, x.father, x.mother) for x in t] == [(c, f, m) for c, f, m, _ in rows]
    text = " ".join(str(w.message) for w in rec)
    assert "1 rows with 3 columns or fewer skipped" in text and "another column count" in text


def test_families_join_siblings_and_generations():
    t = [Trio("c1", "f", "m"), Trio("x", "a", "b"), Trio("c2", "f", "m"), Trio("g", "h", "c1"), Trio("y", "p", "q")]
    assert families(t) == [[0, 2, 3], [1], [4]]
    assert families([Trio(f"c{i}", f"f{i}", f"m{i}") for i in range(4)]) == [[0], [1], [2], [3]]


def _sibships(rng, n_fam, kids, R=0.6):
    """`n_fam` parent pairs with `kids` children each; true values of variance 1, error to reliability R."""
    se = np.sqrt((1 - R) / R)
    vals, trios = {}, []
    for i in range(n_fam):
        hf, hm = rng.normal(0, np.sqrt(0.5), 2), rng.normal(0, np.sqrt(0.5), 2)
        vals[f"f{i}"], vals[f"m{i}"] = hf.sum() + rng.normal(0, se), hm.sum() + rng.normal(0, se)
        for k in range(kids):
            vals[f"c{i}_{k}"] = rng.choice(hf) + rng.choice(hm) + rng.normal(0, se)
            trios.append(Trio(f"c{i}_{k}", f"f{i}", f"m{i}"))
    return vals, trios


def test_the_bootstrap_resamples_families_so_sibships_widen_the_interval():
    rng = np.random.default_rng(3)
    est, width, width_trio, cover = [], [], [], 0
    for rep in range(40):
        vals, trios = _sibships(rng, 100, 4)
        r = transmission(vals, trios, None, n_perm=0, n_boot=200, seed=rep)
        assert r["n_trios"] == 400 and r["n_families"] == 100
        lo, hi = r["reliability_midparent_ci95"]
        est.append(r["reliability_midparent"])
        width.append(hi - lo)
        cover += lo <= 0.6 <= hi
        # the same values with every trio its own family: what resampling trios gave
        v2, t2 = dict(vals), []
        for j, x in enumerate(trios):
            v2[f"F{j}"], v2[f"M{j}"] = vals[x.father], vals[x.mother]
            t2.append(Trio(x.child, f"F{j}", f"M{j}"))
        r2 = transmission(v2, t2, None, n_perm=0, n_boot=200, seed=rep)
        assert r2["n_families"] == 400 and r2["reliability_midparent"] == r["reliability_midparent"]
        width_trio.append(np.subtract(*r2["reliability_midparent_ci95"][::-1]))
    target = 2 * 1.96 * np.std(est, ddof=1)                # the interval's width if it were right
    assert np.mean(width) > 1.15 * np.mean(width_trio)
    assert abs(np.mean(width) / target - 1) < abs(np.mean(width_trio) / target - 1)
    assert cover >= 34                                      # 40 x 0.95 = 38 expected


def test_compare_drops_non_finite_values_and_gives_no_p_for_a_nan_difference():
    rng = np.random.default_rng(1)
    vals, trios = _sibships(rng, 50, 1)
    noisy = {s: v + rng.normal(0, 0.5) for s, v in vals.items()}
    a = dict(vals, c0_0=float("nan"))
    r = compare(a, noisy, trios, None, n_boot=300)
    assert r["n_trios"] == 49 and np.isfinite(r["delta"]) and np.isfinite(r["p_a_better"]) and r["n_families"] == 49
    assert transmission(a, trios, None, n_perm=0, n_boot=0)["n_trios"] == 49
    const = compare({s: 1.0 for s in vals}, noisy, trios, None, n_boot=100)
    assert np.isnan(const["delta"]) and np.isnan(const["p_a_better"])


def test_population_file_header_without_a_group_column_is_refused(tmp_path):
    # a .psam with a sample column but no population column: its second column (SEX) is not a group
    p = tmp_path / "x.psam"
    p.write_text("#IID\tSEX\tPHENO1\nA\t1\t2\nB\t2\t1\n")
    with pytest.raises(ValueError, match="no population/group column"):
        load_population(p)
    # with both, the finer grouping is taken whatever the column order
    q = tmp_path / "y.psam"
    q.write_text("#IID\tSEX\tSuperPop\tPopulation\nA\t1\tAFR\tYRI\nB\t2\tEUR\tCEU\n")
    assert load_population(q) == {"A": "YRI", "B": "CEU"}
