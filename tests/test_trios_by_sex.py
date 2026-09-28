"""Transmission by the sex of parent and child (`trios.by_sex`) on simulated families: an autosomal quantity
passes half of each parent's deviation to every child; a Y-linked one all of the father's to his sons and
none to his daughters, whatever the difference in level between men and women."""
import numpy as np

import pytest

from ngsdose.trios import Trio, by_sex, centre_within_sex


def _families(n=400, seed=7):
    rng = np.random.default_rng(seed)
    trios, sex, pop, auto, ylink = [], {}, {}, {}, {}
    for i in range(n):
        f, m, c = f"F{i}", f"M{i}", f"C{i}"
        child_sex = "M" if i % 2 == 0 else "F"
        sex.update({f: "M", m: "F", c: child_sex})
        pop.update({f: "P1" if i % 3 else "P2", m: "P1" if i % 3 else "P2", c: "P1" if i % 3 else "P2"})
        # autosomal: two haplotype values per parent, one of each parent's passed to the child
        hf, hm = rng.normal(0, 1, 2), rng.normal(0, 1, 2)
        auto[f], auto[m], auto[c] = hf.sum(), hm.sum(), rng.choice(hf) + rng.choice(hm)
        # Y-linked: the father's Y passes whole to his sons; everyone carries a small part elsewhere that is not inherited here
        y = rng.normal(10, 3)
        ylink[f] = y + rng.normal(1, 0.2)
        ylink[m] = rng.normal(1, 0.2)
        ylink[c] = (y if child_sex == "M" else 0.0) + rng.normal(1, 0.2)
        trios.append(Trio(c, f, m, pop[c]))
    return trios, sex, pop, auto, ylink


def test_an_autosomal_quantity_passes_half_to_every_child():
    trios, sex, pop, auto, _ = _families()
    a = by_sex(auto, trios, pop, sex)
    for key in ("father_son", "father_daughter", "mother_son", "mother_daughter"):
        assert a[key]["n"] == 200 and abs(a[key]["slope"] - 0.5) < 0.15, (key, a[key])
    assert abs(a["father_contrast"]["z"]) < 3


def test_a_y_linked_quantity_passes_from_fathers_to_sons_only():
    trios, sex, pop, _, ylink = _families()
    y = by_sex(ylink, trios, pop, sex)
    assert abs(y["father_son"]["slope"] - 1) < 0.05 and y["father_son"]["r"] > 0.95
    assert abs(y["father_daughter"]["slope"]) < 0.1 and abs(y["mother_son"]["slope"]) < 0.5
    assert y["father_contrast"]["z"] > 10


def test_too_few_pairs_gives_nothing():
    trios, sex, pop, auto, _ = _families(n=12)
    assert by_sex(auto, trios, pop, sex) == {}


def test_the_four_pairings_differ_only_for_a_sex_linked_quantity():
    trios, sex, pop, auto, ylink = _families()
    a, y = by_sex(auto, trios, pop, sex, n_perm=200), by_sex(ylink, trios, pop, sex, n_perm=200)
    assert a["heterogeneity"]["p"] > 0.01 and y["heterogeneity"]["p"] < 0.01


def test_a_batch_that_rescales_the_children_moves_the_slope_but_not_the_rescaled_reliability():
    """Every child measured by a later batch reading 1.2x: the midparent slope, and R, read 1.2; R with the
    children on their parents' scale reads 1; and the Mendelian reliability is R + 1 + rho/2 - s^2."""
    from ngsdose.trios import transmission
    trios, _, pop, auto, _ = _families(n=2000, seed=3)
    kids = {t.child for t in trios}
    later = {s: (1.2 * v if s in kids else v) for s, v in auto.items()}
    t = transmission(later, trios, pop, n_perm=0, n_boot=0)
    assert abs(t["midparent_slope"] - 1.2) < 0.06 and abs(t["sd_ratio"] - 1.2) < 0.05
    assert abs(t["reliability_rescaled"] - 1.0) < 0.05
    identity = t["reliability_midparent"] + 1 + t["spousal_r"] / 2 - t["sd_ratio"] ** 2
    assert abs(t["reliability_mendel"] - identity) < 0.03


def test_plink_sex_codes_and_parents_without_a_sex():
    """PLINK's 1/2 read as M/F, and a parent the pedigree gives no sex takes its role's; a sample alone in its
    sex group is centred on everyone, not on itself (which would make it exactly 0)."""
    trios, sex, pop, auto, _ = _families()
    want = by_sex(auto, trios, pop, sex, n_perm=50)
    plink = {s: {"M": "1", "F": "2"}[x] for s, x in sex.items()}
    assert by_sex(auto, trios, pop, plink, n_perm=50) == want
    no_parents = {t.child: sex[t.child] for t in trios}
    assert by_sex(auto, trios, pop, no_parents, n_perm=50) == want
    with pytest.raises(ValueError, match="no child has a sex"):
        by_sex(auto, trios, pop, {})
    v = centre_within_sex({"x": 5, "y": 7, "z": 9, "u": 42}, {}, {"x": "M", "y": "M", "z": "M"}, min_n=3)
    assert v["y"] == 0 and v["u"] == 42 - 63 / 4
