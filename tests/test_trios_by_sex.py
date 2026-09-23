"""Transmission by the sex of parent and child (`trios.by_sex`) on simulated families: an autosomal quantity
passes half of each parent's deviation to every child; a Y-linked one all of the father's to his sons and
none to his daughters, whatever the difference in level between men and women."""
import numpy as np

from ngsdose.trios import Trio, by_sex


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
