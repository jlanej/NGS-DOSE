"""`ngsdose trios` on the pedigrees people have and the columns cohorts have: three pedigree layouts read
the same trios, a column that is constant across the parents (EBV in blood-derived DNA) gets NaN rather
than killing the table, a column with too few complete trios gets a note, and the others are unaffected."""
import subprocess
import sys
from pathlib import Path

import numpy as np

from ngsdose.trios import Trio, load_pedigree, transmission

ROOT = Path(__file__).resolve().parents[1]
N = 10


def _people():
    return [(f"C{i}", f"F{i}", f"M{i}", "POPA" if i % 2 else "POPB") for i in range(N)]


def test_three_pedigree_layouts_read_the_same_trios(tmp_path):
    fams = _people()
    g1k = tmp_path / "1000G.txt"
    g1k.write_text("FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
                   + "".join(f"F{i} {c} {f} {m} 1 {p} SUP\nF{i} {f} 0 0 1 {p} SUP\nF{i} {m} 0 0 2 {p} SUP\n" for i, (c, f, m, p) in enumerate(fams)))
    ped = tmp_path / "plink.fam"                          # no header; column six is a phenotype code, not a population
    ped.write_text("".join(f"F{i} {c} {f} {m} 1 -9\nF{i} {f} 0 0 1 2\nF{i} {m} 0 0 2 1\n" for i, (c, f, m, _) in enumerate(fams)))
    tab = tmp_path / "trios.tsv"
    tab.write_text("child\tfather\tmother\tpopulation\n" + "".join(f"{c}\t{f}\t{m}\t{p}\n" for c, f, m, p in fams))
    bare = tmp_path / "bare.tsv"
    bare.write_text("".join(f"{c}\t{f}\t{m}\n" for c, f, m, _ in fams))
    want = [Trio(c, f, m, p) for c, f, m, p in fams]
    t1, pop1 = load_pedigree(g1k)
    assert t1 == want and pop1["C1"] == "POPA" and pop1["F0"] == "POPB"
    t2, pop2 = load_pedigree(ped)
    assert [(t.child, t.father, t.mother) for t in t2] == [(c, f, m) for c, f, m, _ in fams] and set(pop2.values()) == {""}
    t3, pop3 = load_pedigree(tab)
    assert t3 == want and pop3["M2"] == "POPB"                   # a trios table gives the parents the child's population
    t4, _ = load_pedigree(bare)
    assert [(t.child, t.father, t.mother) for t in t4] == [(c, f, m) for c, f, m, _ in fams]


def test_a_constant_column_gives_nan_and_the_table_survives(tmp_path):
    rng = np.random.default_rng(11)
    fams = _people()
    ped = tmp_path / "ped.txt"
    ped.write_text("child father mother\n" + "".join(f"{c} {f} {m}\n" for c, f, m, _ in fams))
    rows = []
    for c, f, m, _ in fams:
        vf, vm = rng.normal(100, 10), rng.normal(100, 10)
        rows += [(f, vf, 0.0, ""), (m, vm, 0.0, ""), (c, (vf + vm) / 2 + rng.normal(0, 3), 0.0, "")]
    for s, _, _, _ in rows[:6]:                                # `sparse` has values in two trios only
        pass
    table = tmp_path / "table.tsv"
    with open(table, "w") as fh:
        fh.write("sample\tvar\tconst\tsparse\n")
        for i, (s, v, k, _) in enumerate(rows):
            fh.write(f"{s}\t{v:.4f}\t{k}\t{v:.4f}\n" if i < 6 else f"{s}\t{v:.4f}\t{k}\t\n")
    # the library: NaN statistics, no exception
    const = transmission({s: 0.0 for s, _, _, _ in rows}, [Trio(c, f, m) for c, f, m, _ in fams], None, n_perm=0, n_boot=0)
    assert const["n_trios"] == N and all(np.isnan(const[k]) for k in ("reliability_midparent", "midparent_slope", "spousal_r", "error_cv", "reliability_mendel"))
    # the command: every column gets its row
    out = subprocess.run([sys.executable, "-m", "ngsdose", "trios", str(table), "-p", str(ped), "-c", "var", "const", "sparse", "--perm", "0",
                          "--json", str(tmp_path / "t.json")], check=True, cwd=ROOT, capture_output=True, text=True)
    lines = {l.split("\t")[0]: l for l in out.stdout.splitlines() if l and not l.startswith("column")}
    assert float(lines["var"].split("\t")[2]) > 0.5                       # inherited by construction
    assert lines["const"].split("\t")[2] == "nan"
    assert lines["sparse"].split("\t")[2] == "NA" and "# need at least 3 complete trios" in lines["sparse"]
    import json
    j = json.loads((tmp_path / "t.json").read_text())
    assert "error" in j["sparse"] and j["var"]["n_trios"] == N
