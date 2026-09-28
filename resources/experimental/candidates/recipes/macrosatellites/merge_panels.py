#!/usr/bin/env python3
"""Assemble the candidate panel from per-mask builds and remove the shipped panels' k-mers.

    merge_panels.py MANIFEST OUT_PANEL OUT_REPORT PANEL:REPORT:CLASS[ CLASS...] ... -- SHIPPED_PANEL...

Every build was made from the whole manifest (so k-mers shared between two units are already gone from
both), each with its own background masks; each class is taken from the build named for it. A k-mer of a
shipped panel is then removed. Writes the panel in ngs-dose panel format (loaded and checked by the
engine in the gate), the concatenated per-class report rows, and on stdout the per-class accounting:
distinct unit k-mers, dropped as repeated in the unit, shared with another class of this panel, present
in the GRCh38/CHM13 background outside the class's own copies, shared with a shipped panel, kept; then
the pairs of classes whose units share k-mers.
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def read_panel(path):
    header, classes, rows = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##class"):
                classes.append(dict(x.split("=", 1) for x in line.split("\t")[1:]))
            elif line.startswith("##") and not line.startswith("##ngs-dose-panel") and not line.startswith("##k="):
                header.append(line)
            elif line and not line.startswith("#"):
                rows.append(line.split("\t"))
    return header, classes, rows


def unit_kmers(fasta, circular, k=31):
    s = "".join(l.strip().upper() for l in open(fasta) if not l.startswith(">"))
    if circular:
        s += s[: k - 1]
    return {canon(s[i:i + k]) for i in range(len(s) - k + 1) if "N" not in s[i:i + k]}


def main():
    cut = sys.argv.index("--")
    manifest, out_panel, out_rep = sys.argv[1:4]
    specs, shipped_paths = sys.argv[4:cut], sys.argv[cut + 1:]
    order = [l.split("\t")[0] for l in open(manifest) if l.strip() and not l.startswith("#")]
    ship = set()
    for p in shipped_paths:
        ship.update(r[0] for r in read_panel(p)[2])
    source = {}
    for spec in specs:
        panel, rep, names = spec.split(":", 2)
        for c in names.split():
            source[c] = (panel, rep)
    assert sorted(source) == sorted(order), (sorted(source), sorted(order))
    header = ["##max_background_count=0"]
    classdefs, kept_rows, stats = {}, defaultdict(list), {c: defaultdict(int) for c in order}
    rep_out = open(out_rep, "w")
    rep_out.write("#class\tpos\tcount_in_class\tshared\tbackground\tbuild\n")
    done = set()
    for c in order:
        panel, rep = source[c]
        h, cls, rows = read_panel(panel)
        for x in h:
            if x.startswith("##background") and x not in header:
                header.append(x)
        ids = {d["id"]: d for d in cls}
        classdefs[c] = next(d for d in cls if d["name"] == c)
        for r in rows:
            if ids[r[1]]["name"] == c:
                if r[0] in ship:
                    stats[c]["shipped"] += 1
                else:
                    kept_rows[c].append(r)
        build = Path(panel).name.replace("panel.", "").replace(".tsv.gz", "")
        with gzip.open(rep, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] != c:
                    continue
                rep_out.write(line.rstrip("\n") + f"\t{build}\n")
                s = stats[c]
                s["input"] += 1
                if f[3] == "1":
                    s["shared"] += 1
                elif int(f[2]) > 1:
                    s["multi"] += 1
                elif int(f[4]) > 0:
                    s["bg"] += 1
        done.add(c)
    rep_out.close()
    header.append("##filter=k-mers of the shipped panels removed: " + ",".join(Path(p).name for p in shipped_paths))
    header.append("##note=each class built with only its own copies masked in each background (masks/<assembly>.<build>.bed)")
    seen = {}
    with gzip.open(out_panel, "wt") as w:
        w.write("##ngs-dose-panel v1\n##k=31\n")
        for x in header:
            w.write(x + "\n")
        for i, c in enumerate(order):
            d = classdefs[c]
            w.write(f"##class\tid={i}\tname={c}\tkind={d['kind']}\tlength={d['length']}\tcircular={d['circular']}\tsource={d['source']}"
                    f"\tkmers_input={d['kmers_input']}\tkmers_kept={len(kept_rows[c])}\n")
        w.write("#kmer\tclass\tpos\tstrand\n")
        out = []
        for i, c in enumerate(order):
            for r in kept_rows[c]:
                if r[0] in seen:
                    raise SystemExit(f"k-mer {r[0]} kept in {seen[r[0]]} and {c}")
                seen[r[0]] = c
                out.append((i, int(r[2]), r[0], r[3]))
        for i, pos, km, strand in sorted(out):
            w.write(f"{km}\t{i}\t{pos}\t{strand}\n")
    print("class\tbuild_masks\tunit_kmers\tmulticopy_in_unit\tshared_in_panel\tbackground_genomes\tshared_with_shipped\tkept")
    for c in order:
        s = stats[c]
        b = Path(source[c][0]).name.replace("panel.", "").replace(".tsv.gz", "")
        assert s["input"] - s["shared"] - s["multi"] - s["bg"] == len(kept_rows[c]) + s["shipped"], c
        print(f"{c}\t{b}\t{s['input']}\t{s['multi']}\t{s['shared']}\t{s['bg']}\t{s['shipped']}\t{len(kept_rows[c])}")
    base = Path(manifest).resolve().parent
    uk = {}
    for line in open(manifest):
        p = line.rstrip("\n").split("\t")
        uk[p[0]] = unit_kmers(base / p[2], p[3] == "1")
    for i, a in enumerate(order):
        for b in order[i + 1:]:
            n = len(uk[a] & uk[b])
            if n:
                print(f"#units_share\t{a}\t{b}\t{n}")
    for c in order:
        n = len(uk[c] & ship)
        if n:
            print(f"#unit_kmers_in_shipped_panels\t{c}\t{n}")


if __name__ == "__main__":
    main()
