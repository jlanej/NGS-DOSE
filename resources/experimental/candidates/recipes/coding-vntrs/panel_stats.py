#!/usr/bin/env python3
"""panel_stats.py PANEL NOSHIPPED_PANEL REPORT SHIPPED_PANEL...

Per class: candidate k-mers, k-mers dropped as shared with another class of this panel, dropped as
background (GRCh38/CHM13 outside the masks), removed because a shipped panel holds them (kept in the
panel built without the shipped k-mers, absent from the candidate), and kept. Also checks that the
candidate panel shares no k-mer with any shipped panel."""
import gzip
import sys


def load(path):
    names, km = {}, {}
    for line in gzip.open(path, "rt"):
        if line.startswith("##class"):
            f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
            names[f["id"]] = f["name"]
        elif not line.startswith("#"):
            p = line.split("\t")
            km[p[0]] = names[p[1]]
    return km


def main():
    panel, nosh, rep, shipped = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    b, a = load(panel), load(nosh)
    ship = set()
    for p in shipped:
        ship |= set(load(p))
    stats = {}
    for line in open(rep):
        if line.startswith("#"):
            continue
        cls, _pos, _cnt, shared, bg = line.rstrip("\n").split("\t")
        s = stats.setdefault(cls, [0, 0, 0])
        s[0] += 1
        if shared == "1":
            s[1] += 1
        elif int(bg) > 0:
            s[2] += 1
    overlap = sum(1 for k in b if k in ship)
    print("class\tkmers_input\tshared_in_group\tbackground_or_shipped\tremoved_shared_with_shipped\tkept")
    for cls in sorted(stats):
        ka = sum(1 for v in a.values() if v == cls)
        kb = sum(1 for v in b.values() if v == cls)
        removed = sum(1 for k, v in a.items() if v == cls and k in ship)
        print(f"{cls}\t{stats[cls][0]}\t{stats[cls][1]}\t{stats[cls][2]}\t{removed}\t{kb}")
        assert ka - removed == kb, (cls, ka, removed, kb)
    print(f"# candidate k-mers found in the shipped panels: {overlap} (must be 0)")
    if overlap:
        sys.exit(1)


if __name__ == "__main__":
    main()
