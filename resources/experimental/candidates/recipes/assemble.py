#!/usr/bin/env python3
"""Assemble the candidate panels from the group builds, so that any combination of them loads.

    assemble.py --shipped P... --out DIR --shared TSV NAME=PANEL[:CLASS,CLASS...] ...

Each NAME=PANEL becomes DIR/NAME.k31.panel.tsv.gz; ':CLASS,...' keeps only those classes of the panel
(renumbered in the order given), which is how one group build is split into two loadable files.

The engine merges the panels it is given (-p, repeated) and drops a k-mer held by two panels from both
(src/panel.rs, Panel::load_many), and refuses a class name that two panels define. So:
  - a candidate panel must share no k-mer with a shipped panel (bundle, satellite, telomere): loading
    it would otherwise change a shipped class's counts. Such a k-mer is an error here, not a removal:
    every group build already removes them, and this checks it;
  - a k-mer found in two candidate panels would be counted when either is loaded alone and dropped
    when both are loaded. It is removed from both here, so that a class's k-mers, and so its counts,
    do not depend on which other candidate panels a run loads. --shared lists them per class pair;
  - class names must be unique over all panels, shipped and candidate;
  - a k-mer that is periodic with a period of 1-6 bp, allowing two mismatches (a microsatellite, or a
    telomere-repeat variant such as TTAGGG/TCAGGG/TTGGGG runs), is removed. Such k-mers can be absent
    from GRCh38 and CHM13 and so survive the group builds' background filter, yet every genome's reads
    carry them: the telomere-like direct repeats of HHV-7 and HHV-6 gave 283 HHV7 reads at chromosome
    ends in NA12878, which carries no HHV-7. The rule removes 55 viral k-mers and 23 in human classes
    (DXZ4, D4Z4, the 4qB end's telomere junction, EPPK1, DEFB, RASA4).
The ##class lines get their kmers_kept updated and each output gets a ##filter line saying what was
removed. Output is gzip with a zero timestamp, so that a rebuild gives the same bytes (the engine
records the panel's sha256 in every counts file, and the estimator matches experimental panels by it).
"""
from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def read_panel(path):
    """(header lines, class lines, rows as (kmer bytes, class id, rest of line))."""
    head, classes, kmers, cls, rest = [], [], [], [], []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##class\t"):
                classes.append(line)
            elif line.startswith("#"):
                head.append(line)
            elif line:
                k, c, r = line.split("\t", 2)
                kmers.append(k)
                cls.append(int(c))
                rest.append(r)
    return head, classes, np.array(kmers, dtype="S"), np.array(cls, dtype=np.int32), rest


def periodic(k, max_period=6, max_mismatch=2):
    """True for each k-mer (bytes array) that repeats with some period 1..max_period bp, allowing
    max_mismatch positions i where base i differs from base i + period."""
    m = np.frombuffer(k.tobytes(), dtype=np.uint8).reshape(len(k), -1)
    out = np.zeros(len(k), dtype=bool)
    for p in range(1, max_period + 1):
        out |= (m[:, :-p] != m[:, p:]).sum(axis=1) <= max_mismatch
    return out


def class_fields(line):
    return dict(f.split("=", 1) for f in line.split("\t")[1:] if "=" in f)


def set_field(line, key, value):
    return "\t".join(f"{key}={value}" if f.startswith(key + "=") else f for f in line.split("\t"))


def shipped_kmers(paths):
    out = []
    names = []
    for p in paths:
        _, classes, k, _, _ = read_panel(p)
        out.append(k)
        names += [class_fields(c)["name"] for c in classes]
    return np.unique(np.concatenate(out)), names


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--shipped", nargs="+", required=True, help="the shipped panels (bundle, satellites, telomere)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--shared", required=True, help="TSV of the k-mers removed as shared between candidate panels")
    ap.add_argument("panels", nargs="+", metavar="NAME=PANEL[:CLASSES]")
    a = ap.parse_args(argv)
    ship, ship_names = shipped_kmers(a.shipped)
    groups = []
    for spec in a.panels:
        name, _, rest = spec.partition("=")
        path, _, sub = rest.partition(":")
        head, classes, k, c, r = read_panel(path)
        names = [class_fields(x)["name"] for x in classes]
        if sub:
            want = sub.split(",")
            miss = [w for w in want if w not in names]
            if miss:
                sys.exit(f"error: {path} does not define {', '.join(miss)}")
            old = [names.index(w) for w in want]
            remap = np.full(len(names), -1, dtype=np.int32)
            remap[old] = np.arange(len(old), dtype=np.int32)
            keep = remap[c] >= 0
            k, c, r = k[keep], remap[c[keep]], [x for x, y in zip(r, keep) if y]
            classes = [set_field(classes[o], "id", i) for i, o in enumerate(old)]
            names = want
        groups.append(dict(name=name, path=path, head=head, classes=classes, names=names, k=k, c=c, r=r))

    # class names: unique over shipped and candidate panels
    seen = Counter(ship_names + [n for g in groups for n in g["names"]])
    dup = sorted(n for n, x in seen.items() if x > 1)
    if dup:
        sys.exit(f"error: class name(s) defined twice: {', '.join(dup)}")

    # no candidate k-mer in a shipped panel
    bad = []
    for g in groups:
        hit = np.isin(g["k"], ship)
        if hit.any():
            bad.append(f"{g['name']}: {int(hit.sum())}")
    if bad:
        sys.exit("error: candidate k-mers found in a shipped panel (the group build must remove them): " + "; ".join(bad))

    # k-mers in more than one candidate panel: removed from all of them
    allk = np.concatenate([g["k"] for g in groups])
    u, cnt = np.unique(allk, return_counts=True)
    shared = u[cnt > 1]
    pairs: Counter = Counter()
    if len(shared):
        who = {}                                            # shared k-mer -> [(group, class)]
        for g in groups:
            idx = np.nonzero(np.isin(g["k"], shared))[0]
            for i in idx:
                who.setdefault(g["k"][i], []).append((g["name"], g["names"][g["c"][i]]))
        for owners in who.values():
            for i in range(len(owners)):
                for j in range(i + 1, len(owners)):
                    x, y = sorted([owners[i], owners[j]])
                    pairs[(x, y)] += 1
    Path(a.out).mkdir(parents=True, exist_ok=True)
    Path(a.shared).parent.mkdir(parents=True, exist_ok=True)
    with open(a.shared, "w") as fh:
        fh.write("# k-mers found in two candidate panels, removed from both (assemble.py)\n")
        fh.write("group_a\tclass_a\tgroup_b\tclass_b\tkmers\n")
        for ((ga, ca), (gb, cb)), n in sorted(pairs.items(), key=lambda x: (-x[1], x[0])):
            fh.write(f"{ga}\t{ca}\t{gb}\t{cb}\t{n}\n")

    out = Path(a.out)
    print("panel\tclass\tkmers_in\tremoved_shared_with_other_candidate_panels\tremoved_periodic\tkmers_out")
    for g in groups:
        sh = np.isin(g["k"], shared) if len(shared) else np.zeros(len(g["k"]), dtype=bool)
        per = periodic(g["k"]) & ~sh
        drop = sh | per
        n_in = np.bincount(g["c"], minlength=len(g["names"]))
        n_sh = np.bincount(g["c"][sh], minlength=len(g["names"]))
        n_per = np.bincount(g["c"][per], minlength=len(g["names"]))
        n_drop = n_sh + n_per
        classes = [set_field(line, "kmers_kept", int(n_in[i] - n_drop[i])) for i, line in enumerate(g["classes"])]
        # a panel assembled before (a group not rebuilt, taken from DEST) keeps the counts of that assembly
        # added to this one's, so that re-assembling an installed file gives the same bytes
        prev_sh = prev_per = 0
        for h in g["head"]:
            if h.startswith("##filter=assemble.py: "):
                m = re.match(r"##filter=assemble\.py: (\d+) k-mers found .*?; (\d+) k-mers periodic", h)
                prev_sh, prev_per = int(m.group(1)), int(m.group(2))
        head = [h for h in g["head"] if not h.startswith(("#kmer", "##filter=assemble.py"))]
        head.append(f"##filter=assemble.py: {prev_sh + int(sh.sum())} k-mers found in another candidate panel removed (from both; "
                    f"resources/experimental/candidates/shared_between_panels.tsv); {prev_per + int(per.sum())} k-mers periodic "
                    "with a period of 1-6 bp and at most 2 mismatches removed")
        buf = io.StringIO()
        for h in head + classes:
            buf.write(h + "\n")
        buf.write("#kmer\tclass\tpos\tstrand\n")
        keep = ~drop
        for kk, cc, rr in zip(g["k"][keep], g["c"][keep], (x for x, y in zip(g["r"], keep) if y)):
            buf.write(f"{kk.decode()}\t{cc}\t{rr}\n")
        dest = out / f"{g['name']}.k31.panel.tsv.gz"
        with open(dest, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as z:
            z.write(buf.getvalue().encode())
        for i, n in enumerate(g["names"]):
            print(f"{dest.name}\t{n}\t{n_in[i]}\t{n_sh[i]}\t{n_per[i]}\t{n_in[i] - n_drop[i]}")


if __name__ == "__main__":
    main()
