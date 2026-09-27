#!/usr/bin/env python3
"""Build the multicopy-genes candidate panel (see build.sh for the environment it expects).

Every class is positional: one unit taken from GRCh38 with samtools faidx (classes.tsv), and a
k-mer is kept only if it occurs once in the unit, nowhere in GRCh38 (analysis set: primary, alts,
decoys, HLA) or CHM13 outside the copies of that class (masks/<class>.<assembly>.bed), and in no
other class of this group or of the shipped panels.

Stage 1, one engine build per class with that class's own masks (a paralog that belongs to a
different class of the group, e.g. SMN2 for SMN1 or HP for HPR, is background for it).
Stage 2, one engine build of every class together, restricted by keep BEDs to the stage-1 k-mers:
k-mers shared between two classes of the group are dropped from both there (the report's
'shared' column), and, with the shipped k-mers as the only background, so is every k-mer of the
shipped panels (bundle panel, satellites, TEL).

Build modes (classes.tsv column 6):
  plain            the unit as it is
  nwin:A-B[,C-D]   unit offsets A..B set to N for the build (C4: the HERV-K(C4) insertion,
                   which is its own class)
  diffN:COPY       every unit base where the paralog COPY (a classes.tsv row) differs set to N:
                   k-mers common to both copies only (SMN: SMN1 and SMN2 summed)
  psvonly:COPY     every unit base more than k-1 bp from such a difference set to N: k-mers
                   specific to the unit (SMN1 against SMN2)
The unit FASTA shipped in units/ is the real sequence (same length); only the build input is N-masked.
"""
import argparse
import concurrent.futures as cf
import gzip
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
K = 31


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def read_classes(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        rows.append(dict(name=p[0], contig=p[1], start=int(p[2]), end=int(p[3]), circular=p[4], mode=p[5],
                         role=p[6], desc=p[7] if len(p) > 7 else ""))
    return rows


def fasta_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def write_fa(path, name, seq, desc=""):
    with open(path, "w") as fh:
        fh.write(f">{name} {desc}\n".rstrip() + "\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


def diffsites(unit_fa, copy_fa):
    out = subprocess.run([sys.executable, str(HERE / "diffsites.py"), str(unit_fa), str(copy_fa)], capture_output=True, text=True, check=True)
    return [int(x) for x in out.stdout.split()]


def kept_positions(panel):
    pos = []
    with gzip.open(panel, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                pos.append(int(line.split("\t")[2]))
    return sorted(pos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grch38", required=True)
    ap.add_argument("--chm13", required=True)
    ap.add_argument("--engine", required=True)
    ap.add_argument("--shipped", nargs="+", required=True, help="shipped panels whose k-mers must not occur here")
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--masks", default=str(HERE / "masks"), help="masks/<class>.<assembly>.bed (masks.py)")
    ap.add_argument("-j", type=int, default=4)
    ap.add_argument("--units-only", action="store_true", help="extract the unit FASTAs and stop (masks.py needs two of them)")
    a = ap.parse_args()
    out = Path(a.out)
    for d in ("units", "aux", "build", "stage1", "keep"):
        (out / d).mkdir(parents=True, exist_ok=True)
    classes = read_classes(HERE / "classes.tsv")
    # 1. units (real sequence) for every row, panel classes and paralog helpers alike
    for c in classes:
        region = f"{c['contig']}:{c['start'] + 1}-{c['end']}"
        seq = subprocess.run(["samtools", "faidx", a.grch38, region], capture_output=True, text=True, check=True).stdout
        seq = "".join(seq.split("\n")[1:]).upper()
        assert len(seq) == c["end"] - c["start"], c["name"]
        # panel classes in units/; paralogs used for checks, and dropped classes, in aux/
        d = "units" if c["role"] == "class" else "aux"
        write_fa(out / d / f"{c['name']}.fa", c["name"], seq, f"GRCh38 {region} {c['desc']}")
    if a.units_only:
        return
    panel_classes = [c for c in classes if c["role"] == "class"]
    # 2. build inputs
    for c in panel_classes:
        seq = list(fasta_seq(out / "units" / f"{c['name']}.fa"))
        mode = c["mode"]
        if mode.startswith("nwin:"):
            for iv in mode[5:].split(","):
                s, e = map(int, iv.split("-"))
                seq[s:e] = "N" * (e - s)
        elif mode.startswith(("diffN:", "psvonly:")):
            kind, other = mode.split(":")
            sites = diffsites(out / "units" / f"{c['name']}.fa", out / "aux" / f"{other}.fa")
            if kind == "diffN":
                for x in sites:
                    seq[x] = "N"
            else:
                near = [False] * len(seq)
                for x in sites:
                    for y in range(max(0, x - K + 1), min(len(seq), x + K)):
                        near[y] = True
                seq = [b if near[i] else "N" for i, b in enumerate(seq)]
            with open(out / "build" / f"{c['name']}.sites.txt", "w") as fh:
                fh.write("\n".join(map(str, sites)) + "\n")
        elif mode != "plain":
            raise SystemExit(f"unknown mode {mode}")
        write_fa(out / "build" / f"{c['name']}.fa", c["name"], "".join(seq), "build input (N-masked where the mode says)")
    # 3. stage 1
    def stage1(c):
        n = c["name"]
        man = out / "stage1" / f"{n}.manifest.tsv"
        man.write_text(f"{n}\tpositional\t{out / 'build' / (n + '.fa')}\t{c['circular']}\n")
        cmd = [a.engine, "panel", "-m", str(man), "-k", str(K), "--max-bg", "0",
               "-b", f"{a.grch38}:{Path(a.masks) / (n + '.GRCh38.bed')}",
               "-b", f"{a.chm13}:{Path(a.masks) / (n + '.CHM13.bed')}",
               "-o", str(out / "stage1" / f"{n}.panel.tsv.gz"), "--report", str(out / "stage1" / f"{n}.rep.tsv.gz")]
        with open(out / "stage1" / f"{n}.log", "w") as log:
            sh(cmd, stderr=log)
        return n
    todo = [c for c in panel_classes if not (out / "stage1" / f"{c['name']}.panel.tsv.gz").exists()]
    with cf.ThreadPoolExecutor(a.j) as ex:
        for n in ex.map(stage1, todo):
            print(f"[stage1] {n}", file=sys.stderr)
    # 4. keep BEDs from the stage-1 k-mers
    for c in panel_classes:
        n = c["name"]
        with open(out / "keep" / f"{n}.keep.bed", "w") as fh:
            for p in kept_positions(out / "stage1" / f"{n}.panel.tsv.gz"):
                fh.write(f"{n}\t{p}\t{p + 1}\n")
    man = out / "manifest.tsv"
    with open(man, "w") as fh:
        for c in panel_classes:
            n = c["name"]
            fh.write(f"{n}\tpositional\tbuild/{n}.fa\t{c['circular']}\t.\tkeep/{n}.keep.bed\n")
    # 5. stage 2: together, no genome background (every k-mer is already background-clean) ...
    sh([a.engine, "panel", "-m", str(man), "-k", str(K), "--max-bg", "0", "-o", str(out / "stage2.nobg.panel.tsv.gz")],
       stderr=open(out / "stage2.nobg.log", "w"))
    # ... and with the shipped panels' k-mers as the background
    shipped_fa = out / "build" / "shipped_kmers.fa"
    kmers = []
    for p in a.shipped:
        with gzip.open(p, "rt") as fh:
            kmers += [l.split("\t", 1)[0] for l in fh if not l.startswith("#")]
    with open(shipped_fa, "w") as fh:
        fh.write(">shipped_panel_kmers\n")
        for i in range(0, len(kmers), 2):
            fh.write("N".join(kmers[i:i + 2]) + "N\n")
    sh([a.engine, "panel", "-m", str(man), "-k", str(K), "--max-bg", "0", "-b", str(shipped_fa),
        "-o", str(out / "panel.tsv.gz"), "--report", str(out / "rep.tsv.gz")],
       stderr=open(out / "stage2.log", "w"))
    with gzip.open(out / "rep.tsv.gz", "rt") as fi, open(out / "rep.tsv", "w") as fo:
        fo.write(fi.read())
    os.remove(out / "rep.tsv.gz")


if __name__ == "__main__":
    main()
