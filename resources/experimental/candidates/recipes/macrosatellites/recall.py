#!/usr/bin/env python3
"""Read recall of a positional candidate panel on every source copy of each class.

    recall.py PANEL.tsv.gz SOURCES.tsv > recall.tsv

SOURCES.tsv: class, label, fasta, region, mode, [record regex]. region is a samtools region or '.' for
every record of the FASTA (optionally only records whose name matches the regex). mode 'tile' tiles
the whole region; 'locate:ID' first finds the copies of the class's unit in it (minimap2 asm20 of 1-kb
unit chunks, identity >= ID, hits merged across gaps <= 2 kb, loci >= 1 kb or 'locate:ID:MINLEN') and tiles those.
Reads are 150 bp every 10 bp (windows with N skipped). A read counts for a class when it carries >= 4
of that class's k-mers, the engine's rule. Reports the share for the source's own class and, for the
other classes, the share of the same reads that would reach 4 of theirs (cross-talk).
"""
import gzip
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

K, RL, STEP, MIN_HITS = 31, 150, 10, 4
COMP = str.maketrans("ACGT", "TGCA")


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def read_fasta_text(text):
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf).upper()))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf).upper()))
    return out


def get_seqs(fasta, region, rx):
    if region != ".":
        return read_fasta_text(subprocess.run(["samtools", "faidx", fasta, region], check=True, capture_output=True, text=True).stdout)
    recs = read_fasta_text(open(fasta).read())
    return [r for r in recs if not rx or re.search(rx, r[0])]


def locate(unit_fa, seqs, min_id, min_len=1000):
    unit = "".join(l.strip().upper() for l in open(unit_fa) if not l.startswith(">"))
    d = unit + unit
    with tempfile.TemporaryDirectory() as td:
        q, t = os.path.join(td, "q.fa"), os.path.join(td, "t.fa")
        with open(q, "w") as fh:
            for i in range(0, len(unit), 250):
                fh.write(f">c{i}\n{d[i:i + 1000]}\n")
        with open(t, "w") as fh:
            for i, (_, s) in enumerate(seqs):
                fh.write(f">s{i}\n{s}\n")
        paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "5000", "-p", "0.05", "--secondary=yes", "-t", "4", t, q],
                             check=True, capture_output=True, text=True).stdout
    hits = {}
    for line in paf.splitlines():
        p = line.split("\t")
        nm, al = int(p[9]), int(p[10])
        if al >= 300 and nm / al >= min_id:
            hits.setdefault(int(p[5][1:]), []).append((int(p[7]), int(p[8])))
    out = []
    for i, h in hits.items():
        h.sort()
        cs, ce = h[0]
        for s, e in h[1:]:
            if s <= ce + 2000:
                ce = max(ce, e)
            else:
                out.append((i, cs, ce))
                cs, ce = s, e
        out.append((i, cs, ce))
    return [seqs[i][1][s:e] for i, s, e in out if e - s >= min_len]


def main():
    panel, sources = sys.argv[1], sys.argv[2]
    names, kmers = {}, {}
    with gzip.open(panel, "rt") as fh:
        for line in fh:
            if line.startswith("##class"):
                f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
                names[f["id"]] = f["name"]
            elif not line.startswith("#"):
                p = line.split("\t")
                kmers[p[0]] = names[p[1]]
    classes = [names[i] for i in sorted(names, key=int)]
    base = Path(sources).resolve().parent
    print("class\tsource\tcopies_bp\treads\trecall\tmedian_hits\tcross_talk")
    for line in open(sources):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        cls, label, fasta, region, mode = p[:5]
        rx = p[5] if len(p) > 5 else ""
        fasta = str(base / fasta) if not fasta.startswith("/") and not fasta.startswith("$") else os.path.expandvars(fasta)
        seqs = get_seqs(fasta, region, rx)
        if mode.startswith("locate:"):
            m = mode.split(":")
            copies = locate(base / "units" / f"{cls}.fa", seqs, float(m[1]), int(m[2]) if len(m) > 2 else 1000)
        else:
            copies = [s for _, s in seqs]
        n, own, hits, other = 0, 0, [], {c: 0 for c in classes}
        for s in copies:
            for i in range(0, len(s) - RL + 1, STEP):
                w = s[i:i + RL]
                if "N" in w:
                    continue
                cnt = {}
                for j in range(RL - K + 1):
                    c = kmers.get(canon(w[j:j + K]))
                    if c:
                        cnt[c] = cnt.get(c, 0) + 1
                n += 1
                h = cnt.get(cls, 0)
                hits.append(h)
                own += h >= MIN_HITS
                for c, v in cnt.items():
                    if c != cls and v >= MIN_HITS:
                        other[c] += 1
        hits.sort()
        xt = ",".join(f"{c}:{v / n:.3f}" for c, v in other.items() if v) or "none"
        bp = sum(len(s) for s in copies)
        if n:
            print(f"{cls}\t{label}\t{bp}\t{n}\t{own / n:.3f}\t{hits[len(hits) // 2]}\t{xt}")
        else:
            print(f"{cls}\t{label}\t{bp}\t0\tNA\tNA\tnone")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
