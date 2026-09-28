#!/usr/bin/env python3
"""Read recall of a panel on every source copy of its classes: 150-bp reads tiled every 10 bp over
each copy (both assemblies), and the share that carry >= 4 k-mers of the class (the engine's
assignment rule), with the median number of hits. Also, per read, the classes of OTHER panels whose
k-mers it carries (--other), to see whether a read of the new class could touch a shipped class.

    recall_tiled.py PANEL.tsv.gz COPIES.tsv [--other shipped1.tsv.gz ...] > recall.tsv

COPIES.tsv: class <TAB> source label <TAB> FASTA path <TAB> region (samtools faidx) [<TAB> union-group]
"""
import argparse, gzip, subprocess, sys
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('panel'); ap.add_argument('copies')
ap.add_argument('--other', nargs='*', default=[])
ap.add_argument('-R', type=int, default=150); ap.add_argument('--step', type=int, default=10); ap.add_argument('--min-hits', type=int, default=4)
a = ap.parse_args()
COMP = str.maketrans('ACGT', 'TGCA')

def load(path):
    names, km, k = {}, {}, 31
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##k='): k = int(line[4:])
            elif line.startswith('##class'):
                f = dict(x.split('=', 1) for x in line.rstrip('\n').split('\t')[1:]); names[f['id']] = f['name']
            elif not line.startswith('#'):
                p = line.split('\t'); km[p[0]] = names[p[1]]
    return km, k
panel, k = load(a.panel)
others = {}
for o in a.other:
    km, _ = load(o); others.update(km)

def canon(s):
    r = s.translate(COMP)[::-1]; return s if s < r else r

def per_pos(seq, table):
    out = [None] * max(0, len(seq) - k + 1)
    for i in range(len(out)):
        w = seq[i:i + k]
        if 'N' in w: continue
        out[i] = table.get(canon(w))
    return out

rows = [l.rstrip('\n').split('\t') for l in open(a.copies) if l.strip() and not l.startswith('#')]
print('class\tsource\tregion\treads\trecall_ge4\tmedian_hits\tunion_recall\treads_touching_other_panels\tother_classes')
agg = {}
for r in rows:
    cls, src, fa, reg = r[:4]
    group = set(r[4].split(',')) if len(r) > 4 and r[4] else {cls}
    seq = ''.join(subprocess.run(['samtools', 'faidx', fa, reg], capture_output=True, text=True, check=True).stdout.split('\n')[1:]).upper()
    pp = per_pos(seq, panel)
    op = per_pos(seq, others) if others else [None] * len(pp)
    nk = a.R - k + 1
    own = np.array([1 if c == cls else 0 for c in pp]); grp = np.array([1 if c in group else 0 for c in pp])
    oth = np.array([1 if c else 0 for c in op])
    cs, cg, co = np.r_[0, np.cumsum(own)], np.r_[0, np.cumsum(grp)], np.r_[0, np.cumsum(oth)]
    starts = np.arange(0, len(seq) - a.R + 1, a.step)
    starts = np.array([s for s in starts if 'N' not in seq[s:s + a.R]])
    if len(starts) == 0:
        continue
    h = cs[starts + nk] - cs[starts]; g = cg[starts + nk] - cg[starts]; o = co[starts + nk] - co[starts]
    ocls = sorted({c for c in op if c})
    rec, un = float((h >= a.min_hits).mean()), float((g >= a.min_hits).mean())
    print(f'{cls}\t{src}\t{reg}\t{len(starts)}\t{rec:.3f}\t{int(np.median(h))}\t{un:.3f}\t{int((o > 0).sum())}\t{",".join(ocls)}')
    key = (cls, src if 'crosscheck' in src else src.split(':')[0])
    s = agg.setdefault(key, [0, 0, 0]); s[0] += int((h >= a.min_hits).sum()); s[1] += len(starts); s[2] += int((g >= a.min_hits).sum())
print('class\tsource\treads_ge4/reads\trecall\tunion_recall', file=sys.stderr)
for (c, s), (x, n, u) in sorted(agg.items()):
    print(f'{c}\t{s}\t{x}/{n}\t{x / n:.3f}\t{u / n:.3f}', file=sys.stderr)
