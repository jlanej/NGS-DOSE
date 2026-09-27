#!/usr/bin/env python3
"""Copies of a unit in a genome extract: map chunks of the unit with minimap2 (asm20, many
secondaries), cluster the hits per contig (gap <= --gap), and keep clusters holding a real share of
the unit (distinct chunks >= --min-frac of all chunks, and >= 3), which leaves out interspersed-repeat
hits. Prints BED (padded) with the number of distinct chunks, hits and median identity."""
import argparse, statistics, subprocess, sys
ap = argparse.ArgumentParser()
ap.add_argument('unit'); ap.add_argument('genome'); ap.add_argument('name')
ap.add_argument('--size', type=int, default=1000); ap.add_argument('--step', type=int, default=500)
ap.add_argument('--min-id', type=float, default=0.90); ap.add_argument('--gap', type=int, default=5000)
ap.add_argument('--min-frac', type=float, default=0.25); ap.add_argument('--pad', type=int, default=1000)
a = ap.parse_args()
seq = ''.join(l.strip() for l in open(a.unit) if not l.startswith('>')).upper()
n = 0
with open(f'{a.name}.chunks.fa', 'w') as fh:
    for i in range(0, max(1, len(seq) - a.size + 1), a.step):
        c = seq[i:i + a.size]
        if c.count('N') > 0.5 * len(c): continue
        fh.write(f'>c{i}\n{c}\n'); n += 1
paf = subprocess.run(['minimap2', '-c', '-x', 'asm20', '-N', '200', '-p', '0.1', '--secondary=yes', '-t', '8', a.genome, f'{a.name}.chunks.fa'],
                     capture_output=True, text=True).stdout
open(f'{a.name}.paf', 'w').write(paf)
hits = []
for l in paf.splitlines():
    p = l.split('\t')
    if int(p[10]) < 0.5 * min(a.size, len(seq)): continue
    idt = int(p[9]) / int(p[10])
    if idt < a.min_id: continue
    hits.append((p[5], int(p[7]), int(p[8]), p[0], idt, int(p[6])))
hits.sort()
cl = []
for h in hits:
    if cl and cl[-1]['c'] == h[0] and h[1] <= cl[-1]['e'] + a.gap:
        x = cl[-1]; x['e'] = max(x['e'], h[2]); x['q'].add(h[3]); x['id'].append(h[4]); x['n'] += 1
    else:
        cl.append(dict(c=h[0], s=h[1], e=h[2], q={h[3]}, id=[h[4]], n=1, L=h[5]))
need = max(3, a.min_frac * n)
kept = 0
for x in cl:
    if len(x['q']) >= need:
        kept += 1
        print(f"{x['c']}\t{max(0, x['s'] - a.pad)}\t{min(x['L'], x['e'] + a.pad)}\t{a.name}\tchunks={len(x['q'])}/{n};hits={x['n']};med_id={statistics.median(x['id']):.3f};span={x['e']-x['s']}")
    elif len(x['q']) >= 2:
        print(f"#dropped {x['c']}:{x['s']}-{x['e']} chunks={len(x['q'])}/{n} med_id={statistics.median(x['id']):.3f}", file=sys.stderr)
print(f'{a.name}: {kept} loci', file=sys.stderr)
