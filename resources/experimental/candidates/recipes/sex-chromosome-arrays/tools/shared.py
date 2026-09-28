#!/usr/bin/env python3
"""K-mers of a candidate panel that also occur in other (shipped) panels, per candidate class and
shipped class; with --fasta, write every shipped k-mer as a 31-bp FASTA record (a background that
removes them at build time)."""
import gzip, sys, collections
def load(path):
    names, km = {}, {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##class'):
                f = dict(x.split('=', 1) for x in line.rstrip('\n').split('\t')[1:]); names[f['id']] = f['name']
            elif not line.startswith('#'):
                p = line.split('\t'); km[p[0]] = names[p[1]]
    return km
args = sys.argv[1:]
if args[0] == '--fasta':
    out = args[1]; shipped = args[2:]
    with open(out, 'w') as fh:
        n = 0
        for s in shipped:
            for k, c in load(s).items():
                fh.write(f'>{c}_{n}\n{k}\n'); n += 1
    print(f'{n} shipped k-mers written to {out}', file=sys.stderr); sys.exit()
cand = load(args[0]); ship = {}
for s in args[1:]:
    ship.update(load(s))
tab = collections.Counter((c, ship[k]) for k, c in cand.items() if k in ship)
tot = collections.Counter(cand.values())
for c in sorted(tot):
    parts = {s: n for (cc, s), n in tab.items() if cc == c}
    print(c, tot[c], sum(parts.values()), ' '.join(f'{s}:{n}' for s, n in sorted(parts.items())), sep='\t')
