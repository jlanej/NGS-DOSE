#!/usr/bin/env python3
"""Share of each positional unit the estimator can use: position-strands at which a 150-bp read
holds >= 20 panel k-mers (ngsdose.estimate.callable_masks with its default min_kmers=20).

    PYTHONPATH=<NGS-DOSE checkout> usable.py PANEL.tsv.gz
"""
import sys
from ngsdose.io import load_panel
from ngsdose.estimate import callable_masks
p = load_panel(sys.argv[1])
print('class\tlength\tkmers\tusable_fraction_min20')
for c in p.classes.values() if isinstance(p.classes, dict) else p.classes:
    if c.kind != 'positional':
        continue
    f, r = callable_masks(c, p.k, 150, 20)
    print(f'{c.name}\t{c.length}\t{len(c.kmer_pos)}\t{(f.sum() + r.sum()) / (2 * c.length):.3f}')
