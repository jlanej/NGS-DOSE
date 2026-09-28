# fetchplan worked examples

Written by `resources/build/fetch_examples.sh` on 2026-09-28; do not edit by hand.
The READMEs and `docs/DESIGN.md` quote plan figures (MB per genome, intervals kept, expected capture)
from this file: rerun the script after a change to the menu, a sinks or statistics file, the control
regions or `ngsdose fetchplan`. The descriptions in `resources/fetch_menu.tsv` and the headers and
table of `resources/experimental/subsets/` carry rounded costs written when those files were made;
examples 7 and 14 to 16 give those options' own costs (mb_median) to check them against.
Medians of cumulative totals are not additive: the difference between two rows, or between two plans'
totals, is not the median of what an option adds per CRAM.

Costs are what `ngs-dose count -m fetch` reads: the CRAM slices its fetches decode, in MB (1e6 bytes), median over
13 CRAM indexes:
HG00096, HG00706, HG01084, HG01552, HG01871, HG02383, HG02451, HG02466, HG02792, HG03007, HG03072, HG03267, NA12878
(`resources/build/fetch_examples.sh --get DIR` downloads them; sha256 of their concatenation
`a9f6c4d4766a7095`).
These are NYGC 30x CRAMs of the 1000 Genomes high-coverage release (bwa-mem, GRCh38 with decoys
and HLA; contigs from `GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`). The sinks, their statistics and so these costs belong to
that aligner and reference: another pipeline learns its own sinks from its own scans. Percentages are of
each whole CRAM. The engine merges a plan's intervals where they touch or overlap and makes one indexed fetch per
run, and each fetch decodes every slice that overlaps its interval (with its container's compression header),
so a slice under several runs is decoded once per run. In each table, mb_median is the option's own intervals
alone and cum_mb_median the plan up to and including that row, both priced so; cum_mb_floor is the plan with
every slice decoded once, the floor a reader that sorted the plan's slices would reach (the engine does not:
the gap is a few percent in men and about 10-15% in women, whose few sparse chrY slices are decoded once per
chrY truth region). Candidate classes have no sinks yet: a plan lists them for the whole-file scans only.

Engine for count_flags.txt: `ngs-dose` (ngs-dose 0.1.1; the binary's sha256 begins
`4ac2d1fa3c50b819`, as `--version` does not tell builds apart), which takes `count --classes`.
`count --classes` is in engines from the fetch-menu change of 2026-09-26 on; fae1124 lacks it. The engine changes only count_flags.txt, and whether fetchplan accepts a plan whose
panels define classes it does not select (`ngsdose fetchplan --help`, --engine).

| example | MB per genome | % of the CRAM | floor: every slice once |
| --- | ---: | ---: | ---: |
| [1. core: the bundle's positional classes](#example-1) | 503.9 | 3.07 | 482.5 |
| [2. core_tel: core and the telomeric repeat](#example-2) | 546.9 | 3.33 | 522.7 |
| [3. biobank_lite: core_tel and the four cheapest satellites](#example-3) | 600.3 | 3.64 | 573.9 |
| [4. core and the cheap satellites SST1, CER, SATR, ACRO at a capture target of 0.995](#example-4) | 457.0 | 2.77 | 414.7 |
| [5. core with the 200 lite control regions](#example-5) | 369.8 | 2.27 | 338.8 |
| [6. core_tel with a capture target for TEL alone](#example-6) | 528.0 | 3.25 | 508.5 |
| [7. the per-array options DXZ1 and DYZ3](#example-7) | 274.7 | 1.73 | 244.8 |
| [8. core_tel and the tier-A candidates, for the whole-file scans](#example-8) | 546.9 | 3.33 | 522.7 |
| [9. a budget of 600 MB, filled](#example-9) | 596.7 | 3.66 | 572.7 |
| [10. a budget of 1000 MB, filled](#example-10) | 944.9 | 5.81 | 893.5 |
| [11. the ten satellite families and the unmapped bin at a capture target of 0.995](#example-11) | 2840.4 | 17.61 | 2608.2 |
| [12. every option with sinks (each row's mb_median is that option alone)](#example-12) | 3416.7 | 21.59 | 3099.5 |
| [13. core_tel with the lite control regions and a capture target of 0.995](#example-13) | 385.9 | 2.39 | 355.5 |
| [14. biobank_lite at a capture target of 0.995](#example-14) | 567.1 | 3.48 | 544.8 |
| [15. the four per-array options (preset xy_arrays)](#example-15) | 340.7 | 2.05 | 293.9 |
| [16. core_tel and the four per-array options](#example-16) | 623.6 | 3.78 | 581.4 |

<a id="example-1"></a>
## 1. core: the bundle's positional classes

```
ngsdose fetchplan --preset core --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**503.9 MB per genome (3.07% of the CRAM; 482.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- (none)

<a id="example-2"></a>
## 2. core_tel: core and the telomeric repeat

```
ngsdose fetchplan --preset core_tel --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**546.9 MB per genome (3.33% of the CRAM; 522.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- (none)

<a id="example-3"></a>
## 3. biobank_lite: core_tel and the four cheapest satellites

```
ngsdose fetchplan --preset biobank_lite --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**600.3 MB per genome (3.64% of the CRAM; 573.9 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 17.1 | NA | 563.9 | 538.7 | 3.43 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 569.1 | 544.0 | 3.45 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 579.7 | 554.1 | 3.51 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 600.3 | 573.9 | 3.64 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,SST1,CER,SATR,ACRO`

Notes:

- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,SST1,CER,SATR,ACRO (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, CER, SATR, ACRO): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-4"></a>
## 4. core and the cheap satellites SST1, CER, SATR, ACRO at a capture target of 0.995

```
ngsdose fetchplan --preset core --classes SST1 CER SATR ACRO --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**457.0 MB per genome (2.77% of the CRAM; 414.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.9 | 0.0 | 263.3 | 233.1 | 1.64 | per byte | . |
| DJ | shipped | A | 53 | 6 | 0.995 | 0.99580 | 0.00140 | 46.0 | 3.5 | 305.5 | 277.9 | 1.92 | per byte | 6 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for CER), so they cost nothing and add to its capture |
| rDNA45S | shipped | A | 5 | 14 | 0.995 | 0.99816 | 0.00116 | 102.8 | 98.9 | 411.2 | 368.7 | 2.50 | per byte | . |
| SST1 | experimental | C | 28 | 6 | 0.995 | 0.99609 | 0.00346 | 13.7 | 3.3 | 425.8 | 382.3 | 2.58 | per read (fewer bytes than per byte) | . |
| SATR | experimental | D | 20 | 1 | 0.995 | 0.99720 | 0.00178 | 10.6 | 95.4 | 436.8 | 391.5 | 2.65 | per byte | . |
| ACRO | experimental | D | 14 | 8 | 0.995 | 0.99528 | 0.00421 | 19.2 | 97.2 | 454.2 | 409.6 | 2.75 | per byte | . |
| CER | experimental | D | 44 | 5 | 0.995 | 0.99500 | 0.00466 | 27.0 | 2.8 | 457.0 | 414.7 | 2.77 | per read (fewer bytes than per byte) | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for DJ), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,SST1,SATR,ACRO,CER`

Notes:

- capture targets saved 100.3 MB of the plan: 557.3 MB with every option's intervals whole, 457.0 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (rDNA5S, DJ, rDNA45S, SST1, SATR, ACRO, CER) need not add up to the plan's
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,SST1,SATR,ACRO,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, SATR, ACRO, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, SATR, ACRO, CER, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), SST1 (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans): scans not used to learn the sinks

<a id="example-5"></a>
## 5. core with the 200 lite control regions

```
ngsdose fetchplan --preset core --controls resources/GRCh38/controls.lite200.bed --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**369.8 MB per genome (2.27% of the CRAM; 338.8 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 382 | 0 | NA | NA | NA | 120.9 | NA | 120.9 | 93.1 | 0.78 | NA | 382 regions of controls.lite200.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 124.0 | 94.9 | 0.78 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 185.2 | 142.5 | 1.08 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 369.8 | 338.8 | 2.27 | NA | . |

controls.txt (the -c FASTA): controls.lite200.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.lite200.bed, not the menu's controls.bed: the fetch must pass -c controls.lite200.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names

<a id="example-6"></a>
## 6. core_tel with a capture target for TEL alone

```
ngsdose fetchplan --preset core_tel --capture-class TEL=0.99 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**528.0 MB per genome (3.25% of the CRAM; 508.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 55 | 8 | 0.99 | 0.99698 | 0.00092 | 119.6 | 15.2 | 528.0 | 508.5 | 3.25 | per byte | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- capture targets saved 18.8 MB of the plan: 546.9 MB with every option's intervals whole, 528.0 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (TEL) need not add up to the plan's
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for TEL, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans): scans not used to learn the sinks

<a id="example-7"></a>
## 7. the per-array options DXZ1 and DYZ3

```
ngsdose fetchplan --classes DXZ1 DYZ3 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**274.7 MB per genome (1.73% of the CRAM; 244.8 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 5.0 | NA | 264.9 | 236.2 | 1.66 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.8 | NA | 274.7 | 244.8 | 1.73 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=aSatHOR`

Notes:

- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, SATR, ACRO, SST1, CER, not selected: count_flags.txt has --classes=aSatHOR (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-8"></a>
## 8. core_tel and the tier-A candidates, for the whole-file scans

```
ngsdose fetchplan --preset core_tel candidates_A --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**546.9 MB per genome (3.33% of the CRAM; 522.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |
| DXZ4 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| CT47 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| RS447 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| MSR5p | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| FLJ40296 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| D4Z4 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| KIV2 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| C4 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| AMY1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| AMY2B | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SMN | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SMN1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| RHD | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HBA | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| GSTM1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| GSTT1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| UGT2B17 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| LCE3BC | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| APOBEC3B | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HPR | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| CCL3L | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| ORM1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| NPY4R | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| TSPY | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| DAZ | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| BPY2 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| DYZ19 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| OPN1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| TDNA1Q23 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| RNU1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SNORD116 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SNORD115 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HHV6A | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HHV6B | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SMRV | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| VNTR_ACAN | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| VNTR_MUC1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| VNTR_MUC19 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

scan_panels.txt (panels for the whole-file scans): panel.k31.tsv.gz, telomere.k31.panel.tsv.gz, macrosatellites.k31.panel.tsv.gz, multicopy-genes.k31.panel.tsv.gz, sex-chromosome-arrays.k31.panel.tsv.gz, rna-arrays.k31.panel.tsv.gz, nonhuman.k31.panel.tsv.gz, coding-vntrs.k31.panel.tsv.gz

Notes:

- (each of the 38 candidates has this note; the first:) DXZ4 is a candidate: no sinks have been learned for it, so a fetch cannot measure it. Its panel is in scan_panels.txt: load it in the whole-file scans, learn its sinks there (`ngsdose sinks SCANS`), check their capture on held-out scans (`ngsdose sinks HELD_OUT --evaluate BED`), then give it the sinks file in the menu

<a id="example-9"></a>
## 9. a budget of 600 MB, filled

```
ngsdose fetchplan --budget-mb 600 --fill --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**596.7 MB per genome (3.66% of the CRAM; 572.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 5.0 | NA | 551.9 | 525.7 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 555.9 | 534.7 | 3.41 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.8 | NA | 563.9 | 541.9 | 3.48 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 17.1 | NA | 581.0 | 557.4 | 3.57 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 586.1 | 562.6 | 3.60 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 596.7 | 572.7 | 3.66 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,CER,SATR`

Notes:

- 83 candidate option(s) are left out: they have no learned sinks yet
- DYZ2 (tier C) is left out: it would take the plan to 621.8 MB, past the budget of 600 MB
- DYZ1 (tier C) is left out: it would take the plan to 626.0 MB, past the budget of 600 MB
- HSat2 (tier C) is left out: it would take the plan to 877.1 MB, past the budget of 600 MB
- aSatHOR (tier C) is left out: it would take the plan to 1225.4 MB, past the budget of 600 MB
- ACRO (tier D) is left out: it would take the plan to 617.4 MB, past the budget of 600 MB
- HSat1A (tier D) is left out: it would take the plan to 760.9 MB, past the budget of 600 MB
- HSat1B (tier D) is left out: it would take the plan to 1008.6 MB, past the budget of 600 MB
- bSat (tier D) is left out: it would take the plan to 1047.5 MB, past the budget of 600 MB
- HSat3 (tier D) is left out: it would take the plan to 1841.9 MB, past the budget of 600 MB
- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, ACRO, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,CER,SATR (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, SST1, CER, SATR): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-10"></a>
## 10. a budget of 1000 MB, filled

```
ngsdose fetchplan --budget-mb 1000 --fill --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**944.9 MB per genome (5.81% of the CRAM; 893.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 5.0 | NA | 551.9 | 525.7 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 555.9 | 534.7 | 3.41 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.8 | NA | 563.9 | 541.9 | 3.48 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 17.1 | NA | 581.0 | 557.4 | 3.57 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 621.8 | 587.0 | 3.76 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 648.2 | 605.0 | 3.92 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |
| HSat2 | experimental | C | 53 | 0 | NA | 0.99993 | 0.00000 | 374.6 | NA | 911.0 | 858.3 | 5.58 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 916.2 | 863.5 | 5.62 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 927.2 | 873.6 | 5.68 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 944.9 | 893.5 | 5.81 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,HSat1B,HSat3,HSat2,CER,SATR,ACRO`

Notes:

- 83 candidate option(s) are left out: they have no learned sinks yet
- aSatHOR (tier C) is left out: it would take the plan to 1557.1 MB, past the budget of 1000 MB
- HSat1A (tier D) is left out: it would take the plan to 1086.9 MB, past the budget of 1000 MB
- HSat1B (tier D) is left out: it would take the plan to 1276.3 MB, past the budget of 1000 MB
- bSat (tier D) is left out: it would take the plan to 1399.4 MB, past the budget of 1000 MB
- HSat3 (tier D) is left out: it would take the plan to 2114.5 MB, past the budget of 1000 MB
- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- DYZ2 is a named subset of HSat1B's intervals, fetched without the rest of HSat1B's sinks: the fetch counts HSat1B only inside them, so `ngsdose estimate` reports DYZ2.mass_Mb and marks HSat1B itself as not measured (subset_only)
- DYZ1 is a named subset of HSat3's intervals, fetched without the rest of HSat3's sinks: the fetch counts HSat3 only inside them, so `ngsdose estimate` reports DYZ1.mass_Mb and marks HSat3 itself as not measured (subset_only)
- the panels loaded also define HSat1A, bSat, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,HSat1B,HSat3,HSat2,CER,SATR,ACRO (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, SST1, DYZ2, DYZ1, HSat2, CER, SATR, ACRO): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-11"></a>
## 11. the ten satellite families and the unmapped bin at a capture target of 0.995

```
ngsdose fetchplan --preset satellites --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**2840.4 MB per genome (17.61% of the CRAM; 2608.2 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 269.0 | 242.7 | 1.68 | NA | every read without a coordinate (count --unmapped) |
| SST1 | experimental | C | 31 | 3 | 0.995 | 0.99522 | 0.00432 | 15.7 | 1.2 | 286.4 | 258.3 | 1.77 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat3), so they cost nothing and add to its capture |
| HSat2 | experimental | C | 42 | 11 | 0.995 | 0.99539 | 0.00455 | 366.9 | 7.1 | 673.2 | 612.6 | 4.03 | per byte | 6 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat1B, HSat3), so they cost nothing and add to its capture |
| aSatHOR | experimental | C | 208 | 40 | 0.995 | 0.99720 | 0.00158 | 763.2 | 16.0 | 1336.2 | 1235.0 | 8.02 | per byte | 5 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for HSat1A, bSat, HSat1B, HSat3), so they cost nothing and add to its capture |
| SATR | experimental | D | 21 | 0 | 0.995 | 0.99898 | 0.00000 | 104.0 | 0.0 | 1347.2 | 1246.0 | 8.07 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for HSat1B), so they cost nothing and add to its capture |
| ACRO | experimental | D | 16 | 6 | 0.995 | 0.99528 | 0.00421 | 112.1 | 2.4 | 1361.8 | 1258.9 | 8.17 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat1B), so they cost nothing and add to its capture |
| CER | experimental | D | 43 | 6 | 0.995 | 0.99500 | 0.00466 | 26.4 | 3.7 | 1383.5 | 1278.7 | 8.33 | per read (fewer bytes than per byte) | . |
| HSat1A | experimental | D | 17 | 13 | 0.995 | 0.99899 | 0.00088 | 220.8 | 45.4 | 1405.4 | 1296.9 | 8.44 | per byte | 7 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, bSat, HSat1B, HSat3), so they cost nothing and add to its capture |
| bSat | experimental | D | 301 | 107 | 0.995 | 0.99508 | 0.00439 | 469.5 | 79.6 | 1775.7 | 1549.8 | 10.65 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat1B), so they cost nothing and add to its capture |
| HSat1B | experimental | D | 1232 | 0 | 0.995 | 0.97171 | 0.00000 | 507.6 | 0.0 | 2116.2 | 1883.4 | 12.95 | per byte | . |
| HSat3 | experimental | D | 81 | 11 | 0.995 | 0.99545 | 0.00407 | 1016.3 | 344.4 | 2840.4 | 2608.2 | 17.61 | per byte | 12 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for HSat2, aSatHOR, HSat1B), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped`

Notes:

- capture targets saved 442.9 MB of the plan: 3283.3 MB with every option's intervals whole, 2840.4 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (SST1, HSat2, aSatHOR, SATR, ACRO, CER, HSat1A, bSat, HSat1B, HSat3) need not add up to the plan's
- HSat1B: all its intervals together capture 0.9717 (p10 over the scans of its statistics), below the target 0.995: all are kept
- experimental options selected (SST1, HSat2, aSatHOR, SATR, ACRO, CER, HSat1A, bSat, HSat1B, HSat3): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for SST1, HSat2, aSatHOR, SATR, ACRO, HSat1A, bSat, HSat1B, HSat3, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for SST1 (1648 scans), HSat2 (1648 scans), aSatHOR (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans), HSat1A (1648 scans), bSat (1648 scans), HSat1B (1648 scans), HSat3 (1648 scans): scans not used to learn the sinks

<a id="example-12"></a>
## 12. every option with sinks (each row's mb_median is that option alone)

```
ngsdose fetchplan --preset core_tel satellites xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**3416.7 MB per genome (21.59% of the CRAM; 3099.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 5.0 | NA | 551.9 | 525.7 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 555.9 | 534.7 | 3.41 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.8 | NA | 563.9 | 541.9 | 3.48 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 17.1 | NA | 581.0 | 557.4 | 3.57 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 621.8 | 587.0 | 3.76 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 648.2 | 605.0 | 3.92 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |
| HSat2 | experimental | C | 53 | 0 | NA | 0.99993 | 0.00000 | 374.6 | NA | 911.0 | 858.3 | 5.58 | NA | . |
| aSatHOR | experimental | C | 248 | 0 | NA | 0.99878 | 0.00000 | 781.5 | NA | 1557.1 | 1452.2 | 9.42 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 1561.6 | 1458.4 | 9.46 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 1572.6 | 1466.2 | 9.52 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 1589.2 | 1482.5 | 9.63 | NA | . |
| HSat1A | experimental | D | 30 | 0 | NA | 0.99987 | 0.00000 | 268.7 | NA | 1660.4 | 1531.5 | 10.01 | NA | . |
| HSat1B | experimental | D | 1232 | 0 | NA | 0.97171 | 0.00000 | 507.6 | NA | 1979.7 | 1832.7 | 12.13 | NA | . |
| bSat | experimental | D | 408 | 0 | NA | 0.99947 | 0.00000 | 558.3 | NA | 2434.6 | 2085.2 | 14.88 | NA | . |
| HSat3 | experimental | D | 92 | 0 | NA | 0.99952 | 0.00000 | 1360.0 | NA | 3416.7 | 3099.5 | 21.59 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped`

Notes:

- DYZ3, DXZ1 (named subsets of aSatHOR's intervals) and aSatHOR are both fetched: their intervals are read once
- DYZ2 (a named subset of HSat1B's intervals) and HSat1B are both fetched: their intervals are read once
- DYZ1 (a named subset of HSat3's intervals) and HSat3 are both fetched: their intervals are read once
- experimental options selected (DYZ3, DXZ1, SST1, DYZ2, DYZ1, HSat2, aSatHOR, CER, SATR, ACRO, HSat1A, HSat1B, bSat, HSat3): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-13"></a>
## 13. core_tel with the lite control regions and a capture target of 0.995

```
ngsdose fetchplan --preset core_tel --controls resources/GRCh38/controls.lite200.bed --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**385.9 MB per genome (2.39% of the CRAM; 355.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 382 | 0 | NA | NA | NA | 120.9 | NA | 120.9 | 93.1 | 0.78 | NA | 382 regions of controls.lite200.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.9 | 0.0 | 124.0 | 94.9 | 0.78 | per byte | . |
| DJ | shipped | A | 48 | 11 | 0.995 | 0.99501 | 0.00218 | 44.5 | 6.0 | 179.1 | 137.6 | 1.05 | per byte | . |
| rDNA45S | shipped | A | 6 | 13 | 0.995 | 0.99821 | 0.00111 | 208.1 | 3.9 | 358.6 | 330.2 | 2.21 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| TEL | shipped | B | 55 | 8 | 0.995 | 0.99698 | 0.00092 | 119.6 | 15.2 | 385.9 | 355.5 | 2.39 | per byte | . |

controls.txt (the -c FASTA): controls.lite200.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.lite200.bed, not the menu's controls.bed: the fetch must pass -c controls.lite200.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names
- capture targets saved 25.3 MB of the plan: 411.3 MB with every option's intervals whole, 385.9 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (rDNA5S, DJ, rDNA45S, TEL) need not add up to the plan's
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, TEL, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans): scans not used to learn the sinks

<a id="example-14"></a>
## 14. biobank_lite at a capture target of 0.995

```
ngsdose fetchplan --preset biobank_lite --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**567.1 MB per genome (3.48% of the CRAM; 544.8 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.9 | 0.0 | 263.3 | 233.1 | 1.64 | per byte | . |
| DJ | shipped | A | 53 | 6 | 0.995 | 0.99580 | 0.00140 | 46.0 | 3.5 | 305.5 | 277.9 | 1.92 | per byte | 6 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for CER), so they cost nothing and add to its capture |
| rDNA45S | shipped | A | 6 | 13 | 0.995 | 0.99821 | 0.00111 | 208.1 | 3.9 | 496.5 | 475.1 | 3.03 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| TEL | shipped | B | 55 | 8 | 0.995 | 0.99698 | 0.00092 | 119.6 | 15.2 | 520.7 | 501.1 | 3.20 | per byte | . |
| SST1 | experimental | C | 28 | 6 | 0.995 | 0.99609 | 0.00346 | 13.7 | 3.3 | 534.4 | 514.8 | 3.28 | per read (fewer bytes than per byte) | . |
| SATR | experimental | D | 21 | 0 | 0.995 | 0.99898 | 0.00000 | 104.0 | 0.0 | 545.0 | 525.6 | 3.34 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| ACRO | experimental | D | 15 | 7 | 0.995 | 0.99528 | 0.00421 | 111.8 | 2.0 | 564.3 | 540.0 | 3.46 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| CER | experimental | D | 44 | 5 | 0.995 | 0.99500 | 0.00466 | 27.0 | 2.8 | 567.1 | 544.8 | 3.48 | per read (fewer bytes than per byte) | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for DJ), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,SST1,SATR,ACRO,CER`

Notes:

- capture targets saved 33.2 MB of the plan: 600.3 MB with every option's intervals whole, 567.1 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (rDNA5S, DJ, rDNA45S, TEL, SST1, SATR, ACRO, CER) need not add up to the plan's
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,SST1,SATR,ACRO,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, SATR, ACRO, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, TEL, SATR, ACRO, CER, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans), SST1 (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans): scans not used to learn the sinks

<a id="example-15"></a>
## 15. the four per-array options (preset xy_arrays)

```
ngsdose fetchplan --preset xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**340.7 MB per genome (2.05% of the CRAM; 293.9 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 5.0 | NA | 264.9 | 236.2 | 1.66 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.8 | NA | 274.7 | 244.8 | 1.73 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 307.4 | 276.7 | 1.86 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 340.7 | 293.9 | 2.05 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=aSatHOR,HSat1B,HSat3`

Notes:

- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- DYZ2 is a named subset of HSat1B's intervals, fetched without the rest of HSat1B's sinks: the fetch counts HSat1B only inside them, so `ngsdose estimate` reports DYZ2.mass_Mb and marks HSat1B itself as not measured (subset_only)
- DYZ1 is a named subset of HSat3's intervals, fetched without the rest of HSat3's sinks: the fetch counts HSat3 only inside them, so `ngsdose estimate` reports DYZ1.mass_Mb and marks HSat3 itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat2, bSat, SATR, ACRO, SST1, CER, not selected: count_flags.txt has --classes=aSatHOR,HSat1B,HSat3 (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, DYZ2, DYZ1): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-16"></a>
## 16. core_tel and the four per-array options

```
ngsdose fetchplan --preset core_tel xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**623.6 MB per genome (3.78% of the CRAM; 581.4 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 260.0 | NA | 260.0 | 231.4 | 1.63 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.9 | NA | 263.3 | 233.1 | 1.64 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 49.6 | NA | 309.6 | 281.1 | 1.94 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 503.9 | 482.5 | 3.07 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 137.7 | NA | 546.9 | 522.7 | 3.33 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 5.0 | NA | 551.9 | 525.7 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.8 | NA | 560.0 | 533.5 | 3.43 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 598.7 | 560.7 | 3.61 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 623.6 | 581.4 | 3.78 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,HSat1B,HSat3`

Notes:

- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- DYZ2 is a named subset of HSat1B's intervals, fetched without the rest of HSat1B's sinks: the fetch counts HSat1B only inside them, so `ngsdose estimate` reports DYZ2.mass_Mb and marks HSat1B itself as not measured (subset_only)
- DYZ1 is a named subset of HSat3's intervals, fetched without the rest of HSat3's sinks: the fetch counts HSat3 only inside them, so `ngsdose estimate` reports DYZ1.mass_Mb and marks HSat3 itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat2, bSat, SATR, ACRO, SST1, CER, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,HSat1B,HSat3 (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, DYZ2, DYZ1): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
