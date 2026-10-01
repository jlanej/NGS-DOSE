# fetchplan worked examples

Written by `resources/build/fetch_examples.sh` on 2026-10-01; do not edit by hand.
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
each whole CRAM. The engine joins a plan's intervals into one indexed fetch where less than 50,000 bp separate
them (16,384 in a BAM: `count --group-gap`), and each fetch decodes every slice that overlaps it (with its
container's compression header), so a slice under two fetches is decoded twice. In each table, mb_median is the
option's own intervals alone and cum_mb_median the plan up to and including that row, both priced so; cum_mb_floor
is the plan with every slice decoded once (the engine comes within a few percent of it: a slice is read twice only
where it spans two fetches more than 50,000 bp apart, as the few sparse chrY slices of a woman do). Over HTTPS
htslib moves more bytes than it decodes: each fetch opens a range request without an end, and what is in flight
at the next seek is thrown away (measured on the windows' candidates: 1,517 MB moved for 410 MB of slices); a
copy of the containers by exact byte range moves what is priced here. Candidate classes have no sinks yet: a plan
lists them for the whole-file scans only.

Engine for count_flags.txt: `ngs-dose` (ngs-dose 0.3.0; the binary's sha256 begins
`9cc6250da674f522`, as `--version` does not tell builds apart), which takes `count --classes`.
`count --classes` is in engines from the fetch-menu change of 2026-09-26 on; fae1124 lacks it. The engine changes only count_flags.txt, and whether fetchplan accepts a plan whose
panels define classes it does not select (`ngsdose fetchplan --help`, --engine).

| example | MB per genome | % of the CRAM | floor: every slice once |
| --- | ---: | ---: | ---: |
| [1. core: the bundle's positional classes](#example-1) | 485.4 | 2.97 | 482.5 |
| [2. core_tel: core and the telomeric repeat](#example-2) | 524.5 | 3.24 | 522.7 |
| [3. biobank_lite: core_tel and the four cheapest satellites](#example-3) | 580.3 | 3.59 | 573.9 |
| [4. core and the cheap satellites SST1, CER, SATR, ACRO at a capture target of 0.995](#example-4) | 436.0 | 2.66 | 415.1 |
| [5. core with the 200 lite control regions](#example-5) | 360.3 | 2.16 | 338.8 |
| [6. core_tel with a capture target for TEL alone](#example-6) | 512.3 | 3.15 | 508.5 |
| [7. the per-array options DXZ1 and DYZ3](#example-7) | 258.7 | 1.56 | 244.8 |
| [8. core_tel and the tier-A candidates, for the whole-file scans](#example-8) | 524.5 | 3.24 | 522.7 |
| [9. a budget of 600 MB, filled](#example-9) | 586.0 | 3.58 | 572.7 |
| [10. a budget of 1000 MB, filled](#example-10) | 914.3 | 5.65 | 893.5 |
| [11. the ten satellite families and the unmapped bin at a capture target of 0.995](#example-11) | 2773.6 | 17.31 | 2605.4 |
| [12. every option with sinks (each row's mb_median is that option alone)](#example-12) | 3351.6 | 21.06 | 3099.5 |
| [13. core_tel with the lite control regions and a capture target of 0.995](#example-13) | 377.9 | 2.25 | 355.5 |
| [14. biobank_lite at a capture target of 0.995](#example-14) | 552.2 | 3.40 | 544.9 |
| [15. the four per-array options (preset xy_arrays)](#example-15) | 314.6 | 1.87 | 293.9 |
| [16. core_tel and the four per-array options](#example-16) | 599.2 | 3.63 | 581.4 |
| [17. core with the karyotype set: every chromosome at full precision](#example-17) | 521.9 | 3.15 | 515.7 |
| [18. core with the screen set: every chromosome, the sex chromosomes at full precision](#example-18) | 370.0 | 2.21 | 358.0 |
| [19. core_tel with the karyotype set](#example-19) | 560.5 | 3.42 | 555.6 |
| [20. core with every region of the bundle (the controls and the windows)](#example-20) | 646.9 | 3.95 | 640.9 |

<a id="example-1"></a>
## 1. core: the bundle's positional classes

```
ngsdose fetchplan --preset core --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**485.4 MB per genome (2.97% of the CRAM; 482.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- (none)

<a id="example-2"></a>
## 2. core_tel: core and the telomeric repeat

```
ngsdose fetchplan --preset core_tel --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**524.5 MB per genome (3.24% of the CRAM; 522.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- (none)

<a id="example-3"></a>
## 3. biobank_lite: core_tel and the four cheapest satellites

```
ngsdose fetchplan --preset biobank_lite --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**580.3 MB per genome (3.59% of the CRAM; 573.9 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 542.4 | 538.7 | 3.35 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 548.7 | 544.0 | 3.39 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 562.9 | 554.1 | 3.47 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 580.3 | 573.9 | 3.59 | NA | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,SST1,CER,SATR,ACRO`

Notes:

- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,SST1,CER,SATR,ACRO (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, CER, SATR, ACRO): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-4"></a>
## 4. core and the cheap satellites SST1, CER, SATR, ACRO at a capture target of 0.995

```
ngsdose fetchplan --preset core --classes SST1 CER SATR ACRO --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**436.0 MB per genome (2.66% of the CRAM; 415.1 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.7 | 0.0 | 249.7 | 233.1 | 1.48 | per byte | . |
| DJ | shipped | A | 53 | 6 | 0.995 | 0.99528 | 0.00191 | 44.4 | 3.8 | 294.2 | 278.1 | 1.75 | per read (fewer bytes than per byte) | 6 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for CER), so they cost nothing and add to its capture |
| rDNA45S | shipped | A | 5 | 14 | 0.995 | 0.99816 | 0.00116 | 102.8 | 98.9 | 390.9 | 368.9 | 2.36 | per byte | . |
| SST1 | experimental | C | 28 | 6 | 0.995 | 0.99609 | 0.00346 | 13.5 | 3.5 | 404.9 | 382.4 | 2.44 | per read (fewer bytes than per byte) | . |
| SATR | experimental | D | 20 | 1 | 0.995 | 0.99720 | 0.00178 | 10.6 | 95.4 | 417.4 | 391.7 | 2.51 | per byte | . |
| ACRO | experimental | D | 14 | 8 | 0.995 | 0.99530 | 0.00419 | 20.0 | 97.1 | 431.6 | 410.0 | 2.63 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for rDNA45S), so they cost nothing and add to its capture |
| CER | experimental | D | 43 | 6 | 0.995 | 0.99500 | 0.00466 | 26.4 | 4.0 | 436.0 | 415.1 | 2.66 | per read (fewer bytes than per byte) | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,SST1,SATR,ACRO,CER`

Notes:

- capture targets saved 105.0 MB of the plan: 540.9 MB with every option's intervals whole, 436.0 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (rDNA5S, DJ, rDNA45S, SST1, SATR, ACRO, CER) need not add up to the plan's
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,SST1,SATR,ACRO,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, SATR, ACRO, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, SATR, ACRO, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), SST1 (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans): scans not used to learn the sinks

<a id="example-5"></a>
## 5. core with the 200 lite control regions

```
ngsdose fetchplan --preset core --controls resources/GRCh38/controls.lite200.bed --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**360.3 MB per genome (2.16% of the CRAM; 338.8 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 382 | 0 | NA | NA | NA | 112.9 | NA | 112.9 | 93.1 | 0.70 | NA | 382 regions of controls.lite200.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 115.2 | 94.9 | 0.71 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 159.3 | 142.5 | 0.99 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 360.3 | 338.8 | 2.16 | NA | . |

controls.txt (the -c FASTA): controls.lite200.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.lite200.bed, not the menu's controls.base.bed: the fetch must pass -c controls.lite200.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names

<a id="example-6"></a>
## 6. core_tel with a capture target for TEL alone

```
ngsdose fetchplan --preset core_tel --capture-class TEL=0.99 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**512.3 MB per genome (3.15% of the CRAM; 508.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 55 | 8 | 0.99 | 0.99698 | 0.00092 | 118.9 | 15.2 | 512.3 | 508.5 | 3.15 | per byte | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- capture targets saved 12.1 MB of the plan: 524.5 MB with every option's intervals whole, 512.3 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (TEL) need not add up to the plan's
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for TEL, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans): scans not used to learn the sinks

<a id="example-7"></a>
## 7. the per-array options DXZ1 and DYZ3

```
ngsdose fetchplan --classes DXZ1 DYZ3 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**258.7 MB per genome (1.56% of the CRAM; 244.8 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 250.8 | 236.2 | 1.49 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 258.7 | 244.8 | 1.56 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |

controls.txt (the -c FASTA): controls.base.fa.gz

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

**524.5 MB per genome (3.24% of the CRAM; 522.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |
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

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

scan_panels.txt (panels for the whole-file scans): panel.k31.tsv.gz, telomere.k31.panel.tsv.gz, macrosatellites.k31.panel.tsv.gz, multicopy-genes.k31.panel.tsv.gz, sex-chromosome-arrays.k31.panel.tsv.gz, rna-arrays.k31.panel.tsv.gz, nonhuman.k31.panel.tsv.gz, coding-vntrs.k31.panel.tsv.gz

Notes:

- (each of the 38 candidates has this note; the first:) DXZ4 is a candidate: no sinks have been learned for it, so a fetch cannot measure it. Its panel is in scan_panels.txt: load it in the whole-file scans, learn its sinks there (`ngsdose sinks SCANS`), check their capture on held-out scans (`ngsdose sinks HELD_OUT --evaluate BED`), then give it the sinks file in the menu

<a id="example-9"></a>
## 9. a budget of 600 MB, filled

```
ngsdose fetchplan --budget-mb 600 --fill --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**586.0 MB per genome (3.58% of the CRAM; 572.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 529.2 | 525.7 | 3.25 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 540.6 | 534.7 | 3.30 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 547.5 | 541.9 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 565.5 | 557.4 | 3.46 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 571.8 | 562.6 | 3.51 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 586.0 | 572.7 | 3.58 | NA | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--unmapped --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,CER,SATR`

Notes:

- 83 candidate option(s) are left out: they have no learned sinks yet
- DYZ2 (tier C) is left out: it would take the plan to 608.7 MB, past the budget of 600 MB
- DYZ1 (tier C) is left out: it would take the plan to 600.9 MB, past the budget of 600 MB
- HSat2 (tier C) is left out: it would take the plan to 844.9 MB, past the budget of 600 MB
- aSatHOR (tier C) is left out: it would take the plan to 1240.3 MB, past the budget of 600 MB
- ACRO (tier D) is left out: it would take the plan to 603.4 MB, past the budget of 600 MB
- HSat1A (tier D) is left out: it would take the plan to 751.3 MB, past the budget of 600 MB
- HSat1B (tier D) is left out: it would take the plan to 994.0 MB, past the budget of 600 MB
- bSat (tier D) is left out: it would take the plan to 1031.9 MB, past the budget of 600 MB
- HSat3 (tier D) is left out: it would take the plan to 1802.4 MB, past the budget of 600 MB
- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, ACRO, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,CER,SATR (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, SST1, CER, SATR): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-10"></a>
## 10. a budget of 1000 MB, filled

```
ngsdose fetchplan --budget-mb 1000 --fill --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**914.3 MB per genome (5.65% of the CRAM; 893.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 529.2 | 525.7 | 3.25 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 540.6 | 534.7 | 3.30 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 547.5 | 541.9 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 565.5 | 557.4 | 3.46 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 608.7 | 587.0 | 3.68 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 625.5 | 605.0 | 3.78 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |
| HSat2 | experimental | C | 53 | 0 | NA | 0.99993 | 0.00000 | 377.6 | NA | 879.5 | 858.3 | 5.44 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 885.1 | 863.5 | 5.48 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 899.2 | 873.6 | 5.55 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 914.3 | 893.5 | 5.65 | NA | . |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--unmapped --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,HSat1B,HSat3,HSat2,CER,SATR,ACRO`

Notes:

- 83 candidate option(s) are left out: they have no learned sinks yet
- aSatHOR (tier C) is left out: it would take the plan to 1565.6 MB, past the budget of 1000 MB
- HSat1A (tier D) is left out: it would take the plan to 1084.1 MB, past the budget of 1000 MB
- HSat1B (tier D) is left out: it would take the plan to 1265.6 MB, past the budget of 1000 MB
- bSat (tier D) is left out: it would take the plan to 1363.0 MB, past the budget of 1000 MB
- HSat3 (tier D) is left out: it would take the plan to 2072.1 MB, past the budget of 1000 MB
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

**2773.6 MB per genome (17.31% of the CRAM; 2605.4 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 259.7 | 242.7 | 1.53 | NA | every read without a coordinate (count --unmapped) |
| SST1 | experimental | C | 31 | 3 | 0.995 | 0.99522 | 0.00432 | 15.5 | 0.9 | 276.5 | 258.3 | 1.61 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat3), so they cost nothing and add to its capture |
| HSat2 | experimental | C | 44 | 9 | 0.995 | 0.99551 | 0.00442 | 370.9 | 5.4 | 648.5 | 612.9 | 3.88 | per byte | 6 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat1B, HSat3), so they cost nothing and add to its capture |
| aSatHOR | experimental | C | 209 | 39 | 0.995 | 0.99722 | 0.00156 | 790.0 | 14.6 | 1336.7 | 1235.0 | 8.17 | per byte | 5 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for HSat1A, bSat, HSat1B, HSat3), so they cost nothing and add to its capture |
| SATR | experimental | D | 21 | 0 | 0.995 | 0.99898 | 0.00000 | 104.0 | 0.0 | 1347.7 | 1246.0 | 8.23 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for HSat1B), so they cost nothing and add to its capture |
| ACRO | experimental | D | 16 | 6 | 0.995 | 0.99528 | 0.00421 | 112.1 | 3.2 | 1361.9 | 1259.0 | 8.33 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat1B), so they cost nothing and add to its capture |
| CER | experimental | D | 43 | 6 | 0.995 | 0.99500 | 0.00466 | 26.4 | 4.8 | 1383.7 | 1278.7 | 8.46 | per read (fewer bytes than per byte) | . |
| HSat1A | experimental | D | 17 | 13 | 0.995 | 0.99899 | 0.00088 | 220.8 | 45.4 | 1403.9 | 1296.9 | 8.56 | per byte | 7 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, bSat, HSat1B, HSat3), so they cost nothing and add to its capture |
| bSat | experimental | D | 300 | 108 | 0.995 | 0.99502 | 0.00444 | 465.9 | 79.5 | 1750.5 | 1549.8 | 10.50 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for aSatHOR, HSat1B), so they cost nothing and add to its capture |
| HSat1B | experimental | D | 1232 | 0 | 0.995 | 0.97171 | 0.00000 | 516.0 | 0.0 | 2081.5 | 1882.9 | 12.85 | per byte | . |
| HSat3 | experimental | D | 74 | 18 | 0.995 | 0.99505 | 0.00447 | 1009.2 | 347.2 | 2773.6 | 2605.4 | 17.31 | per byte | 12 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for HSat2, aSatHOR, HSat1B), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--unmapped`

Notes:

- capture targets saved 425.7 MB of the plan: 3199.3 MB with every option's intervals whole, 2773.6 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (SST1, HSat2, aSatHOR, SATR, ACRO, CER, HSat1A, bSat, HSat1B, HSat3) need not add up to the plan's
- HSat1B: all its intervals together capture 0.9717 (p10 over the scans of its statistics), below the target 0.995: all are kept
- experimental options selected (SST1, HSat2, aSatHOR, SATR, ACRO, CER, HSat1A, bSat, HSat1B, HSat3): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for SST1, HSat2, aSatHOR, SATR, ACRO, HSat1A, bSat, HSat1B, HSat3, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for SST1 (1648 scans), HSat2 (1648 scans), aSatHOR (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans), HSat1A (1648 scans), bSat (1648 scans), HSat1B (1648 scans), HSat3 (1648 scans): scans not used to learn the sinks

<a id="example-12"></a>
## 12. every option with sinks (each row's mb_median is that option alone)

```
ngsdose fetchplan --preset core_tel satellites xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**3351.6 MB per genome (21.06% of the CRAM; 3099.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 529.2 | 525.7 | 3.25 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 540.6 | 534.7 | 3.30 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 547.5 | 541.9 | 3.36 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 565.5 | 557.4 | 3.46 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 608.7 | 587.0 | 3.68 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 625.5 | 605.0 | 3.78 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |
| HSat2 | experimental | C | 53 | 0 | NA | 0.99993 | 0.00000 | 377.6 | NA | 879.5 | 858.3 | 5.44 | NA | . |
| aSatHOR | experimental | C | 248 | 0 | NA | 0.99878 | 0.00000 | 806.2 | NA | 1565.6 | 1452.2 | 9.49 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.2 | NA | 1571.8 | 1458.4 | 9.53 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 1583.2 | 1466.2 | 9.62 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 1600.8 | 1482.5 | 9.72 | NA | . |
| HSat1A | experimental | D | 30 | 0 | NA | 0.99987 | 0.00000 | 268.7 | NA | 1659.8 | 1531.5 | 10.06 | NA | . |
| HSat1B | experimental | D | 1232 | 0 | NA | 0.97171 | 0.00000 | 516.0 | NA | 1964.6 | 1832.7 | 12.19 | NA | . |
| bSat | experimental | D | 408 | 0 | NA | 0.99947 | 0.00000 | 547.0 | NA | 2377.9 | 2085.2 | 14.57 | NA | . |
| HSat3 | experimental | D | 92 | 0 | NA | 0.99952 | 0.00000 | 1352.2 | NA | 3351.6 | 3099.5 | 21.06 | NA | . |

controls.txt (the -c FASTA): controls.base.fa.gz

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

**377.9 MB per genome (2.25% of the CRAM; 355.5 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 382 | 0 | NA | NA | NA | 112.9 | NA | 112.9 | 93.1 | 0.70 | NA | 382 regions of controls.lite200.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.7 | 0.0 | 115.2 | 94.9 | 0.71 | per byte | . |
| DJ | shipped | A | 49 | 10 | 0.995 | 0.99515 | 0.00204 | 43.2 | 5.5 | 153.4 | 137.6 | 0.96 | per byte | . |
| rDNA45S | shipped | A | 6 | 13 | 0.995 | 0.99821 | 0.00111 | 208.1 | 3.9 | 351.2 | 330.2 | 2.10 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| TEL | shipped | B | 55 | 8 | 0.995 | 0.99698 | 0.00092 | 118.9 | 15.2 | 377.9 | 355.5 | 2.25 | per byte | . |

controls.txt (the -c FASTA): controls.lite200.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.lite200.bed, not the menu's controls.base.bed: the fetch must pass -c controls.lite200.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names
- capture targets saved 21.3 MB of the plan: 399.2 MB with every option's intervals whole, 377.9 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (rDNA5S, DJ, rDNA45S, TEL) need not add up to the plan's
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, TEL, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans): scans not used to learn the sinks

<a id="example-14"></a>
## 14. biobank_lite at a capture target of 0.995

```
ngsdose fetchplan --preset biobank_lite --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**552.2 MB per genome (3.40% of the CRAM; 544.9 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.7 | 0.0 | 249.7 | 233.1 | 1.48 | per byte | . |
| DJ | shipped | A | 53 | 6 | 0.995 | 0.99580 | 0.00140 | 44.3 | 3.5 | 294.1 | 277.9 | 1.75 | per byte | 6 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for CER), so they cost nothing and add to its capture |
| rDNA45S | shipped | A | 6 | 13 | 0.995 | 0.99821 | 0.00111 | 208.1 | 3.9 | 477.8 | 475.1 | 2.91 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| TEL | shipped | B | 55 | 8 | 0.995 | 0.99698 | 0.00092 | 118.9 | 15.2 | 504.2 | 501.1 | 3.10 | per byte | . |
| SST1 | experimental | C | 28 | 6 | 0.995 | 0.99609 | 0.00346 | 13.5 | 3.5 | 518.1 | 514.8 | 3.18 | per read (fewer bytes than per byte) | . |
| SATR | experimental | D | 21 | 0 | 0.995 | 0.99898 | 0.00000 | 104.0 | 0.0 | 531.1 | 525.6 | 3.25 | per byte | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| ACRO | experimental | D | 15 | 7 | 0.995 | 0.99552 | 0.00398 | 112.5 | 1.8 | 548.4 | 540.1 | 3.36 | per byte | 2 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for rDNA45S, TEL), so they cost nothing and add to its capture |
| CER | experimental | D | 44 | 5 | 0.995 | 0.99500 | 0.00466 | 27.0 | 3.2 | 552.2 | 544.9 | 3.40 | per read (fewer bytes than per byte) | 1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for DJ), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,SST1,SATR,ACRO,CER`

Notes:

- capture targets saved 28.1 MB of the plan: 580.3 MB with every option's intervals whole, 552.2 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: an interval the plan's other options fetch anyway is no saving, so the options' savings (rDNA5S, DJ, rDNA45S, TEL, SST1, SATR, ACRO, CER) need not add up to the plan's
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,SST1,SATR,ACRO,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, SATR, ACRO, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, TEL, SATR, ACRO, CER, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans), SST1 (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans): scans not used to learn the sinks

<a id="example-15"></a>
## 15. the four per-array options (preset xy_arrays)

```
ngsdose fetchplan --preset xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**314.6 MB per genome (1.87% of the CRAM; 293.9 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 250.8 | 236.2 | 1.49 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 258.7 | 244.8 | 1.56 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 292.9 | 276.7 | 1.75 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 314.6 | 293.9 | 1.87 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |

controls.txt (the -c FASTA): controls.base.fa.gz

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

**599.2 MB per genome (3.63% of the CRAM; 581.4 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 247.9 | NA | 247.9 | 231.4 | 1.48 | NA | 982 regions of controls.base.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 249.7 | 233.1 | 1.48 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 297.4 | 281.1 | 1.77 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 485.4 | 482.5 | 2.97 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 524.5 | 522.7 | 3.24 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 529.2 | 525.7 | 3.25 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 538.9 | 533.5 | 3.32 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 581.6 | 560.7 | 3.52 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 599.2 | 581.4 | 3.63 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |

controls.txt (the -c FASTA): controls.base.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,HSat1B,HSat3`

Notes:

- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- DYZ2 is a named subset of HSat1B's intervals, fetched without the rest of HSat1B's sinks: the fetch counts HSat1B only inside them, so `ngsdose estimate` reports DYZ2.mass_Mb and marks HSat1B itself as not measured (subset_only)
- DYZ1 is a named subset of HSat3's intervals, fetched without the rest of HSat3's sinks: the fetch counts HSat3 only inside them, so `ngsdose estimate` reports DYZ1.mass_Mb and marks HSat3 itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat2, bSat, SATR, ACRO, SST1, CER, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,HSat1B,HSat3 (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, DYZ2, DYZ1): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-17"></a>
## 17. core with the karyotype set: every chromosome at full precision

```
ngsdose fetchplan --preset core --controls karyotype --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**521.9 MB per genome (3.15% of the CRAM; 515.7 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 2647 | 0 | NA | NA | NA | 277.2 | NA | 277.2 | 266.4 | 1.67 | NA | 2647 regions of controls.karyotype.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 278.8 | 268.2 | 1.68 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 323.8 | 314.1 | 1.97 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 521.9 | 515.7 | 3.15 | NA | . |

controls.txt (the -c FASTA): p17.controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the controls FASTA of controls.karyotype.bed is cut from the bundle's controls.fa.gz when the plan is written (PREFIX.controls.fa.gz)
- the plan is costed on the control regions of controls.karyotype.bed, not the menu's controls.base.bed: the fetch must pass -c controls.karyotype.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names

<a id="example-18"></a>
## 18. core with the screen set: every chromosome, the sex chromosomes at full precision

```
ngsdose fetchplan --preset core --controls screen --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**370.0 MB per genome (2.21% of the CRAM; 358.0 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 977 | 0 | NA | NA | NA | 123.4 | NA | 123.4 | 114.6 | 0.74 | NA | 977 regions of controls.screen.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 125.2 | 116.4 | 0.75 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 172.9 | 164.6 | 1.04 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 370.0 | 358.0 | 2.21 | NA | . |

controls.txt (the -c FASTA): p18.controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the controls FASTA of controls.screen.bed is cut from the bundle's controls.fa.gz when the plan is written (PREFIX.controls.fa.gz)
- the plan is costed on the control regions of controls.screen.bed, not the menu's controls.base.bed: the fetch must pass -c controls.screen.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names

<a id="example-19"></a>
## 19. core_tel with the karyotype set

```
ngsdose fetchplan --preset core_tel --controls karyotype --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**560.5 MB per genome (3.42% of the CRAM; 555.6 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 2647 | 0 | NA | NA | NA | 277.2 | NA | 277.2 | 266.4 | 1.67 | NA | 2647 regions of controls.karyotype.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 278.8 | 268.2 | 1.68 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 323.8 | 314.1 | 1.97 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 521.9 | 515.7 | 3.15 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.7 | NA | 560.5 | 555.6 | 3.42 | NA | . |

controls.txt (the -c FASTA): p19.controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the controls FASTA of controls.karyotype.bed is cut from the bundle's controls.fa.gz when the plan is written (PREFIX.controls.fa.gz)
- the plan is costed on the control regions of controls.karyotype.bed, not the menu's controls.base.bed: the fetch must pass -c controls.karyotype.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names

<a id="example-20"></a>
## 20. core with every region of the bundle (the controls and the windows)

```
ngsdose fetchplan --preset core --controls all --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**646.9 MB per genome (3.95% of the CRAM; 640.9 MB with every slice decoded once).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_mb_floor | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 3247 | 0 | NA | NA | NA | 404.4 | NA | 404.4 | 395.1 | 2.49 | NA | 3247 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 406.0 | 396.9 | 2.50 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 47.7 | NA | 452.5 | 442.8 | 2.78 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 646.9 | 640.9 | 3.95 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.bed, not the menu's controls.base.bed: the fetch must pass -c controls.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names
