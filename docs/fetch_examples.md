# fetchplan worked examples

Written by `resources/build/fetch_examples.sh` on 2026-09-26; do not edit by hand.
The READMEs and `docs/DESIGN.md` quote plan figures (MB per genome, intervals kept, expected capture)
from this file: rerun the script after a change to the menu, a sinks or statistics file, the control
regions or `ngsdose fetchplan`. The descriptions in `resources/fetch_menu.tsv` and the headers and
table of `resources/experimental/subsets/` carry rounded costs written when those files were made;
examples 7 and 14 to 16 give those options' own costs (mb_median) to check them against.
Medians of cumulative totals are not additive: the difference between two rows, or between two plans'
totals, is not the median of what an option adds per CRAM.

Costs are the CRAM slices a fetch reads, in MB (1e6 bytes), median over 13 CRAM indexes:
HG00096, HG00706, HG01084, HG01552, HG01871, HG02383, HG02451, HG02466, HG02792, HG03007, HG03072, HG03267, NA12878
(`resources/build/fetch_examples.sh --get DIR` downloads them; sha256 of their concatenation
`a9f6c4d4766a7095`).
These are NYGC 30x CRAMs of the 1000 Genomes high-coverage release (bwa-mem, GRCh38 with decoys
and HLA; contigs from `GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`). The sinks, their statistics and so these costs belong to
that aligner and reference: another pipeline learns its own sinks from its own scans. Percentages are of
each whole CRAM. In each table, mb_median is the option's own intervals alone and cum_mb_median the plan up to
and including that row, each CRAM slice counted once. Candidate classes have no sinks yet: a plan lists them for
the whole-file scans only.

Engine for count_flags.txt: `ngs-dose` (ngs-dose 0.1.0; the binary's sha256 begins
`7775654435ab863d`, as `--version` does not tell builds apart), which takes `count --classes`.
`count --classes` is in engines from the fetch-menu change of 2026-09-26 on; fae1124 lacks it. The engine changes only count_flags.txt, and whether fetchplan accepts a plan whose
panels define classes it does not select (`ngsdose fetchplan --help`, --engine).

| example | MB per genome | % of the CRAM |
| --- | ---: | ---: |
| [1. core: the bundle's positional classes](#example-1) | 482.5 | 2.85 |
| [2. core_tel: core and the telomeric repeat](#example-2) | 522.7 | 3.09 |
| [3. biobank_lite: core_tel and the four cheapest satellites](#example-3) | 573.9 | 3.39 |
| [4. core and the cheap satellites SST1, CER, SATR, ACRO at a capture target of 0.995](#example-4) | 414.7 | 2.45 |
| [5. core with the 200 lite control regions](#example-5) | 338.8 | 1.99 |
| [6. core_tel with a capture target for TEL alone](#example-6) | 508.5 | 3.00 |
| [7. the per-array options DXZ1 and DYZ3](#example-7) | 244.8 | 1.47 |
| [8. core_tel and the tier-A candidates, for the whole-file scans](#example-8) | 522.7 | 3.09 |
| [9. a budget of 600 MB, filled](#example-9) | 594.0 | 3.58 |
| [10. a budget of 1000 MB, filled](#example-10) | 893.5 | 5.52 |
| [11. the ten satellite families and the unmapped bin at a capture target of 0.995](#example-11) | 2603.5 | 16.20 |
| [12. every option with sinks (each row's mb_median is that option alone)](#example-12) | 3099.5 | 19.45 |
| [13. core_tel with the lite control regions and a capture target of 0.995](#example-13) | 356.9 | 2.13 |
| [14. biobank_lite at a capture target of 0.995](#example-14) | 544.5 | 3.21 |
| [15. the four per-array options (preset xy_arrays)](#example-15) | 293.9 | 1.71 |
| [16. core_tel and the four per-array options](#example-16) | 581.4 | 3.43 |

<a id="example-1"></a>
## 1. core: the bundle's positional classes

```
ngsdose fetchplan --preset core --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**482.5 MB per genome (2.85% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- (none)

<a id="example-2"></a>
## 2. core_tel: core and the telomeric repeat

```
ngsdose fetchplan --preset core_tel --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**522.7 MB per genome (3.09% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- (none)

<a id="example-3"></a>
## 3. biobank_lite: core_tel and the four cheapest satellites

```
ngsdose fetchplan --preset biobank_lite --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**573.9 MB per genome (3.39% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 538.7 | 3.19 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.0 | NA | 544.0 | 3.23 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 554.1 | 3.29 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 573.9 | 3.39 | NA | . |

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

**414.7 MB per genome (2.45% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.7 | 0.0 | 233.1 | 1.40 | per byte | . |
| DJ | shipped | A | 53 | 6 | 0.995 | 0.99580 | 0.00140 | 43.6 | 3.2 | 277.9 | 1.66 | per byte | 9 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for CER), so they cost nothing and add to its capture |
| rDNA45S | shipped | A | 5 | 14 | 0.995 | 0.99816 | 0.00116 | 102.8 | 98.9 | 368.7 | 2.17 | per byte | . |
| SST1 | experimental | C | 28 | 6 | 0.995 | 0.99609 | 0.00346 | 13.5 | 3.3 | 382.3 | 2.25 | per read (fewer bytes than per byte) | . |
| SATR | experimental | D | 20 | 1 | 0.995 | 0.99720 | 0.00178 | 10.6 | 95.4 | 391.5 | 2.31 | per byte | . |
| ACRO | experimental | D | 14 | 8 | 0.995 | 0.99528 | 0.00421 | 19.2 | 97.2 | 409.6 | 2.42 | per byte | . |
| CER | experimental | D | 44 | 5 | 0.995 | 0.99500 | 0.00466 | 27.0 | 2.8 | 414.7 | 2.45 | per read (fewer bytes than per byte) | 1 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for DJ), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,SST1,SATR,ACRO,CER`

Notes:

- capture targets saved 118.5 MB of the plan: 533.2 MB with every option's intervals whole, 414.7 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: slices another option reads are no saving, so the options' savings (rDNA5S, DJ, rDNA45S, SST1, SATR, ACRO, CER) need not add up to the plan's
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,SST1,SATR,ACRO,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, SATR, ACRO, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, SATR, ACRO, CER, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), SST1 (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans): scans not used to learn the sinks

<a id="example-5"></a>
## 5. core with the 200 lite control regions

```
ngsdose fetchplan --preset core --controls resources/GRCh38/controls.lite200.bed --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**338.8 MB per genome (1.99% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 382 | 0 | NA | NA | NA | 93.1 | NA | 93.1 | 0.55 | NA | 382 regions of controls.lite200.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 94.9 | 0.57 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 142.5 | 0.85 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 338.8 | 1.99 | NA | . |

controls.txt (the -c FASTA): controls.lite200.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.lite200.bed, not the menu's controls.bed: the fetch must pass -c controls.lite200.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names

<a id="example-6"></a>
## 6. core_tel with a capture target for TEL alone

```
ngsdose fetchplan --preset core_tel --capture-class TEL=0.99 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**508.5 MB per genome (3.00% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 55 | 8 | 0.99 | 0.99698 | 0.00092 | 118.9 | 15.2 | 508.5 | 3.00 | per byte | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- capture targets saved 14.1 MB of the plan: 522.7 MB with every option's intervals whole, 508.5 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: slices another option reads are no saving, so the options' savings (TEL) need not add up to the plan's
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for TEL, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans): scans not used to learn the sinks

<a id="example-7"></a>
## 7. the per-array options DXZ1 and DYZ3

```
ngsdose fetchplan --classes DXZ1 DYZ3 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**244.8 MB per genome (1.47% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 236.2 | 1.40 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 244.8 | 1.47 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |

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

**522.7 MB per genome (3.09% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |
| DXZ4 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| CT47 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| RS447 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| MSR5p | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| FLJ40296 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| D4Z4 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| KIV2 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| C4 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| AMY1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| AMY2B | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SMN | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SMN1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| RHD | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HBA | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| GSTM1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| GSTT1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| UGT2B17 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| LCE3BC | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| APOBEC3B | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HPR | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| CCL3L | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| ORM1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| NPY4R | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| TSPY | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| DAZ | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| BPY2 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| DYZ19 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| OPN1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| TDNA1Q23 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| RNU1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SNORD116 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SNORD115 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HHV6A | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| HHV6B | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| SMRV | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| VNTR_ACAN | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| VNTR_MUC1 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |
| VNTR_MUC19 | candidate | A | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | scan only: no sinks learned yet (panel in scan_panels.txt) |

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

**594.0 MB per genome (3.58% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 525.7 | 3.10 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 534.7 | 3.15 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 541.9 | 3.24 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 557.4 | 3.34 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 587.0 | 3.54 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.0 | NA | 594.0 | 3.58 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,HSat1B,CER`

Notes:

- 83 candidate option(s) are left out: they have no learned sinks yet
- DYZ1 (tier C) is left out: it would take the plan to 605.0 MB, past the budget of 600 MB
- HSat2 (tier C) is left out: it would take the plan to 831.1 MB, past the budget of 600 MB
- aSatHOR (tier C) is left out: it would take the plan to 1199.2 MB, past the budget of 600 MB
- SATR (tier D) is left out: it would take the plan to 604.5 MB, past the budget of 600 MB
- ACRO (tier D) is left out: it would take the plan to 612.0 MB, past the budget of 600 MB
- HSat1A (tier D) is left out: it would take the plan to 735.2 MB, past the budget of 600 MB
- bSat (tier D) is left out: it would take the plan to 854.9 MB, past the budget of 600 MB
- HSat1B (tier D) is left out: it would take the plan to 949.1 MB, past the budget of 600 MB
- HSat3 (tier D) is left out: it would take the plan to 1773.1 MB, past the budget of 600 MB
- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- DYZ2 is a named subset of HSat1B's intervals, fetched without the rest of HSat1B's sinks: the fetch counts HSat1B only inside them, so `ngsdose estimate` reports DYZ2.mass_Mb and marks HSat1B itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat2, HSat3, bSat, SATR, ACRO, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,HSat1B,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, SST1, DYZ2, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-10"></a>
## 10. a budget of 1000 MB, filled

```
ngsdose fetchplan --budget-mb 1000 --fill --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**893.5 MB per genome (5.52% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 525.7 | 3.10 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 534.7 | 3.15 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 541.9 | 3.24 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 557.4 | 3.34 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 587.0 | 3.54 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 605.0 | 3.58 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |
| HSat2 | experimental | C | 53 | 0 | NA | 0.99993 | 0.00000 | 369.6 | NA | 858.3 | 5.33 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.0 | NA | 863.5 | 5.37 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 873.6 | 5.42 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 893.5 | 5.52 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,SST1,HSat1B,HSat3,HSat2,CER,SATR,ACRO`

Notes:

- 83 candidate option(s) are left out: they have no learned sinks yet
- aSatHOR (tier C) is left out: it would take the plan to 1452.2 MB, past the budget of 1000 MB
- HSat1A (tier D) is left out: it would take the plan to 1042.2 MB, past the budget of 1000 MB
- bSat (tier D) is left out: it would take the plan to 1154.0 MB, past the budget of 1000 MB
- HSat1B (tier D) is left out: it would take the plan to 1239.1 MB, past the budget of 1000 MB
- HSat3 (tier D) is left out: it would take the plan to 2025.3 MB, past the budget of 1000 MB
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

**2603.5 MB per genome (16.20% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 242.7 | 1.43 | NA | every read without a coordinate (count --unmapped) |
| SST1 | experimental | C | 30 | 4 | 0.995 | 0.99609 | 0.00346 | 15.0 | 0.7 | 257.3 | 1.53 | per read (fewer bytes than per byte) | 2 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for aSatHOR, HSat3), so they cost nothing and add to its capture |
| HSat2 | experimental | C | 45 | 8 | 0.995 | 0.99562 | 0.00432 | 364.5 | 4.1 | 612.8 | 3.77 | per byte | 7 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for aSatHOR, HSat1B, HSat3), so they cost nothing and add to its capture |
| aSatHOR | experimental | C | 209 | 39 | 0.995 | 0.99746 | 0.00132 | 720.7 | 10.3 | 1228.8 | 7.49 | per read (fewer bytes than per byte) | 8 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for HSat2, HSat1A, bSat, HSat1B, HSat3), so they cost nothing and add to its capture |
| SATR | experimental | D | 21 | 0 | 0.995 | 0.99898 | 0.00000 | 104.0 | 0.0 | 1239.8 | 7.56 | per byte | 1 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for HSat1B), so they cost nothing and add to its capture |
| ACRO | experimental | D | 16 | 6 | 0.995 | 0.99528 | 0.00421 | 112.1 | 2.3 | 1252.7 | 7.65 | per byte | 2 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for aSatHOR, HSat1B), so they cost nothing and add to its capture |
| CER | experimental | D | 43 | 6 | 0.995 | 0.99500 | 0.00466 | 26.1 | 3.6 | 1272.4 | 7.77 | per read (fewer bytes than per byte) | . |
| HSat1A | experimental | D | 19 | 11 | 0.995 | 0.99913 | 0.00074 | 221.4 | 39.0 | 1291.4 | 7.87 | per byte | 8 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for aSatHOR, bSat, HSat1B, HSat3), so they cost nothing and add to its capture |
| bSat | experimental | D | 371 | 37 | 0.995 | 0.99766 | 0.00181 | 362.2 | 14.9 | 1544.6 | 9.29 | per byte | 3 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for aSatHOR, HSat1B), so they cost nothing and add to its capture |
| HSat1B | experimental | D | 1232 | 0 | 0.995 | 0.97171 | 0.00000 | 488.4 | 0.0 | 1877.6 | 11.46 | per byte | . |
| HSat3 | experimental | D | 81 | 11 | 0.995 | 0.99544 | 0.00408 | 1008.7 | 343.6 | 2603.5 | 16.20 | per byte | 12 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for HSat2, aSatHOR, HSat1B), so they cost nothing and add to its capture |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped`

Notes:

- capture targets saved 352.3 MB of the plan: 2955.7 MB with every option's intervals whole, 2603.5 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: slices another option reads are no saving, so the options' savings (SST1, HSat2, aSatHOR, SATR, ACRO, CER, HSat1A, bSat, HSat1B, HSat3) need not add up to the plan's
- HSat1B: all its intervals together capture 0.9717 (p10 over the scans of its statistics), below the target 0.995: all are kept
- experimental options selected (SST1, HSat2, aSatHOR, SATR, ACRO, CER, HSat1A, bSat, HSat1B, HSat3): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for SST1, HSat2, aSatHOR, SATR, ACRO, HSat1A, bSat, HSat1B, HSat3, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for SST1 (1648 scans), HSat2 (1648 scans), aSatHOR (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans), HSat1A (1648 scans), bSat (1648 scans), HSat1B (1648 scans), HSat3 (1648 scans): scans not used to learn the sinks

<a id="example-12"></a>
## 12. every option with sinks (each row's mb_median is that option alone)

```
ngsdose fetchplan --preset core_tel satellites xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**3099.5 MB per genome (19.45% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 525.7 | 3.10 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| unmapped | shipped | C | * | 0 | NA | NA | NA | 8.5 | NA | 534.7 | 3.15 | NA | every read without a coordinate (count --unmapped) |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 541.9 | 3.24 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| SST1 | experimental | C | 34 | 0 | NA | 0.99955 | 0.00000 | 16.5 | NA | 557.4 | 3.34 | NA | . |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 587.0 | 3.54 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 605.0 | 3.58 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |
| HSat2 | experimental | C | 53 | 0 | NA | 0.99993 | 0.00000 | 369.6 | NA | 858.3 | 5.33 | NA | . |
| aSatHOR | experimental | C | 248 | 0 | NA | 0.99878 | 0.00000 | 731.8 | NA | 1452.2 | 8.88 | NA | . |
| CER | experimental | D | 49 | 0 | NA | 0.99967 | 0.00000 | 30.0 | NA | 1458.4 | 8.92 | NA | . |
| SATR | experimental | D | 21 | 0 | NA | 0.99898 | 0.00000 | 104.0 | NA | 1466.2 | 8.97 | NA | . |
| ACRO | experimental | D | 22 | 0 | NA | 0.99950 | 0.00000 | 114.3 | NA | 1482.5 | 9.07 | NA | . |
| HSat1A | experimental | D | 30 | 0 | NA | 0.99987 | 0.00000 | 264.1 | NA | 1531.5 | 9.41 | NA | . |
| bSat | experimental | D | 408 | 0 | NA | 0.99947 | 0.00000 | 378.9 | NA | 1779.3 | 11.04 | NA | . |
| HSat1B | experimental | D | 1232 | 0 | NA | 0.97171 | 0.00000 | 488.4 | NA | 2085.2 | 12.98 | NA | . |
| HSat3 | experimental | D | 92 | 0 | NA | 0.99952 | 0.00000 | 1348.7 | NA | 3099.5 | 19.45 | NA | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--unmapped`

Notes:

- DYZ3, DXZ1 (named subsets of aSatHOR's intervals) and aSatHOR are both fetched: their intervals are read once
- DYZ2 (a named subset of HSat1B's intervals) and HSat1B are both fetched: their intervals are read once
- DYZ1 (a named subset of HSat3's intervals) and HSat3 are both fetched: their intervals are read once
- experimental options selected (DYZ3, DXZ1, SST1, DYZ2, DYZ1, HSat2, aSatHOR, CER, SATR, ACRO, HSat1A, bSat, HSat1B, HSat3): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet

<a id="example-13"></a>
## 13. core_tel with the lite control regions and a capture target of 0.995

```
ngsdose fetchplan --preset core_tel --controls resources/GRCh38/controls.lite200.bed --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**356.9 MB per genome (2.13% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 382 | 0 | NA | NA | NA | 93.1 | NA | 93.1 | 0.55 | NA | 382 regions of controls.lite200.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.7 | 0.0 | 94.9 | 0.57 | per byte | . |
| DJ | shipped | A | 49 | 10 | 0.995 | 0.99533 | 0.00186 | 42.7 | 3.3 | 138.9 | 0.83 | per read (fewer bytes than per byte) | . |
| rDNA45S | shipped | A | 6 | 13 | 0.995 | 0.99821 | 0.00111 | 208.1 | 3.9 | 331.6 | 1.94 | per byte | 1 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| TEL | shipped | B | 55 | 8 | 0.995 | 0.99698 | 0.00092 | 118.9 | 15.2 | 356.9 | 2.13 | per byte | . |

controls.txt (the -c FASTA): controls.lite200.fa.gz

count_flags.txt: (empty: the fetch needs no count flags)

Notes:

- the plan is costed on the control regions of controls.lite200.bed, not the menu's controls.bed: the fetch must pass -c controls.lite200.fa.gz (written to controls.txt); a fetch with another controls file reads other regions than costed here, and `ngsdose estimate` accepts only the bundle's controls or a subset the bundle names
- capture targets saved 23.3 MB of the plan: 380.2 MB with every option's intervals whole, 356.9 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: slices another option reads are no saving, so the options' savings (rDNA5S, DJ, rDNA45S, TEL) need not add up to the plan's
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, rDNA45S, TEL, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans): scans not used to learn the sinks

<a id="example-14"></a>
## 14. biobank_lite at a capture target of 0.995

```
ngsdose fetchplan --preset biobank_lite --capture 0.995 --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**544.5 MB per genome (3.21% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | 0.995 | 0.99970 | 0.00000 | 1.7 | 0.0 | 233.1 | 1.40 | per byte | . |
| DJ | shipped | A | 53 | 6 | 0.995 | 0.99528 | 0.00191 | 43.7 | 3.4 | 278.1 | 1.66 | per read (fewer bytes than per byte) | 9 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for CER), so they cost nothing and add to its capture |
| rDNA45S | shipped | A | 6 | 13 | 0.995 | 0.99821 | 0.00111 | 208.1 | 3.9 | 474.8 | 2.80 | per byte | 1 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| TEL | shipped | B | 55 | 8 | 0.995 | 0.99698 | 0.00092 | 118.9 | 15.2 | 500.8 | 2.95 | per byte | . |
| SST1 | experimental | C | 28 | 6 | 0.995 | 0.99609 | 0.00346 | 13.5 | 3.3 | 514.5 | 3.03 | per read (fewer bytes than per byte) | . |
| SATR | experimental | D | 21 | 0 | 0.995 | 0.99898 | 0.00000 | 104.0 | 0.0 | 525.3 | 3.10 | per byte | 1 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| ACRO | experimental | D | 15 | 7 | 0.995 | 0.99528 | 0.00421 | 111.8 | 2.0 | 539.7 | 3.18 | per byte | 1 interval(s) its target alone would drop are kept: the plan reads their CRAM slices anyway (for TEL), so they cost nothing and add to its capture |
| CER | experimental | D | 43 | 6 | 0.995 | 0.99500 | 0.00466 | 26.1 | 3.7 | 544.5 | 3.21 | per read (fewer bytes than per byte) | . |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,SST1,SATR,ACRO,CER`

Notes:

- capture targets saved 29.5 MB of the plan: 573.9 MB with every option's intervals whole, 544.5 MB as planned (medians). An option's mb_saved is what trimming it saved with the rest of the plan as it is: slices another option reads are no saving, so the options' savings (rDNA5S, DJ, rDNA45S, TEL, SST1, SATR, ACRO, CER) need not add up to the plan's
- the panels loaded also define HSat1A, HSat1B, HSat2, HSat3, bSat, aSatHOR, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,SST1,SATR,ACRO,CER (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (SST1, SATR, ACRO, CER): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
- expected capture is the 10th percentile over the scans of its statistics of the capture of the intervals kept; for rDNA5S, DJ, rDNA45S, TEL, SATR, ACRO, kept by share per byte (order 'per byte') or keeping intervals the plan reads for other options, it is a lower bound: the larger of the capture of all the class's intervals less the largest share each dropped interval held in any one scan, and the capture of the longest run of the statistics' own order kept whole; held-out for rDNA5S (1375 scans), DJ (1375 scans), rDNA45S (1375 scans), TEL (1375 scans), SST1 (1648 scans), SATR (1648 scans), ACRO (1648 scans), CER (1648 scans): scans not used to learn the sinks

<a id="example-15"></a>
## 15. the four per-array options (preset xy_arrays)

```
ngsdose fetchplan --preset xy_arrays --crai CRAI... --contigs GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --engine ENGINE
```

**293.9 MB per genome (1.71% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 236.2 | 1.40 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 244.8 | 1.47 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 276.7 | 1.64 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 293.9 | 1.71 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |

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

**581.4 MB per genome (3.43% of the CRAM).**

| option | status | tier | intervals | dropped | capture_target | expected_capture | capture_lost | mb_median | mb_saved | cum_mb_median | cum_pct_median | order | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| controls | shipped | A | 982 | 0 | NA | NA | NA | 231.4 | NA | 231.4 | 1.39 | NA | 982 regions of controls.bed, padded by 600 bp as the engine reads them |
| rDNA5S | shipped | A | 2 | 0 | NA | 0.99970 | 0.00000 | 1.7 | NA | 233.1 | 1.40 | NA | . |
| DJ | shipped | A | 59 | 0 | NA | 0.99719 | 0.00000 | 46.6 | NA | 281.1 | 1.68 | NA | . |
| rDNA45S | shipped | A | 19 | 0 | NA | 0.99932 | 0.00000 | 211.5 | NA | 482.5 | 2.85 | NA | . |
| TEL | shipped | B | 63 | 0 | NA | 0.99790 | 0.00000 | 136.4 | NA | 522.7 | 3.09 | NA | . |
| DYZ3 | experimental | C | 4 | 0 | NA | NA | NA | 3.3 | NA | 525.7 | 3.10 | NA | a named subset of aSatHOR's learned intervals (DYZ3.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DYZ3.mass_Mb |
| DXZ1 | experimental | C | 4 | 0 | NA | NA | NA | 9.5 | NA | 533.5 | 3.18 | NA | a named subset of aSatHOR's learned intervals (DXZ1.bed), fetched whole and counted as aSatHOR: `ngsdose estimate` reports DXZ1.mass_Mb |
| DYZ2 | experimental | C | 78 | 0 | NA | NA | NA | 34.3 | NA | 560.7 | 3.40 | NA | a named subset of HSat1B's learned intervals (DYZ2.bed), fetched whole and counted as HSat1B: `ngsdose estimate` reports DYZ2.mass_Mb |
| DYZ1 | experimental | C | 1 | 0 | NA | NA | NA | 39.1 | NA | 581.4 | 3.43 | NA | a named subset of HSat3's learned intervals (DYZ1.bed), fetched whole and counted as HSat3: `ngsdose estimate` reports DYZ1.mass_Mb |

controls.txt (the -c FASTA): controls.fa.gz

count_flags.txt: `--classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,HSat1B,HSat3`

Notes:

- DYZ3, DXZ1 are named subsets of aSatHOR's intervals, fetched without the rest of aSatHOR's sinks: the fetch counts aSatHOR only inside them, so `ngsdose estimate` reports DYZ3.mass_Mb, DXZ1.mass_Mb and marks aSatHOR itself as not measured (subset_only)
- DYZ2 is a named subset of HSat1B's intervals, fetched without the rest of HSat1B's sinks: the fetch counts HSat1B only inside them, so `ngsdose estimate` reports DYZ2.mass_Mb and marks HSat1B itself as not measured (subset_only)
- DYZ1 is a named subset of HSat3's intervals, fetched without the rest of HSat3's sinks: the fetch counts HSat3 only inside them, so `ngsdose estimate` reports DYZ1.mass_Mb and marks HSat3 itself as not measured (subset_only)
- the panels loaded also define HSat1A, HSat2, bSat, SATR, ACRO, SST1, CER, not selected: count_flags.txt has --classes=rDNA5S,DJ,rDNA45S,TEL,aSatHOR,HSat1B,HSat3 (ENGINE count --help lists it), so the fetch counts and reports only the selected classes and reads only their sinks. An engine without --classes (fae1124, 7772e32) refuses the flag: make the plan with --engine naming the engine that runs the fetch
- experimental options selected (DYZ3, DXZ1, DYZ2, DYZ1): their sinks come from scans, but no fetch has been compared with a scan of the same genome yet
