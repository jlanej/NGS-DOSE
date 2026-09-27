# Candidate panels: scan first, then learn the sinks

These are 83 new classes in seven k-mer panels, for **whole-file scans**. None of them has sinks
yet, so no fetch can measure them. Where a class's reads land is learned from the scans that carry its panel
(`ngsdose sinks`), checked on scans the sinks were not learned from, and only then used by a fetch.
The coordinates the builds used (a class's copies in GRCh38 and CHM13) only keep a class's own copies from
deleting its k-mers. They are not where the reads are assumed to be.

| panel file | classes | k-mers | what |
| --- | --- | --- | --- |
| `macrosatellites.k31.panel.tsv.gz` | 12 | 45,362 | tandem macrosatellites and RNA-gene arrays (DXZ4, CT47, RS447, MSR5p, FLJ40296, RNU2, D4Z4, ZAV, REXO1L1) and the three D4Z4 distal haplotypes |
| `multicopy-genes.k31.panel.tsv.gz` | 27 | 544,098 | multi-copy and multi-allelic genes (LPA KIV-2, C4, SMN, AMY1, RHD, the common deletions, ...) |
| `sex-chromosome-arrays.k31.panel.tsv.gz` | 13 | 126,346 | chrY ampliconic genes and DYZ19; the chrX opsin, GAGE, CT45 and SPANXB arrays |
| `rna-arrays.k31.panel.tsv.gz` | 5 | 148,026 | the 1q23 tRNA-gene array, U1 and U3 snRNA genes, SNORD116 and SNORD115 |
| `nonhuman.k31.panel.tsv.gz` | 6 | 395,211 | HHV-6A, HHV-6B, HHV-7, SMRV, phiX, EBV type 2 |
| `nonhuman-myco.k31.panel.tsv.gz` | 1 | 3,731,523 | culture Mycoplasma (five species, one class); a file of its own because it is 74% of the candidate k-mers and about 340 MB of a scan's memory |
| `coding-vntrs.k31.panel.tsv.gz` | 19 | 86,162 | long-unit coding VNTRs (ACAN, MUC1, MUC19, FLG, NEB, ...), one compositional class each |

`candidates.tsv` has one row per class: its panel, kind, unit, k-mers, recall, tier, what it measures, the
truth available for 1000 Genomes samples, the reference copies, and the build's caveats. `units/` holds the
sequence each class's k-mer positions refer to. `ngsdose estimate` reads positional classes from here
(`resources.ExperimentalUnits`), matching the panel by the sha256 the counts file records.
`reference_copies.GRCh38.bed` gives the reference copies widened to the engine's placement bins (1 kb for a
positional class, 10 kb for a compositional one), for comparison with learned sinks. `recipes/` and `resources/build/build_candidate_panels.sh` rebuild the panels. Every class has the
row status `candidate` in `resources/fetch_menu.tsv`.

## Why scan first

A fetch reads the control regions and each loaded class's sink intervals, and nothing else. A read the
aligner put somewhere else is lost without any sign. Where the reads go depends on the aligner, the
reference (alt, decoy and unplaced contigs) and the person, and it can only be seen in a whole-file scan,
which classifies every read wherever it is placed, the unmapped ones included.

One whole-file scan of NA12878 (NYGC 30x CRAM, female) with every candidate panel loaded shows why this
matters (`na12878_at_reference_copies` in `candidates.tsv`: the share of a class's reads placed in
`reference_copies.GRCh38.bed`, scored by `ngsdose sinks --evaluate`, which counts a placement bin only when all of it
lies inside). Of the 72 classes with at least 25
reads, 58 had at least 98% of their reads at their reference copies. The other 14 did not:

- `GSTT1` had 65%. GRCh38 has the gene only on `chr22_KI270879v1_alt`, and ALT-aware bwa moves reads off it.
- `D4Z4_4qA` had 81%, `CCL3L` 85%, `D4Z4_10q` 91%, `SPANXB` 92%, `RHD` 93%, `VNTR_MUC1` 95%, `GSTM1` 96%,
  `D4Z4_4qB` 96%, `VNTR_NCAPG2` 97% and `APOBEC3B` 97.6%. For `CCL3L`, `SPANXB` and `RHD` the rest are reads
  of paralog alleles the reference does not show (below).
- `UGT2B17` had 0%. NA12878 carries the deletion (estimated copy number 0.04), and its 174 reads are
  off-target hits: single-substitution alleles at other loci (below).
- In a female, `DAZ` got 109 reads and `RBMY` 47, stray reads outside chrY. These set the floor for the
  chrY classes.
- `HHV7` had 283 reads at chromosome ends in the first assembly, and `HHV6A` 21, from simple-repeat k-mers
  since removed (below); now both have 0.

This is one genome's observation, not a sink. In NA12878 most classes' reads sat at their reference copies,
but some classes did not, the share will differ between people (deletions, alleles the aligner places
elsewhere) and between pipelines, and only scans of many genomes show by how much. Nothing here sets a
sink.

## Loading them in scans

```bash
C=resources/experimental/candidates
ngs-dose count -m scan -i sample.cram -T ref.fa -c resources/GRCh38/controls.fa.gz \
    -p resources/GRCh38/panel.k31.tsv.gz -p resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz \
    -p resources/experimental/telomere.k31.panel.tsv.gz \
    -p $C/macrosatellites.k31.panel.tsv.gz -p $C/multicopy-genes.k31.panel.tsv.gz \
    -p $C/sex-chromosome-arrays.k31.panel.tsv.gz -p $C/rna-arrays.k31.panel.tsv.gz \
    -p $C/nonhuman.k31.panel.tsv.gz -p $C/nonhuman-myco.k31.panel.tsv.gz -p $C/coding-vntrs.k31.panel.tsv.gz \
    -o sample.scan.json.gz
```

Any subset of these files can be loaded. No k-mer is shared between two of them or with a shipped panel,
so a class has the same k-mers whichever other panels a run loads. To choose by class, tier or group, use
`ngsdose fetchplan --preset core candidates -o PREFIX` (or `candidates_A`, or a group: `cand_macrosatellites`,
`cand_multicopy_genes`, `cand_sex_arrays`, `cand_rna_arrays`, `cand_nonhuman`, `cand_coding_vntrs`; or
`--classes` naming them). It writes every selected panel, candidates included, to `PREFIX.scan_panels.txt`
for the scans, and the fetch it plans leaves candidates out. A plan must fetch at least one class with
sinks, so the candidates go with a fetchable preset such as `core` or `core_tel` (example 8 of
[docs/fetch_examples.md](../../../docs/fetch_examples.md)); candidates alone are refused and nothing is
written. In NGS-DOSE-1000G the files go into `EXTRA_PANELS` (`pipeline/config.sh`).

Loading all seven changes no count of any shipped class (checked in `tests/test_candidate_panels.py`).
It does change the panel list and `panel_sha256` that the counts files record, so if they are loaded part
way through a cohort, the cohort tools see mixed panels (`ngsdose sinks --allow-mixed-panels`, and the
cohort check).

**What a scan pays.** On the whole NA12878 CRAM (768.6 M records, 8 threads, one run each):

| engine | time | peak memory |
| --- | --- | --- |
| fetch-menu build | 98.9 s without the candidates, 107.0-122.9 s with them | 1.35 GB without, 1.48-1.71 GB with |
| fae1124 | 100.3-106.1 s without, 109.6-119.5 s with | 1.30-1.43 GB without, 1.49-1.81 GB with |

The runs were repeated after each round of review changes; the spread is run-to-run variation on a shared
machine. With the final panels (sha256 below), fae1124 took 106.1 s and 1.43 GB without the candidates and
119.5 s and 1.65 GB with them, and the fetch-menu build 122.9 s and 1.55 GB with them.

The counts file grows from 237 kB to 329 kB (gzip), or from 1.69 MB to 2.29 MB of JSON.

On the test fixture (198,848 records, 4 threads, median of 3 runs), time goes from 1.2 s to 3.7-3.8 s and
peak memory from 370 MB to 800 MB with all three engines (fae1124, 7772e32 and the fetch-menu build). Almost all of
that is loading 5.08 M k-mers. The fixture holds almost no reads of these classes, so no per-read cost
shows there. The fixture is 1/3,865 of the 30x file, and extrapolating from it gives about +2.5 s and
+430 MB per genome, a fixed cost. The whole-file measurement with the final panels adds about 13 s with fae1124
(106.1 to 119.5 s), so about 11 s, or 10%, is per read.
Without `nonhuman-myco`, the fixture scan takes 1.8 s and 447 MB (+0.6 s and +72 MB).

## Learning and checking the sinks

1. Scan a training set with the panels loaded, then learn the sinks. Positional classes are learned by
   default; compositional ones (`VNTR_*`, `DYZ19`, `MYCO`, `EBV2`) must be named.
   ```bash
   ngsdose sinks TRAIN/*.scan.json.gz --classes VNTR_ACAN VNTR_MUC1 ... DYZ19 MYCO EBV2 -o sinks.candidates.bed
   ```
   The rule is the bundle's. The engine records placements in 1-kb bins (10-kb bins for a compositional
   class). In each training scan, a bin that holds reads of the class is kept when its 10-kb neighbourhood
   holds at least 25 of them and at least 1e-5 of the class's reads in that scan; the bins kept in any
   training scan are merged and padded by 1 kb.
2. Check capture on held-out scans, and write the statistics that `fetchplan --capture` trims by:
   ```bash
   ngsdose sinks HELD_OUT/*.scan.json.gz --evaluate sinks.candidates.bed --held-out --stats sinks.candidates.stats.tsv
   ngsdose sinks HELD_OUT/*.scan.json.gz --evaluate resources/experimental/candidates/reference_copies.GRCh38.bed
   ```
   `--held-out` records in the statistics that these scans were not used to learn the sinks (`# held-out:
   yes`), and the file carries each interval's largest share in any one scan (`share_max`), which
   `fetchplan --crai` needs to trim per byte of the CRAM slices. The second command measures how much each
   class would lose if its reference coordinates were fetched instead.
3. Look at what is not captured. Reads fully unmapped (the viruses, MYCO, part of any class) are fetched
   only with `--unmapped`, a separate menu option. Classes with a floor of off-target reads (`UGT2B17`,
   `GSTM1`, `APOBEC3B`, the chrY classes in females) will show sinks that hold noise, which the class's `hit_frac` or a
   per-sex view can separate.
4. Then change the class's row in `resources/fetch_menu.tsv`: status `experimental`, the sinks BED and its
   statistics. `ngsdose fetchplan` then costs it and can fetch it.

A class whose capture stays low stays scan-only, or, like HSat1B among the satellites (96.6% in the worst
held-out scan), is fetched with that caveat recorded in its menu row.

## What they would cost to fetch

These are planning numbers from one genome, not sinks: until sinks are learned from the cohort's scans, a
candidate is scan-only and `fetchplan` does not cost it. The fetch cost here is for the intervals that
`ngsdose sinks` learns from that one scan of NA12878 (the rule above). It is counted in bytes beyond
today's fetch (the controls and `core_tel`, 522.7 MB), from the CRAM slices of 13 NYGC bwa-mem indexes
(median; `ngsdose.cost`); another aligner would place the reads, and so price them, differently.
The chrY classes have no bins in a female, and the viruses other than EBV2 had too few reads to place.

| selection | classes | MB beyond today's fetch | % of a 16.5-GB CRAM |
| --- | --- | --- | --- |
| all candidates | 83 | 108.9 | 0.66 |
| tier A | 38 | 33.0 | 0.20 |
| tiers A and B | 58 | 57.1 | 0.35 |
| all but `D4Z4`, its three distal classes and `TBC1D3` | 78 | 60.6 | 0.37 |
| macrosatellites / multicopy-genes / coding-vntrs | 12 / 27 / 19 | 49.7 / 37.0 / 7.7 | |
| rna-arrays / sex-chromosome-arrays (chrX only) / nonhuman (EBV2 only) | 5 / 13 / 7 | 5.8 / 4.1 / 2.2 | |
| the unmapped bin (the viruses, MYCO) | | 8.5 | 0.05 |

The most expensive are `D4Z4_10q` (21 MB), `TBC1D3` (10 MB), `D4Z4_4qA` (10 MB), `D4Z4` (9 MB), `RASA4`,
`REXO1L1` and `NOTCH2NL` (about 5 MB each). The CRAM slices at 4q35 and 10q26 are dense with dispersed
beta-satellite reads. Most other classes cost 0.2-1.5 MB each. For comparison, what today's options
cost is in [`docs/fetch_examples.md`](../../../docs/fetch_examples.md). Once sinks are learned, any subset these files carry is cheap to
fetch; the larger price is scanning enough genomes, with the panels loaded, to learn and check the sinks.
The group builds' own estimates, which assume the reads sit at the reference copies (macrosatellites
40 MB, multicopy-genes 33, coding-vntrs 14, sex-chromosome-arrays 10.5 with chrY, rna-arrays 5.5), agree in
size.

## The classes

Recall is the share of 150-bp reads cut from the class's unit (every 1/1500 of its length, alternately
reverse-complemented) that a scan assigns to the class, with every shipped and candidate panel loaded (gate
(e) below). The builds' own recall on each reference copy, other haplotypes and paralogs is in
`candidates.tsv`. The next three columns come from the NA12878 scan: reads (duplicates included), the share
at the reference copies (`-` under 25 reads), and the fetch cost described above. Tier is the group's
rating: A per-sample truth for 1000 Genomes samples; B partial or indirect truth; C a relative or technical
measure; D exploratory: the three D4Z4 distal classes and `SPANXB`, kept below 50% recall on the argument in
their caveats, and the coding VNTRs with no known phenotype.

| class | panel | kind | tier | k-mers | recall | NA12878 reads | at reference copies | fetch MB | measures |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `DXZ4` | macrosatellites | pos | A | 2,762 | 1.00 | 42,786 | 99.99% | 0.78 | DXZ4 macrosatellite units (Xq23; 3.0 kb) |
| `CT47` | macrosatellites | pos | A | 4,227 | 1.00 | 8,732 | 99.99% | 0.31 | CT47 cancer-testis gene-array units (Xq24; 4.9 kb) |
| `RS447` | macrosatellites | pos | A | 3,161 | 1.00 | 53,507 | 99.96% | 1.54 | RS447/USP17L megasatellite units (4p16; 4.7 kb) |
| `MSR5p` | macrosatellites | pos | A | 2,936 | 1.00 | 50,877 | 99.98% | 1.11 | MSR5p/TAF11-like macrosatellite units (5p15; 3.4 kb) |
| `FLJ40296` | macrosatellites | pos | A | 6,438 | 1.00 | 12,921 | 99.67% | 0.95 | FLJ40296/PRR20 array units (13q21; 6.6 kb) |
| `RNU2` | macrosatellites | pos | B | 5,115 | 0.99 | 56,474 | 99.94% | 1.40 | U2 snRNA gene units (17q21, next to BRCA1; 6.1 kb) |
| `D4Z4` | macrosatellites | pos | A | 1,708 | 0.92 | 42,481 | 99.34% | 9.36 | D4Z4 units, 4q35 + 10q26 together (the FSHD macrosatellite) |
| `ZAV` | macrosatellites | pos | B | 5,262 | 1.00 | 16,114 | 99.86% | 0.85 | ZAV/MSat10 units (9q32; 5.4 kb) |
| `REXO1L1` | macrosatellites | pos | C | 10,330 | 1.00 | 259,658 | 99.91% | 4.68 | REXO1L1 gene-array units (8q21; 12.2 kb), relative |
| `D4Z4_4qA` | macrosatellites | pos | D | 1,411 | 0.33 | 406 | 81.28% | 10.13 | number of 4qA-type distal ends (FSHD-permissive haplotype) |
| `D4Z4_4qB` | macrosatellites | pos | D | 1,263 | 0.20 | 571 | 96.32% | 0.58 | number of 4qB-type distal ends |
| `D4Z4_10q` | macrosatellites | pos | D | 749 | 0.21 | 643 | 90.51% | 21.37 | number of 10qA-type distal ends |
| `KIV2` | multicopy-genes | pos | A | 3,234 | 0.69 | 18,033 | 99.99% | 0.95 | LPA kringle IV-2 repeats (Lp(a)) |
| `C4` | multicopy-genes | pos | A | 14,163 | 0.69 | 6,274 | 99.71% | 0.94 | C4 genes (C4A + C4B, long and short) |
| `HERVC4` | multicopy-genes | pos | B | 5,185 | 1.00 | 1,585 | 100.00% | 0.90 | long C4 genes (the HERV-K(C4) insertion) |
| `CYP21A2` | multicopy-genes | pos | B | 744 | 0.52 | 522 | 99.62% | 0.46 | CYP21A2 gene copies (not the pseudogene) |
| `AMY1` | multicopy-genes | pos | A | 6,514 | 0.99 | 6,449 | 99.94% | 0.65 | salivary amylase AMY1 genes |
| `AMY2B` | multicopy-genes | pos | A | 20,689 | 0.98 | 5,822 | 99.62% | 0.34 | AMY2B copies |
| `SMN` | multicopy-genes | pos | A | 20,775 | 0.97 | 13,677 | 99.52% | 1.92 | SMN1 + SMN2 |
| `SMN1` | multicopy-genes | pos | A | 332 | 0.06 | 330 | 99.70% | 0.75 | SMN1 copies (SMA carriers), from the SMN1/SMN2 sites |
| `RHD` | multicopy-genes | pos | A | 16,311 | 0.67 | 5,350 | 92.67% | 0.94 | RHD copies (RhD-negative = 0) |
| `HBA` | multicopy-genes | pos | A | 1,861 | 1.00 | 1,178 | 100.00% | 0.26 | alpha-globin genes (HBA1 + HBA2) |
| `GSTM1` | multicopy-genes | pos | A | 3,430 | 0.92 | 763 | 95.81% | 0.19 | GSTM1 copies (null deletion) |
| `GSTT1` | multicopy-genes | pos | A | 6,720 | 0.95 | 2,197 | 65.41% | 0.40 | GSTT1 copies (null deletion) |
| `UGT2B17` | multicopy-genes | pos | A | 21,945 | 0.92 | 174 | 0.00% | 0.29 | UGT2B17 copies (116-kb deletion) |
| `LCE3BC` | multicopy-genes | pos | A | 29,459 | 0.99 | 3,987 | 98.14% | 0.53 | LCE3B/LCE3C segment (psoriasis-risk deletion) |
| `APOBEC3B` | multicopy-genes | pos | A | 21,227 | 0.94 | 7,767 | 97.59% | 0.67 | APOBEC3B deletion segment |
| `HPR` | multicopy-genes | pos | A | 11,287 | 0.97 | 3,252 | 99.81% | 0.32 | haptoglobin-related gene copies |
| `HPdup` | multicopy-genes | pos | B | 990 | 0.84 | 758 | 100.00% | 0.16 | HP exon 3-4 duplication (Hp1/Hp2 alleles) |
| `CCL3L` | multicopy-genes | pos | A | 16,464 | 0.93 | 3,495 | 85.09% | 1.43 | CCL3L-CCL4L block copies |
| `ORM1` | multicopy-genes | pos | A | 1,373 | 0.69 | 715 | 99.86% | 0.17 | ORM1 copies |
| `DEFB` | multicopy-genes | pos | B | 89,172 | 0.99 | 47,048 | 99.67% | 2.42 | 8p23.1 beta-defensin repeat |
| `FCGR3` | multicopy-genes | pos | B | 7,459 | 1.00 | 3,758 | 99.84% | 0.43 | FCGR3A + FCGR3B copies |
| `DEFA1A3` | multicopy-genes | pos | B | 14,944 | 0.89 | 8,308 | 99.11% | 0.55 | DEFA1/DEFA3 tandem units |
| `NPY4R` | multicopy-genes | pos | A | 109,570 | 1.00 | 45,815 | 99.77% | 2.08 | 10q11.22 NPY4R segmental duplication |
| `SULT1A` | multicopy-genes | pos | B | 12,724 | 0.76 | 5,149 | 99.46% | 0.56 | SULT1A2-SULT1A1 duplicon pair |
| `RASA4` | multicopy-genes | pos | B | 31,096 | 0.73 | 62,953 | 98.92% | 5.53 | 7q22.1 RASA4 tandem unit |
| `NOTCH2NL` | multicopy-genes | pos | B | 73,895 | 0.99 | 94,688 | 99.89% | 4.57 | NOTCH2NL-type copies (1q21.1) |
| `TBC1D3` | multicopy-genes | pos | C | 2,535 | 0.76 | 13,454 | 98.99% | 10.43 | 17q12 TBC1D3 paralogs |
| `TSPY` | sex-chromosome-arrays | pos | A | 16,690 | 0.97 | 12 | - | - | TSPY array units (+ TSPY2) |
| `RBMY` | sex-chromosome-arrays | pos | B | 8,661 | 0.96 | 47 | 0.00% | - | RBMY1 genes |
| `DAZ` | sex-chromosome-arrays | pos | A | 23,072 | 0.47 | 109 | 0.92% | 0.27 | DAZ genes (gr/gr, b2/b3 deletions) |
| `BPY2` | sex-chromosome-arrays | pos | A | 20,189 | 1.00 | 1 | - | - | BPY2 genes |
| `CDY1` | sex-chromosome-arrays | pos | B | 802 | 0.80 | 0 | - | - | CDY1 genes |
| `CDY2` | sex-chromosome-arrays | pos | B | 837 | 0.80 | 0 | - | - | CDY2 genes |
| `DYZ19` | sex-chromosome-arrays | comp | A | 3,223 | 0.95 | 3 | - | - | DYZ19 array mass (Yq11) |
| `OPN1` | sex-chromosome-arrays | pos | A | 30,816 | 0.96 | 27,048 | 99.60% | 0.84 | opsin array units (OPN1LW + OPN1MW) |
| `OPN1LW` | sex-chromosome-arrays | pos | C | 161 | 0.54 | 119 | 100.00% | 0.26 | LW-type opsin exon 5 |
| `OPN1MW` | sex-chromosome-arrays | pos | C | 130 | 0.54 | 202 | 100.00% | 0.63 | MW-type opsin exon 5 |
| `GAGE` | sex-chromosome-arrays | pos | B | 7,590 | 0.97 | 46,263 | 99.89% | 1.10 | GAGE array units (Xp11) |
| `CT45` | sex-chromosome-arrays | pos | B | 10,976 | 0.92 | 21,434 | 99.75% | 0.84 | CT45 array units (Xq26) |
| `SPANXB` | sex-chromosome-arrays | pos | D | 3,199 | 0.51 | 2,410 | 92.07% | 1.01 | SPANXB copies |
| `TDNA1Q23` | rna-arrays | pos | A | 6,304 | 1.00 | 40,377 | 99.98% | 1.54 | 1q23.3 tRNA-gene array units |
| `RNU1` | rna-arrays | pos | A | 400 | 0.22 | 5,972 | 99.93% | 2.11 | canonical U1 snRNA units (1p36) |
| `SNORD3` | rna-arrays | pos | B | 2,649 | 0.90 | 4,900 | 99.92% | 1.19 | U3 snoRNA modules (17p11.2 REPA/REPB) |
| `SNORD116` | rna-arrays | pos | A | 54,013 | 0.98 | 13,890 | 99.70% | 0.34 | SNORD116 cluster dosage profile (15q11.2) |
| `SNORD115` | rna-arrays | pos | A | 84,660 | 0.94 | 25,517 | 99.90% | 0.61 | SNORD115 cluster dosage profile (15q11.2) |
| `HHV6A` | nonhuman | pos | A | 115,912 | 0.98 | 0 | - | - | HHV-6A (iciHHV-6A; whole genome vs solo DR) |
| `HHV6B` | nonhuman | pos | A | 119,799 | 0.99 | 0 | - | - | HHV-6B (iciHHV-6B) |
| `HHV7` | nonhuman | pos | C | 137,144 | 0.98 | 0 | - | - | HHV-7 load (and a competitor for HHV-6 DR reads) |
| `SMRV` | nonhuman | pos | A | 7,903 | 0.93 | 0 | - | - | squirrel monkey retrovirus (LCL culture) |
| `PHIX` | nonhuman | pos | C | 5,386 | 1.00 | 10 | - | - | phiX spike-in (index hopping) |
| `MYCO` | nonhuman-myco | comp | C | 3,731,523 | 1.00 | 0 | - | - | culture Mycoplasma (five species) |
| `EBV2` | nonhuman | comp | B | 9,067 | 0.98 | 34 | 100.00% | 2.23 | EBV type 2 |
| `VNTR_ACAN` | coding-vntrs | comp | A | 626 | 1.00 | 535 | 99.63% | 0.19 | ACAN exon-12 VNTR length (height) |
| `VNTR_MUC1` | coding-vntrs | comp | A | 1,366 | 1.00 | 920 | 95.33% | 0.39 | MUC1 VNTR length (urea, gout) |
| `VNTR_MUC19` | coding-vntrs | comp | A | 11,795 | 1.00 | 3,651 | 99.92% | 0.33 | MUC19 VNTR length (introgressed haplotype) |
| `VNTR_FLG` | coding-vntrs | comp | B | 7,106 | 1.00 | 2,702 | 100.00% | 0.25 | filaggrin repeats |
| `VNTR_NEB` | coding-vntrs | comp | B | 12,744 | 0.99 | 7,851 | 99.66% | 0.41 | nebulin triplicate blocks |
| `VNTR_MUC4` | coding-vntrs | comp | C | 7,272 | 1.00 | 4,052 | 99.63% | 0.86 | MUC4 VNTR length |
| `VNTR_DMBT1` | coding-vntrs | comp | C | 4,886 | 0.99 | 4,121 | 99.59% | 0.33 | DMBT1 SRCR-SID repeats |
| `VNTR_HRNR` | coding-vntrs | comp | C | 4,783 | 1.00 | 4,268 | 99.93% | 0.40 | hornerin repeat length |
| `VNTR_MUC12` | coding-vntrs | comp | D | 3,905 | 0.97 | 7,852 | 99.96% | 0.47 | MUC12 repeat length |
| `VNTR_MUC6` | coding-vntrs | comp | D | 7,812 | 1.00 | 2,673 | 99.74% | 1.02 | MUC6 repeat length |
| `VNTR_AHNAK2` | coding-vntrs | comp | D | 9,226 | 1.00 | 3,041 | 99.93% | 0.23 | AHNAK2 repeat length |
| `VNTR_EPPK1` | coding-vntrs | comp | D | 2,268 | 0.97 | 3,158 | 100.00% | 0.47 | epiplakin repeat length |
| `VNTR_C2orf78` | coding-vntrs | comp | D | 7,004 | 1.00 | 8,389 | 99.53% | 0.52 | C2orf78 repeat length |
| `VNTR_PLIN4` | coding-vntrs | comp | D | 2,617 | 1.00 | 946 | 99.89% | 0.43 | PLIN4 repeat length |
| `VNTR_TRIOBP` | coding-vntrs | comp | D | 805 | 1.00 | 397 | 98.74% | 0.40 | TRIOBP exon-7 repeat length |
| `VNTR_NCAPG2` | coding-vntrs | comp | D | 484 | 1.00 | 596 | 96.81% | 0.37 | NCAPG2 repeat length |
| `VNTR_ZNF512B` | coding-vntrs | comp | D | 292 | 1.00 | 164 | 99.39% | 0.23 | ZNF512B repeat length |
| `VNTR_TUBGCP2` | coding-vntrs | comp | D | 438 | 1.00 | 113 | 100.00% | 0.21 | TUBGCP2 repeat length |
| `VNTR_UVSSA` | coding-vntrs | comp | D | 733 | 1.00 | 302 | 100.00% | 0.21 | UVSSA repeat length |

A first estimate from that one scan, by `ngsdose.estimate` with the units here (positional classes; not
validated), gives:

- `C4` 3.0, `HERVC4` 2.0, `AMY1` 5.9, `SMN` 4.0, `RHD` 1.1, `HBA` 4.1, `GSTM1` 1.1, `GSTT1` 2.2 and
  `UGT2B17` 0.04 copies;
- `LCE3BC` 1.0, `APOBEC3B` 2.0, `DEFB` 4.0, `NOTCH2NL` 9.8, `KIV2` 38, `RNU1` 18.3 and `SNORD3` 9.4
  (reference 10);
- `ORM1` 1.8 (with the trimmed unit);
- the chrX arrays of a female: `OPN1` 5.8 units, with 2.1 LW-type (`OPN1LW`) and 4.1 MW-type (`OPN1MW`)
  exon-5 copies, and `DXZ4` 97;
- `REXO1L1` 165 units (relative; `cn_se_rel` 0.012), and `RNU2` 68;
- about 0 for every chrY class.

`ngsdose estimate` gives these as it gives any positional candidate whose unit is in `units/`: status
`experimental`, the all-window single-sample estimate, with `cn_se_rel`.

### Paralog pairs, and paralog alleles the reference does not show

Two reviews of the first build found these leaks, all fixed or documented before these panels go into a scan.

- **One mask for every class.** The sex-array build masked the loci of all its classes in one BED, so a
  k-mer that one class shares with a paralog copy inside another class's masked locus was never background.
  `OPN1LW` kept 31 k-mers (unit 778-808, in the intron after exon 5) that also occur in OPN1MW2 and OPN1MW3,
  and counted their reads (NA12878: 39 reads at MW2/MW3 besides 121 at LW); `CDY1` kept 20 k-mers (unit
  2720-2750) that also occur at both CDY2 loci, about 30 bp outside the `CDY2` unit, and took about 6.5% of
  each CDY2 copy's reads. `CDY1`, `CDY2`, `OPN1LW` and `OPN1MW` are now each taken from a build that masks
  only that class's own copies (`recipes/sex-chromosome-arrays/masks/*.own_copies.bed`), as the
  macrosatellite recipe already did for the D4Z4 classes; `OPN1LW` and `OPN1MW` also keep only k-mers found
  in every copy of their type (62 of `OPN1MW`'s were in OPN1MW but not MW2/MW3). Reads tiled over each
  whole LW and MW gene +-500 bp, and over the whole CDY loci, now give 0 reads to the other class of the
  pair in GRCh38 and CHM13. In NA12878, `OPN1LW` has 102 reads, all at the LW gene, and `OPN1MW` 194, all at
  the MW genes. `tests/test_candidate_panels.py` checks the four classes against the reference when it is
  present (`NGSDOSE_REF_DIR`).
- **`ORM1`.** Its first unit (gene +-1 kb) had 355 k-mers in the upstream flank that the ORM2 reference lacks
  but ORM2-locus alleles carry: NA12878 had 59 reads with `ORM1` k-mers placed at ORM2, 58 from that flank,
  raising those windows by about 30%. The unit is now the gene and 1 kb downstream (1,373 k-mers); NA12878
  has 644 reads at ORM1 and 1 at ORM2.
- **`RHD`.** Gene-converted RHCE alleles carry RHD-derived sequence (RHCE*C carries RHD exon 2). In NA12878,
  278 of the 4,704 reads with `RHD` k-mers around the two genes (6%, 238 of them at MAPQ >= 30) were placed at
  RHCE, spread along the unit. That cannot be removed by building against references, so it is a caveat:
  `RHD` windows read high in carriers of such alleles.
- **`CCL3L`.** No `CCL3L` k-mer occurs at the CCL3/CCL4 block (94.5% identical) in GRCh38 or CHM13, and reads
  tiled from those references give it nothing. Yet in NA12878, 491 of its 3,495 reads (14%) were placed
  at the CCL3/CCL4 block (chr17:36.08-36.10 Mb), 398 of the 433 reads there with >= 4 class k-mers at MAPQ
  >= 30, with a median of 17 class k-mers against 78 at the CCL3L3 block. They cluster at unit offsets 0-6,
  12, 16 and 20-24 kb, the gene-conversion pattern of `RHD`: this person's CCL3/CCL4 alleles carry CCL3L-like
  sequence. It adds about 0.3 copies at a copy number of 2 and depends on the CCL3/CCL4 haplotype, so
  those windows read high in carriers.
- **`SPANXB`.** 191 of its 2,410 reads (7.9%) were placed at the paralogs: 78 at SPANXC (chrX:141.23-141.25 Mb)
  and about 89 at SPANXA1/A2 (chrX:141.57-141.59 Mb), where a SPANXB-like segment (chrX:141,576,021-141,598,394)
  is left unmasked by design. The class reads SPANXB-type sequence, not SPANXB copies alone.

The general point: a class's leakage measured against the reference sequences (0.000 for every paralog
pair above) says nothing about paralog alleles that the references do not carry. The per-window
reliability of the paralog classes has to be learned from the cohort's scans before integer calls are
made. Among the other paralog pairs, NA12878's reads were clean: `SMN1` 1 read at SMN2 (293 at SMN1),
`CYP21A2` 1 at CYP21A1P (448), `AMY1` 2 at AMY2A/2B (6,112 at AMY1), `RS447` 10 at the chr8p23 paralogs
(53,507 in the scan), `KIV2` none outside the KIV-2 array.

### Off-target floors: one-substitution alleles and simple repeats

A group build drops every k-mer found in the background assemblies, but not k-mers one substitution away
from them. A person's allele at a background locus can then match a run of class k-mers: each such read
hits only 4-5 contiguous class k-mers, one allele's worth. In NA12878:

- `UGT2B17`: 174 reads in a genome with the deletion, none at the gene (35 at chr11:48.53 Mb, 16 at
  chr18:59.36 Mb, 14 at UGT2B15 on chr4:68.66 Mb, 13 at chr6:66.26 Mb, 9 at chr4:115.95 Mb): about
  0.05 copies.
- `GSTM1`: 32 of 763 reads (4%) away from the gene, 1-2 per locus: about 0.06 copies.
- `APOBEC3B`: 187 of 7,767 (2.4%), 70 of them among the other APOBEC3 genes at chr22:38.99-39.09 Mb: about
  0.07 copies.
- `DAZ`: 108 of 109 reads on other chromosomes in this woman, 33 at chr8:64.10 Mb: about 0.03 copies. `RBMY`:
  47, all off chrY (14 at chrX:1.37 Mb, 11 at chr1:147.18 Mb).

These matter for the 0/1/2 deletion calls (`UGT2B17`, `GSTM1`, `APOBEC3B`) and for the chrY classes in
women, and the floor depends on each person's alleles at those loci. The rna-arrays build removes such
k-mers with a one-substitution background filter (`recipes/rna-arrays/scripts/neighbors.py`); the other
groups do not, because for the paralog-specific classes (`SMN1`, `CYP21A2`, `RHD`, `OPN1LW`/`OPN1MW`,
`CDY1`/`CDY2`, the D4Z4 ends, `SPANXB`) the diagnostic k-mers are by construction one substitution from a
paralog, and the filter would remove them. The floors are measured instead: female scans give the chrY
floors, and the deletion classes' 0/1 thresholds should come from the cohort's scans.

A second kind of floor came from simple repeats. A k-mer of a telomere-repeat variant (TTAGGG, TCAGGG,
TTGGGG runs) or of a microsatellite can be absent from GRCh38 and CHM13, survive the background filter, and
still be in every genome's reads. The first assembly kept 27 such k-mers in `HHV7`'s direct repeat, 18 in
`HHV6A` (among them a (TATG)n run) and 10 in `HHV6B`; NA12878, which carries neither virus, got 283 `HHV7`
reads at chromosome ends (175 at chr5:10-13 kb) and 21 `HHV6A` reads. `recipes/assemble.py` now removes
from every panel each k-mer that is periodic with a period of 1-6 bp at <= 2 mismatches: 55 viral k-mers
and 23 in human classes (10 `DXZ4`, 4 `D4Z4` poly-C, 6 at the `D4Z4_4qB` end's telomere junction, 1 each in
`VNTR_EPPK1`, `DEFB` and `RASA4`). The same NA12878 scan then gives 0 `HHV7` and 0 `HHV6A` reads; the other
classes lost 0-3 reads, and the tiled-read recall changed by at most 0.001.

A third floor scales with a virus's load. `EBV2` is built from the two blocks where EBV type 2 differs most
from type 1, with six type-1 genomes as background, so no type-1 genome shares a k-mer with it. But a
type-1 read with a sequencing error or a minor variant at a site where the types differ by one base can
carry a run of `EBV2` k-mers. NA12878's LCL was transformed with B95-8 (type 1): 34 of its 1,777,780 chrEBV
records (27 not duplicate-flagged), 1.9e-5, carry >= 4 `EBV2` k-mers (4-31 each), all at chrEBV 79,941-87,429
(the EBNA-3 block) and none in the EBNA-2 block. In the 1,748 counted scans chrEBV records range from 126,194
to 30,872,944 (median 1,089,570), so at that rate the floor is about 21 reads in the median genome, 96 at the
95th percentile and up to about 590. It is proportional to the chrEBV load, as a real type-2 reading is too,
so type 2 is called from the `EBV2`/chrEBV ratio well above 2e-5, or from reads in the EBNA-2 block, not from
a count threshold. The test fixture's 4,621 chrEBV reads give 0 `EBV2` reads.

### Units that are one period, and a class's own copies among the decoys

A third review found three kinds of fault, all fixed in the files listed under "The final files".

- **`REXO1L1` was not one period.** Its first unit, chr8:85,772,425-85,784,618, was the first 12,194 bp of an
  atypical 15,275-bp unit that holds a 2.9-kb segment twice and lacked the last 3,080 bp of a normal unit.
  Its circular join made 30 k-mers that occur in neither assembly (1 and 3 were found), the duplicated
  segment's 2,711 k-mers were dropped as repeated within the unit, and recall was 0.61. The unit is now
  chr8:85,760,231-85,772,424: the unit's first 31-mer recurs exactly 12,194 bp on, all 30 junction 31-mers
  occur in the GRCh38 and CHM13 arrays, and minimap2 aligns it end to end to 62 CHM13 units at 99.3%
  identity (its self-alignment shows no internal duplication, only a 2.4-kb stretch tandem with a 138-bp period). It keeps 10,330
  k-mers, recall is 0.98-0.99 on both arrays and 0.999 on reads cut from the unit, and NA12878's estimate
  moves from 150 to 165 units with `cn_se_rel` 0.055 to 0.012 (the old unit's windows over the duplicated
  segment expected every unit to carry it).
- **Circular units whose last bases repeated their first.** Every circular unit was checked against the
  array it came from: the unit's first 31-mer must recur at the base after the unit, and the 30 31-mers
  across the circular join must occur in the array. `TSPY` was 21 bp too long (20,329 bp, the period at the
  array start, where this unit's next copy starts after 20,308), `OPN1` 22 bp (37,160 to 37,138), `CT45`
  16 bp (17,268 to 17,252; its units are 17,252-17,287 bp apart) and `DEFA1A3` 1 bp (it ended on the next
  unit's first base). `RNU2`, the GenBank clone U57614.1, is a HindIII fragment whose last 6 bp repeat its
  first 6 (AAGCTT); the unit is now bases 1-6,126, and 17 of its 30 junction 31-mers occur in both arrays
  (the others cover a small difference of this clone 17-19 bp before its end). All are trimmed; the other
  circular units (`DXZ4`, `CT47`, `RS447`, `MSR5p`, `FLJ40296`, `D4Z4`, `ZAV`, `KIV2`, `GAGE`, `TDNA1Q23`, and
  `PHIX`, a complete circular genome) were already exact. Each trimmed class lost 1-22 artificial k-mers.
  In NA12878, reads across the join now carry the real junction k-mers: `DEFA1A3` gained 17 reads and
  `RNU2` 1, the others none, and no estimate moved by more than 0.1 (`DEFA1A3` 3.93 to 3.94, `RNU2` 67.9 to
  68.0). `tests/test_candidate_panels.py` checks
  that no circular unit ends with its first bases and, with the reference present, that each unit cut from
  an assembly is followed there by the next unit.
- **Own copies among the hs38d1 decoys.** Counting every unit's 31-mers in all 3,341 non-primary contigs of
  the analysis set found two decoys that are pieces of a coding VNTR: chrUn_JTFH01000840v1_decoy (1,107 bp,
  1,081/1,095 identical to the CHM13 MUC19 array) and chrUn_JTFH01000899v1_decoy (1,055 bp, 1,047/1,055 to
  the CHM13 MUC6 array). They were not masked, so the 925 `VNTR_MUC19` and 916 `VNTR_MUC6` k-mers they share
  with the classes were dropped as background (the caveats had put the loss down to background elsewhere and,
  for MUC6, the MUC2/MUC5 cluster). Both are now masked as the class's own copies
  (`recipes/coding-vntrs/sources/loci.GRCh38_decoys.tsv`) and listed in `reference_copies.GRCh38.bed`, since
  reads may be placed on them. CHM13 recall rose from 0.964 to 1.000 (`VNTR_MUC19`) and from 0.937 to 1.000
  (`VNTR_MUC6`). Every other non-primary contig that holds a whole copy of a class, or is mostly made of
  one, was already masked (the eleven DYZ19 decoys, the chr4 alt with the inverted D4Z4 copy, the alt
  haplotypes of C4, CYP21A2, CCL3L, TBC1D3, DEFB, GSTT1, MUC4 and MUC6, and the RNU1 and DYZ19 randoms). A twelfth DYZ19-like decoy,
  chrUn_JTFH01001423v1_decoy, is 89% identical to the array but only 513 of its 1,406 distinct 31-mers (36%)
  are in it, against 89-100% for the eleven, so it stays background. Three unplaced contigs
  (chrUn_KI270363v1, KI270373v1, KI270508v1) are 93-96% identical to parts of the `DAZ` unit and also stay
  background.

### Known overlaps, by design

- **`C4` and `HERVC4`:** a read over the HERV-K(C4) insertion of a long C4 gene counts as `HERVC4`; the
  rest of the gene counts as `C4`.
- **`SMN` and `SMN1`:** a read over an SMN1-specific site counts for both, because positional classes are
  scored independently.
- **`OPN1`, `OPN1LW` and `OPN1MW`:** `OPN1`'s panel has no k-mer in OPN1MW exon 5, so exon-5 reads count
  as `OPN1LW` or `OPN1MW`. In the simulated reads, 15 reads of the real `OPN1` unit went to `OPN1MW`.
- **The D4Z4 classes:** the D4Z4 unit class and the three distal classes are built with only their own
  copies masked, which brought cross-talk on the reference copies to 0.
- **`HHV6A` and `HHV6B`:** the two species are 90% identical in U, so a carrier's reads hit both. The
  species is called from the ratio, not from presence.
- **The shipped families:** reads from beta-satellite pieces (`D4Z4_4qA`, `BPY2`) and from telomere-like
  repeats (the HHV direct repeats) also count as `bSat` or `TEL`, as they do today.

## Gates

**(a) k-mers.** No candidate k-mer is in a shipped panel (bundle, satellites, telomere), and none is periodic
with a period of 1-6 bp at <= 2 mismatches (see "Off-target floors"). No k-mer is in two candidate panels: `recipes/assemble.py` would remove such a k-mer from both, and it found none
(`shared_between_panels.tsv`, header only), because a group's k-mers occur nowhere in GRCh38 or CHM13
outside its own class copies. A separate count of the final files, recomputing each k-mer's canonical form,
found 5,076,728 candidate k-mers, no duplicate row, no k-mer in two of the ten files (seven candidate, three
shipped) and every `kmers_kept` equal to its class's rows. The 97 class names are unique, well under the
engine's limit of 255.

**(b) Shipped counts.** The fixture was scanned with the shipped panels alone and with all seven candidate
files added, three runs each, with fae1124 (the engine counting the cohort), main 7772e32 and the fetch-menu build.
In every case the classes[] entry of each of the 14 shipped classes, their 2,041 placement rows, the
controls, the regions, the GC tables and the contig tallies were identical. Only `below_threshold_reads`
moved (2,111 to 2,115). The same comparison on the whole NA12878 CRAM (768.6 M records, fae1124) was also
identical: 16,392 shipped placement rows, `ambiguous_reads` 3,214 both ways, `below_threshold_reads`
1,019,674 to 1,040,620. This branch, scanning the same CRAM with every candidate, gave the same entries and
placements for every class as fae1124. All of this was rerun on the final files.

**(c) Loading.** Each file loads on its own next to the shipped panels in all three engines, without a
shared-k-mer message or a k-mer-count warning, and changes no shipped count.

**(d) Cost.** See "What a scan pays" above.

**(e) Recognition.** Reads cut from every unit (98,846 reads, 139-1,498 per class) were added to the fixture
as reads without a coordinate. A model of the engine's rule (`src/count.rs` classify, applied over every
loaded panel) predicted each class's count, and the engine's count matched it exactly for all 83 classes in
all three engines. Shipped counts on this read set were again identical with and without the candidates.
Every class got its own reads. The only reads a class passed to another class are those listed under
"Known overlaps". `tests/test_candidate_panels.py` repeats (a)-(c), and (e) with five reads per class; set
`NGSDOSE_EXTRA_BINS=path/to/engine` to run it with the cohort's engine as well.

## What was left out, and why

**Macrosatellites**

- The 4qA beta-satellite block is part of the shipped `bSat` class.
- Arm-specific 4q and 10q D4Z4 unit classes are left out: the units are 98-99% identical, so as two classes
  their shared k-mers would drop from both.
- The inverted D4Z4 copy (D4S2463) is not single-copy, so it is not a control region.
- RNU2-2 and the U2 pseudogenes are single-copy loci that share the snRNA k-mers with the array.

**Multicopy genes**

- `AMY2A` keeps 1,046 k-mers and has a recall of 0.22-0.24, because of the AMY2Ap partial copies.
- Olduvai CON1 has no primary consensus reachable here.
- C4A versus C4B would rest on about 16 reads per copy, and gene conversion moves the isotype.
- The FCGR3A/FCGR3B and DEFB REPD/REPP splits would be biased, so the sums are measured instead.
- The RASA4 Y731C variant was not requested.

**Sex-chromosome arrays**

- DYZ1, DYZ2, DYZ17/18 and cenX/cenY are already in the shipped `HSat3`, `HSat1B` and `aSatHOR` classes;
  they are analyses of those classes' placements.
- A TSPY class without the TSPY2 unit has a recall of 0.52.
- Whole-gene OPN1LW/OPN1MW typing is left out: outside exon 5 the differences are mostly private to single
  units.
- A SPANX family class would mix the variable SPANXB with invariant copies.
- The DAZ internal repeats would need a class of their own.

**RNA arrays**

- The variant U1 genes (vU1) are not the canonical array.
- The per-isodecoder tRNA classes have a recall of 0.81, too few reads per copy.
- The U3 gene alone is contained in `SNORD3`.
- A compositional REPA/REPB class would compete with the satellite families for reads.

**Non-human**

- The anelloviruses are too diverse for exact 31-mers: mean recall 0.017.
- *Acholeplasma* is in the background instead of being a class.
- Per-species Mycoplasma classes would each take a small share of an off-panel species' reads.
- HHV-6 DR and U are one unit profile, not two classes.
- The chrEBV wild-type window is a region, not a class.
- HIV-1 and HBV have a recall of 0.2-0.6 on other subtypes with one reference each.

**Coding VNTRs**

- LPA KIV-2 is `KIV2` in the multicopy-genes panel.
- The short-unit or short-array VNTRs (TENT5A, TCHH, TMCO1, GP1BA, SBSN) vary below what depth can resolve.
- EIF3H's k-mers occur elsewhere in GRCh38.
- The segmental-duplication families lose their k-mers to their paralogs.
- The loci below Mukamel's IBD2R 0.7 and 1-kb cuts, and the non-coding VNTRs, are left for another group.

The group summaries give each reason in full.

## Rebuilding

`resources/build/build_candidate_panels.sh` runs each group's recipe (`recipes/<group>/build.sh`, with its
scripts and small inputs) into `work/candidates/<group>/`. It then runs `recipes/assemble.py` and copies
the panels and units here.

- **What the recipes use:** the engine's panel builder, the GRCh38 analysis set and CHM13v2.0 in
  `work/ref`, and the shipped panels (every recipe removes their k-mers).
- **Downloads:** GenBank records from NCBI E-utilities for RNU2, D4Z4 and the non-human group.
- **Assembly:** `assemble.py` checks that no candidate k-mer is in a shipped panel, removes k-mers shared
  between candidate panels and simple-repeat k-mers (period 1-6 bp, <= 2 mismatches), splits MYCO into a file of its own, and writes gzip with no timestamp, so that
  a rebuild gives the same bytes and the same sha256. A panel taken from here rather than rebuilt (a group
  left out of `CANDIDATE_GROUPS`) keeps the removal counts of its `##filter` line, so re-assembling it gives
  the same bytes.
- **Options:** `CANDIDATE_GROUPS` rebuilds some groups only; `RECALL=1` also measures each group's recall
  (more downloads; the HPRC haplotype rows of the macrosatellite recall need local extracts);
  `DERIVE_MASKS=1` (multicopy-genes) redoes the minimap2 search for gene copies instead of using the
  shipped masks.
- **Time:** about 30 minutes, peak memory 3.2 GB.

The final files were rebuilt by the script from scratch on 2026-09-26 (all six groups, 29 minutes, peak
3.2 GB; the non-human group's GenBank records from a local copy of the same accessions, the others fetched
again), and the rebuild reproduced all seven panel files, `shared_between_panels.tsv` and all 83 units byte
for byte. Separate runs of the coding-vntrs and sex-chromosome-arrays recipes gave the same k-mer rows.

## The final files

These are the files to load in the cohort's remaining scans. Each counts file records the sha256 of the
panels it was counted with (`panel_sha256`), and `ngsdose estimate` matches experimental units by it.

| panel file | k-mers | sha256 |
| --- | --- | --- |
| `coding-vntrs.k31.panel.tsv.gz` | 86,162 | 46bf3ddd0b48031d2d968bf433fe2a36d2cd7bdde8c71c4c1ec2acae39917a29 |
| `macrosatellites.k31.panel.tsv.gz` | 45,362 | 385bf0d44373a563770e3e4cb6cfc9bb272d202bc582c9a87d3b06be5c5eaf6f |
| `multicopy-genes.k31.panel.tsv.gz` | 544,098 | a7dafd323ae5c23dc825bd9f8375be9123d7883f906c2e20208fca4b5cfd2aa8 |
| `nonhuman-myco.k31.panel.tsv.gz` | 3,731,523 | ad14926aec70efba19cf87e60f795699b6a8487b58fba76df8d19db0a644450f |
| `nonhuman.k31.panel.tsv.gz` | 395,211 | 8c4a653669fc45979785203d27880e5d4ec89b58491e11d91912b9aff8070f4f |
| `rna-arrays.k31.panel.tsv.gz` | 148,026 | c33c4a131c8c5bee0b4d260619f94a265226f49af12193401c31d89ff1b9d1d3 |
| `sex-chromosome-arrays.k31.panel.tsv.gz` | 126,346 | d103151637f200c147dedf90aa6a2b5640b65e316f542bbc2a554745d1ff1aa8 |

Total 5,076,728 k-mers in 83 classes. Against the files of the previous round, `nonhuman`, `nonhuman-myco`
and `rna-arrays` are unchanged; `macrosatellites` (`REXO1L1`, `RNU2`), `multicopy-genes` (`DEFA1A3`),
`sex-chromosome-arrays` (`TSPY`, `OPN1`, `CT45`) and `coding-vntrs` (`VNTR_MUC19`, `VNTR_MUC6`) changed, as did
the units of the six re-cut or trimmed classes (see "Units that are one period").
