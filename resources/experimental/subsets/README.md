# Sub-options: named arrays within a satellite family's sinks

Each `NAME.bed` here is a subset of the learned sink intervals of one class in
`../sinks.satellites.bed`, with that class in the fourth column. The fetch menu
(`../../fetch_menu.tsv`, kind `subset`, preset `xy_arrays`) offers each one as an
option. A fetch reads its intervals and counts them as the class. `ngsdose estimate` then reports
`NAME.reads`, `NAME.mass_Mb` and `NAME.status` for every counts file, scan or fetch. The value
comes from the class's reads placed inside the intervals. A placement bin counts only when all
of it, clipped at the contig's end, lies inside. The mass is those reads times the class's mass
per read.

The intervals come from where the NYGC bwa-mem alignments put each family's reads in
whole-file scans, not from reference coordinates. Like the sinks they are cut from, they belong
to that aligner and reference.

| option | class | intervals | MB of a 30x CRAM | measures | truth |
|---|---|---|---|---|---|
| DXZ1 | aSatHOR | 4 (4.21 Mb), chrX:58.0-63.0 Mb | 9.8 (7.2-15.8) | DXZ1 array mass over all X chromosomes | HPRC r2 chrX HOR array, both haplotypes: n = 124, r 0.992; women / men 2.02 |
| DYZ3 | aSatHOR | 4 (0.38 Mb), chrY:10.0-10.7 Mb | 5.0 (2.6-7.2) | DYZ3 array mass (chrY centromere) | sex; HPRC r2 chrY HOR array, men: n = 60, r 0.835 |
| DYZ1 | HSat3 | 1 (0.11 Mb), chrY:56.67-56.78 Mb | 39.1 (1.0-75.5) | Yq12 HSat3 mass | sex only |
| DYZ2 | HSat1B | 78 (0.91 Mb), mostly autosomes, chrX, decoys (2 on chrY) | 34.3 (22.6-52.1) | an index of male-specific (Y-derived) HSat1B: the part in 78 male-only intervals, 21.7% of men's HSat1B, about a quarter of the Y-derived part, not the whole of it | sex only |

MB is the median (range) over 13 NYGC CRAM indexes of the bytes of the slices the intervals
overlap, with their containers' compression headers, as `ngsdose fetchplan` prices them: what
`ngs-dose count -m fetch` reads, one fetch per run of touching or overlapping intervals, a slice
under several runs decoded once per run (DXZ1's four intervals and DYZ3's four share slices, so
they cost more than the slices read once: 9.5 and 3.3 MB). DYZ2's
intervals were chosen on the 100 training scans of the satellite sinks: men's median at least
100 HSat1B reads, and women's median under 2% of men's. Y-derived HSat1B mostly does not land on chrY:
the 13 chrY HSat1B sinks hold 0.16% of it. Of the 78 intervals, 66 are autosomal, 7 on chrX, 3 on
decoys and 2 on chrY. The HPRC release-2 assemblies rarely close Yq12, so
DYZ1 and DYZ2 have no assembly truth: against them r is 0.23 (DYZ1, n = 60) and 0.07 (DYZ2, n =
27). The measured shares and separations by sex, over the 1,648 cohort scans held out from the
satellite sinks, are in each file's header.

A fetch that reads a sub-option without the rest of its class's sinks counts the class only
there. `ngsdose estimate` then marks the class `subset_only` and does not measure it. It can tell
this from the fetch's sinks BED when it knows that file: the bundle's and these are known, and a
fetchplan `PREFIX.sinks.bed` is passed with `--fetch-sinks`. Without the BED, the sub-options of
a fetch are `unverified` (reads given, mass NaN), and so are the three families with sub-options,
`aSatHOR`, `HSat1B` and `HSat3` (NaN), since the fetch may have read them only at their sub-options;
`--fetch-sinks PREFIX.sinks.bed` settles both.

The MB column is each sub-option's cost alone over the 13 NYGC indexes, written when the files were
made; [`docs/fetch_examples.md`](../../../docs/fetch_examples.md) (examples 7 and 15) prices plans with
them and is the figure to re-run. The headers of the BED files themselves still carry the costs of the
model before 2026-09-28 (every slice once); `ngsdose estimate` knows these files by their hash, so
they are not rewritten for a comment.
