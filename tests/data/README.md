# Test data

`NA12878.subsample.bam` (+ `.csi`): a 2% template subsample (`samtools view -s 1.02`) of NA12878
from the 1000 Genomes 30× resource (NYGC; ENA ERR3239334; Byrska-Bishop et al., *Cell* 2022),
restricted to the GRCh38-v1 bundle's control regions and its rDNA45S, rDNA5S and DJ sinks as of
2026-09-19. Read names, base qualities and tags were removed and contigs without reads dropped from
the header; flags, positions, CIGARs, mate fields and sequences are as released. 198,848 reads;
~0.75× over the 800 control regions (duplicates included; ~0.67× without) and ~1.5× over all
retained regions, which include the rDNA and DJ sinks. It holds no read without a coordinate.
The `TEL` sinks added on 2026-09-22 are not covered, except the one that coincides with an rDNA45S
sink (chr2:32,909,000-32,921,000); the other 62 hold no fixture reads. `make_fixture.sh` rebuilds
it from the full CRAM (needed only when the bundle's control regions or sinks change, and to cover
the `TEL` sinks). 1000 Genomes data are open access (https://www.internationalgenome.org/IGSR_disclaimer).

It lets CI run the real engine and estimator on real reads and assert the method's own claims:
known-copy-number sequence comes out at its known copy number (`tests/test_real_data.py`).

`ngspca_1000G.singularvalues.txt`: the 200 singular values NGS-PCA kept of its 1000 Genomes run (the PCs and bins themselves are in NGS-DOSE-1000G/meta/ngspca); the Marchenko-Pastur edge test reads them.
