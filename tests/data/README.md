# Test data

`NA12878.subsample.bam` (+ `.csi`): a 2% template subsample (`samtools view -s 1.02`) of NA12878
from the 1000 Genomes 30× resource (NYGC; ENA ERR3239334; Byrska-Bishop et al., *Cell* 2022),
restricted to the GRCh38-v1 bundle's control regions and class sinks. Read names, base qualities
and tags were removed and contigs without reads dropped from the header; flags, positions,
CIGARs, mate fields and sequences are as released. 198,848 reads, ~0.75× over the retained
regions. `make_fixture.sh` rebuilds it from the full CRAM (needed only when the bundle's control
regions or sinks change). 1000 Genomes data are open access (https://www.internationalgenome.org/IGSR_disclaimer).

It lets CI run the real engine and estimator on real reads and assert the method's own claims:
known-copy-number sequence comes out at its known copy number (`tests/test_real_data.py`).
