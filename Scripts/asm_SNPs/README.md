This is the pipeline for determining the molecular mechanism of DNA methylation evolution.

1. First install [SNPsplit](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/).

2. Make sure that you have masked your marker loci in your reference genome (see instructions for SNPsplit for more details)

3. Then run asm_snp.sh

From here on we use the R package [BiSeq](https://www.bioconductor.org/packages/release/bioc/html/BiSeq.html) to test for statistical differences between sample groups and subgenomes of the hybrids

4. BiSeq_regulatory_mechanism_part1.R, can be run on a cluster (see BiSeq_regulatory_mechanism_part1.sh). Here is where the statistical differences are being assessed.

5. BiSeq_regulatory_mechanism_part2.R, takes some time but I run directly on my desktop. This script is used to classify the molecular mechanism and follow-up plots and stats.

--

ams_snp_CXparser.sh and ams_snp_parentals_meth.sh are not necessarily needed but can be used to calculate methylation levels.
