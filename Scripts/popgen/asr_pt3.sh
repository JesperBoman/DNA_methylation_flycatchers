#!/bin/bash -l


module load bioinfo-tools vcftools/0.1.16 

vcftools --gzvcf OC.OP.comb.vcf.gz --counts --stdout | awk 'NR>1{print $1 "\t" $2}' > OC.OP.SNP_list

perl translateScaffToChrom.pl -in=OC.OP.SNP_list -version=20130221 -level=strict -scafCol=1 -posCol=2 -out=OC.OP.SNP_list.chrom
