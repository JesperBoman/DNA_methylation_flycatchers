#!/bin/bash -l


module load bioinfo-tools bcftools/1.12

Oland_collared="gt_OC.vcf.gz"
Oland_pied="gt_OP.vcf.gz"


bcftools merge -O z $Oland_collared $Oland_pied > OC.OP.comb.vcf.gz

bcftools query -l $Oland_collared > COL_samples
bcftools query -l $Oland_pied > PIE_samples

module load vcftools/0.1.16 

bedfile="genes_5n3p_5kb_flanks_100bp_seg_with_header.bed"
pop1="COL_samples"
pop2="PIE_samples"

vcftools --gzvcf --weir-fst-pop --bed --fst-window-size 100
