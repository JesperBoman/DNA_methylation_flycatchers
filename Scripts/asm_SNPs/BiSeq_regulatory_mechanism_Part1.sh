#!/bin/bash -l

ml bioinfo-tools R_packages/4.0.4

annot="fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_plus_minus_100bp.bed"
tissue="Heart Brain Liver Testis"

Rscript --vanilla BiSeq_regulatory_mechanism_Part1.R $annot $tissue
