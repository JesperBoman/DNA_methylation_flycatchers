#!/bin/bash -l

#LiftOver of fixed differences from fAlb15scaf to fAlb15chrom

fixd="temp/fixed_differences_Oland_alexn_vcf_Jesper.txt"

#This is a script made by Linnéa Smeds, it does a liftover between scaffold and chromosome assemblies
perl translateScaffToChrom.pl -in="$fixd" -version=20130221 -level=strict -scafCol=1 -posCol=2 -out=fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper

awk '{OFS="\t"; print $1, $2-1, $2}' fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper > fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper.bed

#wc -l fixed_diff_chr_all_OC_OP_alexn_vcf_ludoscript_Jesper
#38923

#Bed of of poly-C and poly-G tracts
#According to Linnéa Smeds, the assembler of the fAlb15 genome put long poly-G/C tracts at scaffold breaks
awk 'NR>3{if(($10  == "(G)n" || $10  == "(C)n") && ($7-$6 >100)){OFS="\t"; print $5, $6-1, $7}}' fAlb15.chrom.fa.out > poly_C_and_G_100bp_plus.bed

#Masking the reference

ml bioinfo-tools BEDTools/2.29.2

bedtools maskfasta -fi fAlb15.chrom.fa -bed <(cat fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper.bed poly_C_and_G_100bp_plus.bed) -fo fAlb15.chrom_Olandfix_polyCorG100bpthres_hardmasked.fa

bedtools maskfasta -fi fAlb15.chrom.fa -bed fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper.bed -fo tmp.fa

bedtools maskfasta -fi tmp.fa -bed poly_C_and_G_100bp_plus.bed -soft -fo fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa
