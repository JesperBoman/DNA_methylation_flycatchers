#!/bin/bash -l


#This script creates an ancestral genome sequence (.fasta) based on a list of inferred ancestral alleles 

GENOME="fAlb15.chrom.fa"


ml bioinfo-tools BEDTools/2.29.2

SNP_list="shared.ancestral.states.chrom"
grep -w --no-group-separator "A" $SNP_list > snps_A

n=1
awk '{print $1 "\t" $2-1 "\t" $2}' snps_A > snpsA.bed


bedtools maskfasta -fi $GENOME \
 -bed snpsA.bed \
 -fo fAlb_ancestral${n}.fa \
 -mc A

rm snps_A
rm snpsA.bed

o=$((n+1))
grep -w --no-group-separator "T" $SNP_list > snps_T
awk '{print $1 "\t" $2-1 "\t" $2}' snps_T > snpsT.bed

bedtools maskfasta -fi fAlb_ancestral${n}.fa \
 -bed snpsT.bed \
 -fo fAlb_ancestral${o}.fa \
 -mc T

rm fAlb_ancestral${n}.fa
rm snps_T
rm snpsT.bed

n=$o
o=$((n+1))
grep -w --no-group-separator "C" $SNP_list > snps_C
awk '{print $1 "\t" $2-1 "\t" $2}' snps_C > snpsC.bed

bedtools maskfasta -fi fAlb_ancestral${n}.fa \
 -bed snpsC.bed \
 -fo fAlb_ancestral${o}.fa \
 -mc C

rm fAlb_ancestral${n}.fa
rm snps_C
rm snpsC.bed

n=$o
grep -w --no-group-separator "G" $SNP_list > snps_G
awk '{print $1 "\t" $2-1 "\t" $2}' snps_G > snpsG.bed

bedtools maskfasta -fi fAlb_ancestral${n}.fa \
 -bed snpsG.bed \
 -fo fAlb_ancestral.fa \
 -mc G

rm fAlb_ancestral${n}.fa
rm snps_G
rm snpsG.bed
