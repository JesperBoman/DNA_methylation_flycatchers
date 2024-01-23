#!/bin/bash -l

GENOME="fAlb_ancestral.fa"
SNP_list="OC.OP.SNP_list.chrom"

ml bioinfo-tools samtools/1.12

#This code parallellizes the process of finding ancestral CpG sites
#It basically splits the SNP list into ten chunks and runs triNuc_filter.awk separately on those chunks

nrows=$(wc -l $SNP_list | cut -f1 -d ' ')

for ((i=0; i<=9; i++))
do
lower=$((i* nrows/ 10))
upper=$(((i+1) * nrows /10))
awk -v genome=$GENOME -f triNuc_filter.awk <(awk -v lower=$lower -v upper=$upper '{if (NR > lower && NR <= upper)print $0}' $SNP_list ) > temp${i} &
done
wait

for ((i=0; i<=9; i++))
do
awk '{print $0}' temp${i} >> ${SNP_list}_CpG_annot
rm temp${i}
done
