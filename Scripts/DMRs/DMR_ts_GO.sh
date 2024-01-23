#!/bin/bash -l



mkdir ts_GO

module load bioinfo-tools BEDTools/2.29.2

dirC="collared/tissue_specific_results"
dirP="pied/tissue_specific_results"
dirH="hybrid/tissue_specific_results"

geneflanksbed="../annotation/maker/run4/fAlb15.maker.output/mRNA_parent_name_mirror.sp.prom.gene5kbflanks.bed"

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for tissue in "${arr[@]}"
do
bedtools intersect -wa -b $dirC/${tissue}_uniq0.25fr_all4comp -a $geneflanksbed | awk '{if($6 != "NA"){print $6}}' > ts_GO/COL_${tissue}_sp.list
bedtools intersect -wa -b $dirP/${tissue}_uniq0.25fr_all4comp -a $geneflanksbed | awk '{if($6 != "NA"){print $6}}' > ts_GO/PIE_${tissue}_sp.list
bedtools intersect -wa -b $dirH/${tissue}_uniq0.25fr_all4comp -a $geneflanksbed | awk '{if($6 != "NA"){print $6}}' > ts_GO/HYB_${tissue}_sp.list

bedtools intersect -wa -b conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons_inHyb -a $geneflanksbed | awk '{if($6 != "NA"){print $6}}' > ts_GO/CvP_cons_inHyb_${tissue}_sp.list
bedtools intersect -wa -b conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons_inHyb -a $geneflanksbed | awk '{if($5 ~ /0:/){split($5, a, ":"); print a[3]} else if($5 ~ /maker/){print $6 } else{print $5}}' > ts_GO/CvP_cons_inHyb_${tissue}_gene_name.list
done
