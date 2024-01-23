#!/bin/bash -l


mkdir conserved_per_tissue

ml bioinfo-tools BEDTools/2.29.2
dirC="collared"
dirP="pied"
dirH="hybrid"


declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")
declare -a arr2=("Heart" "Kidney" "Liver" "Brain" "Testis")

p=0
for foc in "${arr[@]}"
do

#Remove the focal tissue from this array
unset 'arr2[$p]'
p=$(($p+1))


for i in "${arr2[@]}"
do

Tissue1=$foc
Tissue2=$i


bedtools intersect -header -f 0.25 -r -wa -a $dirC/${Tissue1}v${Tissue2}/dmrs.txt -b $dirP/${Tissue1}v${Tissue2}/dmrs.txt > conserved_per_tissue/${Tissue1}v${Tissue2}_0.25fr_CvP_cons

bedtools intersect -header -f 0.25 -r -wo -a $dirH/${Tissue1}v${Tissue2}/dmrs.txt -b conserved_per_tissue/${Tissue1}v${Tissue2}_0.25fr_CvP_cons > conserved_per_tissue/${Tissue1}v${Tissue2}_0.25fr_CvP_cons_inHyb

#Count, mean number of CpGs, mean width, meanDiff
awk -v tiss1="$Tissue1" -v tiss2="$Tissue2" 'NR>1{CpG+=$7; MW+=$8; MD+=$12} END{print tiss1 "v" tiss2 "\t" NR-1 "\t" CpG/(NR-1) "\t" MW/(NR-1) "\t" MD/(NR-1)}' conserved_per_tissue/${Tissue1}v${Tissue2}_0.25fr_CvP_cons >> conserved_per_tissue/stats_cons_p_t_CvP

awk -v tiss1="$Tissue1" -v tiss2="$Tissue2" 'NR>1{CpG+=$7; MW+=$8; MD+=$12} END{print tiss1 "v" tiss2 "\t" NR-1 "\t" CpG/(NR-1) "\t" MW/(NR-1) "\t" MD/(NR-1)}' conserved_per_tissue/${Tissue1}v${Tissue2}_0.25fr_CvP_cons_inHyb >> conserved_per_tissue/stats_cons_p_t_CvP_inHyb

done
done
