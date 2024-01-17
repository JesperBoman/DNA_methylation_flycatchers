#!/bin/bash -l

ml bioinfo-tools R_packages/4.0.4






#Between-species comparisons ####



declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for i in "${arr[@]}"
do

Tissue=$i


args=$(ls bsseq_formatted_data | grep -E "$Tissue" | grep -E 'COL|HYB04' | awk '{printf $1 "\t"}')

#Loading data and creating BSseq object 
Rscript DMR_analysis.R $args

mv BS_data.rda BS_data_${Tissue}.rda
mv BS_data.fit.rda BS_data_${Tissue}.fit.rda 

mv BS_data_${Tissue}.rda bsseq_read_rda
mv BS_data_${Tissue}.fit.rda  bsseq_fit_rda

#Obtaining t-stats and finding DMRs
cd $Tissue
Rscript ../DMR_analysis_findDMRs-${Tissue}.R



done




#Between-tissue comparisons (done per sample group) ####

dir="."
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

cd $dir

Tissue1=$foc
Tissue2=$i




#Collared
#argsT1=$(ls bsseq_formatted_data | grep -E "$Tissue1" | grep -E 'COL|HYB04' | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
#argsT2=$(ls bsseq_formatted_data | grep -E "$Tissue2" | grep -E 'COL|HYB04' | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')

#Pied
argsT1=$(ls bsseq_formatted_data | grep -E "$Tissue1" | grep -E 'PIE' | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
argsT2=$(ls bsseq_formatted_data | grep -E "$Tissue2" | grep -E 'PIE' | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')

#Hybrid
#argsT1=$(ls bsseq_formatted_data | grep -E "$Tissue1" | grep -E 'HYB' | grep -v 'HYB04' | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
#argsT2=$(ls bsseq_formatted_data | grep -E "$Tissue2" | grep -E 'HYB' | grep -v 'HYB04' | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')

Rscript DMR_analysis.R $argsT1 $argsT2

#Collared
#mv BS_data.rda BS_data_C_${Tissue1}v${Tissue2}.rda
#mv BS_data.fit.rda BS_data_C_${Tissue1}v${Tissue2}.fit.rda

#mv BS_data_C_${Tissue1}v${Tissue2}.rda bsseq_read_rda
#mv BS_data_C_${Tissue1}v${Tissue2}.fit.rda bsseq_fit_rda

#cd $dir/collared/${Tissue1}v${Tissue2}

#Rscript ../../DMR_analysis_findDMRs-Collared.R $dir/bsseq_fit_rda/BS_data_C_${Tissue1}v${Tissue2}.fit.rda $dir/bsseq_read_rda/BS_data_C_${Tissue1}v${Tissue2}.rda

#Pied
mv BS_data.rda BS_data_P_${Tissue1}v${Tissue2}.rda
mv BS_data.fit.rda BS_data_P_${Tissue1}v${Tissue2}.fit.rda

mv BS_data_P_${Tissue1}v${Tissue2}.rda bsseq_read_rda
mv BS_data_P_${Tissue1}v${Tissue2}.fit.rda bsseq_fit_rda

cd $dir/pied/${Tissue1}v${Tissue2}

Rscript ../../DMR_analysis_findDMRs-Pied.R $dir/bsseq_fit_rda/BS_data_P_${Tissue1}v${Tissue2}.fit.rda $dir/bsseq_read_rda/BS_data_P_${Tissue1}v${Tissue2}.rda

#Hybrid
#mv BS_data.rda BS_data_H_${Tissue1}v${Tissue2}.rda
#mv BS_data.fit.rda BS_data_H_${Tissue1}v${Tissue2}.fit.rda

#mv BS_data_H_${Tissue1}v${Tissue2}.rda bsseq_read_rda
#mv BS_data_H_${Tissue1}v${Tissue2}.fit.rda bsseq_fit_rda

#cd $dir/hybrid/${Tissue1}v${Tissue2}

#Rscript ../../DMR_analysis_findDMRs-Hybrid.R $dir/bsseq_fit_rda/BS_data_H_${Tissue1}v${Tissue2}.fit.rda $dir/bsseq_read_rda/BS_data_H_${Tissue1}v${Tissue2}.rda $Tissue1 $Tissue2

done
done




#NOTES
#10 cores not enough when smoothing 14 samples. Use 20 cores instead.
