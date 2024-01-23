#!/bin/bash -l
#SBATCH -J tsDMR_tsEXPR_association
#SBATCH -o tsDMR_tsEXPR_association.output
#SBATCH -e tsDMR_tsEXPR_association.error
#SBATCH --mail-user jesper.boman@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 00-12:00:00
#SBATCH -A p2018002
#SBATCH -p core
#SBATCH -n 3




ml bioinfo-tools BEDTools/2.29.2  R_packages/4.0.4
faidx="../methylseq_pipe/fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa.fai"
dir="."
PEM_data="../Meth_x_RNA/PEM.data_species_mean"

dirC="collared/tissue_specific_results"
dirP="pied/tissue_specific_results"
dirH="hybrid/tissue_specific_results"
declare -a tissues=("Heart" "Kidney" "Liver" "Brain" "Testis")
declare -a species_list=("COL" "PIE" "HYB")

while IFS= read -r annot_double
do

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annot=$(cut -f2 <(echo $annot_double) -d ",")

awk '{print $1 "\t" $2 "\t" $3}' $annot > tmp_annot
annot=tmp_annot


bedtools intersect -v -wa -a $annot -b $dirC/Heart_uniq0.25fr_all4comp_sameDir  -b $dirC/Kidney_uniq0.25fr_all4comp_sameDir -b $dirC/Liver_uniq0.25fr_all4comp_sameDir  -b $dirC/Brain_uniq0.25fr_all4comp_sameDir -b $dirC/Testis_uniq0.25fr_all4comp_sameDir | awk '{print $1 "_" $2 "_" $3 }' > tmpCOL_nonTSproms
bedtools intersect -v -wa -a $annot -b $dirP/Heart_uniq0.25fr_all4comp_sameDir  -b $dirP/Kidney_uniq0.25fr_all4comp_sameDir -b $dirP/Liver_uniq0.25fr_all4comp_sameDir  -b $dirP/Brain_uniq0.25fr_all4comp_sameDir -b $dirP/Testis_uniq0.25fr_all4comp_sameDir | awk '{print $1 "_" $2 "_" $3 }' > tmpPIE_nonTSproms
bedtools intersect -v -wa -a $annot -b $dirH/Heart_uniq0.25fr_all4comp_sameDir  -b $dirH/Kidney_uniq0.25fr_all4comp_sameDir -b $dirH/Liver_uniq0.25fr_all4comp_sameDir  -b $dirH/Brain_uniq0.25fr_all4comp_sameDir -b $dirH/Testis_uniq0.25fr_all4comp_sameDir | awk '{print $1 "_" $2 "_" $3 }' > tmpHYB_nonTSproms


echo $shortAnnot

for tissue in "${tissues[@]}"
do

bedtools intersect -wa -a $annot -b $dirC/${tissue}_uniq0.25fr_all4comp_sameDir | awk '{print $1 "_" $2 "_" $3 }' > tmpCOL_proms
bedtools intersect -wa -a $annot -b $dirP/${tissue}_uniq0.25fr_all4comp_sameDir | awk '{print $1 "_" $2 "_" $3 }' > tmpPIE_proms
bedtools intersect -wa -a $annot -b $dirH/${tissue}_uniq0.25fr_all4comp_sameDir | awk '{print $1 "_" $2 "_" $3 }' > tmpHYB_proms

for species in "${species_list[@]}"
do

awk 'NR==FNR{a[$1]} NR!=FNR && FNR>1{if($4 in a){print $0 "\t" "nonTS"}}' tmp${species}_nonTSproms <(grep "$species" $PEM_data) >  tmp_nonTSproms_PEM
awk 'NR==FNR{a[$1]} NR!=FNR && FNR>1{if($4 in a){print $0 "\t" "TS"}}' tmp${species}_proms <(grep "$tissue" $PEM_data | grep "$species" ) > tmp_TSproms_PEM
awk 'NR==FNR{a[$1]} NR!=FNR && FNR>1{if($4 in a){print $0 "\t" "REF"}}' tmp${species}_proms <(grep -v "$tissue" $PEM_data | grep "$species" ) > ref_test

echo $tissue $species
Rscript tsDMR_tsEXPR_association_stats.R $tissue >> tsDMR_tsEXPR_association_stats.results_TPM0plus

rm tmp${species}_proms

done
done

rm tmp${species}_nonTSproms

rm tmp_annot

done < "annotation.list.csv"
