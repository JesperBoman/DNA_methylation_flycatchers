#!/bin/bash -l


ml bioinfo-tools R_packages/4.0.4


#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr"
annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot_ExIn"

declare -a tiss=("Heart" "Kidney" "Liver" "Brain" "Testis")
declare -a spec=("COL" "PIE" "HYB")

randnum=$(echo $RANDOM)

for species in "${spec[@]}"
do
for tissue in "${tiss[@]}"
do

if [[ $species == "COL" ]] ; then
echo $species
#awk -v tissue="$tissue"  '{if(($13 ~/COL/ || $13 == "HYB04") && $14 == tissue)print $0} '$annot_seg > temp_${species}_${tissue}
grep -E 'COL|HYB04' $annot_seg | grep "$tissue" > temp_${species}_${tissue}_${randnum}

fi

if [[ $species == "PIE" ]] ; then
echo $species
#awk -v tissue="$tissue"  '{if($13 ~/PIE/ && $14 == tissue)print $0}' $annot_seg > temp_${species}_${tissue}
grep "$species" $annot_seg | grep "$tissue" > temp_${species}_${tissue}_${randnum}

fi

if [[ $species == "HYB" ]] ; then
echo $species
#awk -v tissue="$tissue"  '{if($13 ~/HYB/ && $13 != HYB04 && $14 == tissue)print $0}' $annot_seg > temp_${species}_${tissue}
grep "$species" $annot_seg | grep "$tissue" | grep -v "HYB04" > temp_${species}_${tissue}_${randnum}

fi

p=$(($p+1))

#Rscript --vanilla Gene_profile_plots_comb_2021.R temp_${species}_${tissue} $species $tissue
#Rscript --vanilla Gene_profile_comb_2021_promoter_type.R temp_${species}_${tissue}_${randnum} $species $tissue &
#Rscript --vanilla Gene_profile_Meth_x_RNA_cor.R temp_${species}_${tissue}_${randnum} $species $tissue &
#Rscript --vanilla Gene_profile_Meth_x_RNA_LMH.R temp_${species}_${tissue}_${randnum} $species $tissue &
#Rscript --vanilla Gene_profile_Meth_x_RNA_Inpat.R temp_${species}_${tissue}_${randnum} $species $tissue &
Rscript --vanilla Gene_profile_Exin_DE.R temp_${species}_${tissue}_${randnum} $species $tissue &


pids[${p}]=$!

if [[ ${#pids[@]} -gt 4 ]] ; then

wait ${pids[$p-4]}
unset pids[$p-4]

fi


done
done

wait

declare -a spec=("COL" "PIE" "HYB")

for species in "${spec[@]}"
do
for tissue in "${tiss[@]}"
do

mv *.png gene_profile_plots
mv *.rda rdata_storage
rm temp_${species}_${tissue}_${randnum}

done
done
