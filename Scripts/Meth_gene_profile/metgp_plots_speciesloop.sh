#!/bin/bash -l


ml bioinfo-tools R_packages/4.0.4


#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot"
annot_seg="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot_ExIn"

declare -a spec=("COL" "HYB" "PIE")


randnum=$(echo $RANDOM)

for species in "${spec[@]}"
do


if [[ $species == "COL" ]] ; then
echo $species
grep -E 'COL|HYB04' $annot_seg  > temp_${species}_${randnum}

fi

if [[ $species == "PIE" ]] ; then
echo $species
grep "$species" $annot_seg  > temp_${species}_${randnum}

fi

if [[ $species == "HYB" ]] ; then
echo $species
grep "$species" $annot_seg | grep -v "HYB04" > temp_${species}_${randnum}

fi


#Exclude brain alternative
grep -v "Brain" temp_${species}_${randnum} > tmpi
mv tmpi temp_${species}_${randnum}
#

Rscript --vanilla Gene_profile_acrossTissCor.R temp_${species}_${randnum} $species


done



for species in "${spec[@]}"
do

mv *.png gene_profile_plots
mv *.rda rdata_storage

rm temp_${species}_${randnum}

done
