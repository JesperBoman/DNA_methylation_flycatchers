#!/bin/bash -l

ml bioinfo-tools R_packages/4.0.4


#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME"
#annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot"
annot_seg="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot_ExIn_HDMR"


declare -a tiss=("Heart" "Kidney" "Liver" "Brain" "Testis")


randnum=$(echo $RANDOM)

for tissue in "${tiss[@]}"
do

grep "$tissue" $annot_seg > temp_${tissue}_${randnum}


p=$(($p+1))


#Rscript --vanilla Gene_profile_CvP_DMR.R temp_${tissue}_${randnum} $tissue &
#Rscript --vanilla Gene_profile_CvP_DMR_l2fcCOLvPIE.R temp_${tissue}_${randnum} $tissue &
#Rscript --vanilla Gene_profile_CvP_DMR_DEstat.R temp_${tissue}_${randnum} $tissue &
#Rscript --vanilla Gene_profile_Fst.R temp_${tissue}_${randnum} $tissue &
#Rscript --vanilla Gene_profile_HDMR.R temp_${tissue}_${randnum} $tissue &
#Rscript --vanilla Gene_profile_HybFst.R temp_${tissue}_${randnum} $tissue &
Rscript --vanilla Gene_profile_DE_status_CGI_all.R temp_${tissue}_${randnum} $tissue &

pids[${p}]=$!

if [[ ${#pids[@]} -gt 4 ]] ; then

wait ${pids[$p-4]}
unset pids[$p-4]

fi


done

wait


for tissue in "${tiss[@]}"
do

mv *.png gene_profile_plots
mv *.rda rdata_storage

rm temp_${tissue}_${randnum}

done
