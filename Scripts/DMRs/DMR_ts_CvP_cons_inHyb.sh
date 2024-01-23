#!/bin/bash -l


mkdir conserved_ts

ml bioinfo-tools BEDTools/2.29.2
dirC="collared/tissue_specific_results"
dirP="pied/tissue_specific_results"
dirH="hybrid/tissue_specific_results"

faidx="../methylseq_pipe/fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa.fai"
genome_len=$(awk '{sum+=$2} END{print sum}' $faidx)

resamples=1000

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for tissue in "${arr[@]}"
do

bedtools intersect -header -f 0.25 -r -wa -a $dirC/${tissue}_uniq0.25fr_all4comp -b $dirP/${tissue}_uniq0.25fr_all4comp > conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons
bedtools intersect -header -f 0.25 -r -wa -a $dirC/${tissue}_uniq0.25fr_all4comp -b $dirH/${tissue}_uniq0.25fr_all4comp > conserved_ts/${tissue}_uniq0.25fr_all4comp_CvH_cons
bedtools intersect -header -f 0.25 -r -wa -a $dirP/${tissue}_uniq0.25fr_all4comp -b $dirH/${tissue}_uniq0.25fr_all4comp > conserved_ts/${tissue}_uniq0.25fr_all4comp_PvH_cons
bedtools intersect -header -f 0.25 -r -wo -a $dirH/${tissue}_uniq0.25fr_all4comp -b conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons > conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons_inHyb


RANDOM=$(date +%N | cut -b4-9)



for i in $(seq 1 $resamples)
do

bedtools intersect -header -f 0.25 -r -wa -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirC/${tissue}_uniq0.25fr_all4comp ) -b $dirP/${tissue}_uniq0.25fr_all4comp | wc -l  >> tsCP_overlap_tmp &
bedtools intersect -header -f 0.25 -r -wa -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirC/${tissue}_uniq0.25fr_all4comp ) -b $dirH/${tissue}_uniq0.25fr_all4comp | wc -l  >> tsCH_overlap_tmp &
bedtools intersect -header -f 0.25 -r -wa -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirP/${tissue}_uniq0.25fr_all4comp ) -b $dirH/${tissue}_uniq0.25fr_all4comp | wc -l  >> tsPH_overlap_tmp &
bedtools intersect -header -f 0.25 -r -wa -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirH/${tissue}_uniq0.25fr_all4comp ) -b conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons | wc -l  >> tsCPH_overlap_tmp &
wait
done

tsCP_overlap=$(wc -l conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons | cut -f1 -d ' ' | awk '{print $1-1}')
tsCH_overlap=$(wc -l conserved_ts/${tissue}_uniq0.25fr_all4comp_CvH_cons | cut -f1 -d ' ' | awk '{print $1-1}')
tsPH_overlap=$(wc -l conserved_ts/${tissue}_uniq0.25fr_all4comp_PvH_cons | cut -f1 -d ' ' | awk '{print $1-1}')
tsCPH_overlap=$(wc -l conserved_ts/${tissue}_uniq0.25fr_all4comp_CvP_cons_inHyb | cut -f1 -d ' ' | awk '{print $1-1}')

tsCP_ln=$(cat tsCP_overlap_tmp <(echo $tsCP_overlap) | sort -n | grep -n -w "$tsCP_overlap" | cut -f1 -d ":" | head -n1 )
tsCH_ln=$(cat tsCH_overlap_tmp <(echo $tsCH_overlap) | sort -n | grep -n -w "$tsCH_overlap" | cut -f1 -d ":" | head -n1 )
tsPH_ln=$(cat tsPH_overlap_tmp <(echo $tsPH_overlap) | sort -n | grep -n -w "$tsPH_overlap" | cut -f1 -d ":" | head -n1 )
tsCPH_ln=$(cat tsCPH_overlap_tmp <(echo $tsCPH_overlap) | sort -n | grep -n -w "$tsCPH_overlap" | cut -f1 -d ":" | head -n1 )

tsCP_pval=$(awk -v tsCP_ln=$tsCP_ln -v resamples=$resamples 'BEGIN{r=resamples+1-tsCP_ln; print r/(resamples)}')
tsCH_pval=$(awk -v tsCH_ln=$tsCH_ln -v resamples=$resamples 'BEGIN{r=resamples+1-tsCH_ln; print r/(resamples)}')
tsPH_pval=$(awk -v tsPH_ln=$tsPH_ln -v resamples=$resamples 'BEGIN{r=resamples+1-tsPH_ln; print r/(resamples)}')
tsCPH_pval=$(awk -v tsCPH_ln=$tsCPH_ln -v resamples=$resamples 'BEGIN{r=resamples+1-tsCPH_ln; print r/(resamples)}')

tsCP_resamp_mean=$(awk '{sum+=$1} END{print sum/NR}' tsCP_overlap_tmp)
tsCH_resamp_mean=$(awk '{sum+=$1} END{print sum/NR}' tsCH_overlap_tmp)
tsPH_resamp_mean=$(awk '{sum+=$1} END{print sum/NR}' tsPH_overlap_tmp)
tsCPH_resamp_mean=$(awk '{sum+=$1} END{print sum/NR}' tsCPH_overlap_tmp)

awk  -v tissue=$tissue -v tsCP_overlap=$tsCP_overlap -v tsCH_overlap=$tsCH_overlap -v tsPH_overlap=$tsPH_overlap \
	 -v tsCPH_overlap=$tsCPH_overlap -v tsCP_pval=$tsCP_pval -v tsCH_pval=$tsCH_pval -v tsPH_pval=$tsPH_pval -v tsCPH_pval=$tsCPH_pval  \
	 -v tsCP_resamp_mean=$tsCP_resamp_mean -v tsCH_resamp_mean=$tsCH_resamp_mean -v tsPH_resamp_mean=$tsPH_resamp_mean -v tsCPH_resamp_mean=$tsCPH_resamp_mean \
	'BEGIN{OFS="\t"; print tissue, tsCP_overlap, tsCH_overlap, tsPH_overlap,  tsCPH_overlap, tsCP_resamp_mean, tsCH_resamp_mean, tsPH_resamp_mean, tsCPH_resamp_mean, tsCP_pval, tsCH_pval, tsPH_pval,  tsCPH_pval  }' >> conserved_ts/stats


rm tsCP_overlap_tmp
rm tsCH_overlap_tmp
rm tsPH_overlap_tmp
rm tsCPH_overlap_tmp

done
