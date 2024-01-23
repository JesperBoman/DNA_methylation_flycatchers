#!/bin/bash -l


## DMR overlap analysis ##

mkdir annotation_overlap

ml bioinfo-tools BEDTools/2.29.2
faidx="../methylseq_pipe/fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa.fai"
dir="."


dirC="collared/tissue_specific_results"
dirP="pied/tissue_specific_results"
dirH="hybrid/tissue_specific_results"


while IFS= read -r annot_double
do

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annot=$(cut -f2 <(echo $annot_double) -d ",")

#Bedtools merge here is a convenient heuristic since annotation features can overlap but are bound to a specific gene such as is the case for UTRs
#Genes can also overlap and it is probably fine to consider them as one long gene region instead of counting them as a double overlap.
awk '{print $1 "\t" $2 "\t" $3}' $annot | sort -k1,1 -k2,2n | bedtools merge > tmp_annot
annot=tmp_annot

annot_len=$(awk '{sum+=($3-$2)} END{print sum}' $annot)
genome_len=$(awk '{sum+=$2} END{print sum}' $faidx)


declare -a tissues=("Heart" "Kidney" "Liver" "Brain" "Testis")

RANDOM=$(date +%N | cut -b4-9)

for tissue in "${tissues[@]}"
do

C_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dirC/${tissue}_uniq0.25fr_all4comp)
P_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dirP/${tissue}_uniq0.25fr_all4comp)
H_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dirH/${tissue}_uniq0.25fr_all4comp)

C_overlap=$(bedtools intersect -header -wo -a $dirC/${tissue}_uniq0.25fr_all4comp -b $annot | cut -f22 | awk '{sum+=$1}END{print sum}')
P_overlap=$(bedtools intersect -header -wo -a $dirP/${tissue}_uniq0.25fr_all4comp -b $annot | cut -f22 | awk '{sum+=$1}END{print sum}')
H_overlap=$(bedtools intersect -header -wo -a $dirH/${tissue}_uniq0.25fr_all4comp -b $annot | cut -f22 | awk '{sum+=$1}END{print sum}')

resamples=1000

for i in $(seq 1 $resamples)
do
bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirC/${tissue}_uniq0.25fr_all4comp ) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> C_resample_overlap &

bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirP/${tissue}_uniq0.25fr_all4comp) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> P_resample_overlap &

bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dirH/${tissue}_uniq0.25fr_all4comp) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> H_resample_overlap &
	wait
done

C_ln=$(cat C_resample_overlap <(echo $C_overlap) | sort -n | grep -n -w "$C_overlap" | cut -f1 -d ":" | head -n1 )
P_ln=$(cat P_resample_overlap <(echo $P_overlap) | sort -n | grep -n -w "$P_overlap" | cut -f1 -d ":" | head -n1 )
H_ln=$(cat H_resample_overlap <(echo $H_overlap) | sort -n | grep -n -w "$H_overlap" | cut -f1 -d ":" | head -n1 )

C_pval=$(awk -v C_ln=$C_ln -v resamples=$resamples 'BEGIN{r=resamples+1-C_ln; print r/(resamples)}')
P_pval=$(awk -v P_ln=$P_ln -v resamples=$resamples 'BEGIN{r=resamples+1-P_ln; print r/(resamples)}')
H_pval=$(awk -v H_ln=$H_ln -v resamples=$resamples 'BEGIN{r=resamples+1-H_ln; print r/(resamples)}')

#Output: tissue, length of all DMRs, length of overlap between DMR and annotation, fraction of genome length of annotation, odds ratio (enrichment)
awk -v tissue=$tissue -v annot_len=$annot_len -v genome_len=$genome_len -v C_DMR_len=$C_DMR_len -v C_overlap=$C_overlap -v C_pval=$C_pval -v shortAnnot=$shortAnnot 'BEGIN{print tissue "\t" C_DMR_len "\t" C_overlap  "\t" annot_len/genome_len "\t" (C_overlap/C_DMR_len)/(annot_len/genome_len) "\t" C_pval "\t" "COL" "\t" shortAnnot}' >> annotation_overlap/stats_COL_ts_${shortAnnot}
awk -v tissue=$tissue -v annot_len=$annot_len -v genome_len=$genome_len -v P_DMR_len=$P_DMR_len -v P_overlap=$P_overlap -v P_pval=$P_pval -v shortAnnot=$shortAnnot 'BEGIN{print tissue "\t" P_DMR_len "\t" P_overlap  "\t" annot_len/genome_len "\t" (P_overlap/P_DMR_len)/(annot_len/genome_len) "\t" P_pval "\t" "PIE" "\t" shortAnnot}' >> annotation_overlap/stats_PIE_ts_${shortAnnot}
awk -v tissue=$tissue -v annot_len=$annot_len -v genome_len=$genome_len -v H_DMR_len=$H_DMR_len -v H_overlap=$H_overlap -v H_pval=$H_pval -v shortAnnot=$shortAnnot 'BEGIN{print tissue "\t" H_DMR_len "\t" H_overlap  "\t" annot_len/genome_len "\t" (H_overlap/H_DMR_len)/(annot_len/genome_len) "\t" H_pval "\t" "HYB" "\t" shortAnnot}' >> annotation_overlap/stats_HYB_ts_${shortAnnot}

rm C_resample_overlap
rm P_resample_overlap
rm H_resample_overlap
done

rm tmp_annot

done < "annotation.list.csv"
