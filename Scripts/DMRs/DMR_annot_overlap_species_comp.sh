#!/bin/bash -l

## DMR overlap analysis ##

#mkdir annotation_overlap

ml bioinfo-tools BEDTools/2.29.2
faidx="../methylseq_pipe/fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa.fai"
dir="."


while IFS= read -r annot_double
do

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annot=$(cut -f2 <(echo $annot_double) -d ",")

#Bedtools merge here is a convenient heuristic since annotation features can overlap but are bound to a specific gene such as is the case for UTRs
#Genes can also overlap and it is probably fine to consider them as one long gene region instead of counting them as a double overlap.
awk '{print $1 "\t" $2 "\t" $3}' $annot | sort -k1,1 -k2,2n | bedtools merge > tmp_annot_sp
annot=tmp_annot_sp

annot_len=$(awk '{sum+=($3-$2)} END{print sum}' $annot)
genome_len=$(awk '{sum+=$2} END{print sum}' $faidx)


declare -a tissues=("heart" "kidney" "liver" "brain" "testis")

RANDOM=$(date +%N | cut -b4-9)

for tissue in "${tissues[@]}"
do

CvP_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dir/$tissue/CvP_${tissue}_dmrs.bed)
CvH_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dir/$tissue/CvH_${tissue}_dmrs.bed)
PvH_DMR_len=$(awk 'NR>1{sum+=($3-$2)} END{print sum}' $dir/$tissue/PvH_${tissue}_dmrs.bed)

CvP_overlap=$(bedtools intersect -header -wo -a $dir/$tissue/CvP_${tissue}_dmrs.bed -b $annot | cut -f20 | awk '{sum+=$1}END{print sum}')
CvH_overlap=$(bedtools intersect -header -wo -a $dir/$tissue/CvH_${tissue}_dmrs.bed -b $annot | cut -f20 | awk '{sum+=$1}END{print sum}')
PvH_overlap=$(bedtools intersect -header -wo -a $dir/$tissue/PvH_${tissue}_dmrs.bed -b $annot | cut -f20 | awk '{sum+=$1}END{print sum}')

resamples=1000

for i in $(seq 1 $resamples)
do

bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/$tissue/CvP_${tissue}_dmrs.bed ) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> CvP_resample_overlap &

bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/$tissue/CvH_${tissue}_dmrs.bed) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> CvH_resample_overlap &

bedtools intersect -header -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/$tissue/PvH_${tissue}_dmrs.bed) \
	-b $annot | cut -f7 | awk '{sum+=$1}END{print sum}' >> PvH_resample_overlap &
	wait
done

CvP_ln=$(cat CvP_resample_overlap <(echo $CvP_overlap) | sort -n | grep -n -w "$CvP_overlap" | cut -f1 -d ":" | head -n1 )
CvH_ln=$(cat CvH_resample_overlap <(echo $CvH_overlap) | sort -n | grep -n -w "$CvH_overlap" | cut -f1 -d ":" | head -n1 )
PvH_ln=$(cat PvH_resample_overlap <(echo $PvH_overlap) | sort -n | grep -n -w "$PvH_overlap" | cut -f1 -d ":" | head -n1 )


CvP_pval=$(awk -v CvP_ln=$CvP_ln -v resamples=$resamples 'BEGIN{r=resamples+1-CvP_ln; print r/(resamples)}')
CvH_pval=$(awk -v CvH_ln=$CvH_ln -v resamples=$resamples 'BEGIN{r=resamples+1-CvH_ln; print r/(resamples)}')
PvH_pval=$(awk -v PvH_ln=$PvH_ln -v resamples=$resamples 'BEGIN{r=resamples+1-PvH_ln; print r/(resamples)}')


#Count, mean number of CpGs, mean width, meanDiff (abs value)
awk -v tissue=$tissue -v annot_len=$annot_len -v genome_len=$genome_len -v CvP_DMR_len=$CvP_DMR_len -v CvP_overlap=$CvP_overlap -v CvP_pval=$CvP_pval -v shortAnnot=$shortAnnot 'BEGIN{print tissue "\t" CvP_DMR_len "\t" CvP_overlap  "\t" annot_len/genome_len "\t" (CvP_overlap/CvP_DMR_len)/(annot_len/genome_len) "\t" CvP_pval "\t" "CvP" "\t" shortAnnot}' >> annotation_overlap/stats_CvP_${shortAnnot}
awk -v tissue=$tissue -v annot_len=$annot_len -v genome_len=$genome_len -v CvH_DMR_len=$CvH_DMR_len -v CvH_overlap=$CvH_overlap -v CvH_pval=$CvH_pval -v shortAnnot=$shortAnnot 'BEGIN{print tissue "\t" CvH_DMR_len "\t" CvH_overlap  "\t" annot_len/genome_len "\t" (CvH_overlap/CvH_DMR_len)/(annot_len/genome_len) "\t" CvH_pval "\t" "CvH" "\t" shortAnnot}' >> annotation_overlap/stats_CvH_${shortAnnot}
awk -v tissue=$tissue -v annot_len=$annot_len -v genome_len=$genome_len -v PvH_DMR_len=$PvH_DMR_len -v PvH_overlap=$PvH_overlap -v PvH_pval=$PvH_pval -v shortAnnot=$shortAnnot 'BEGIN{print tissue "\t" PvH_DMR_len "\t" PvH_overlap  "\t" annot_len/genome_len "\t" (PvH_overlap/PvH_DMR_len)/(annot_len/genome_len) "\t" PvH_pval "\t" "PvH" "\t" shortAnnot}' >> annotation_overlap/stats_PvH_${shortAnnot}


rm CvP_resample_overlap
rm CvH_resample_overlap
rm PvH_resample_overlap
done

rm tmp_annot_sp

done < "annotation.list.csv"
