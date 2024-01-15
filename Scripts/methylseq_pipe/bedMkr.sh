#!/bin/bash -l

#Transform stranded CpG format into .bed-files

declare -a arr=("Liver")

for tissue in "${arr[@]}"; do

dir="/results/${tissue}/bismark_methylation_calls/stranded_CpG_report"
ls "$dir" > sample_list

	while IFS= read -r samp
	do
		shortname=$(echo $samp | cut -f 1 -d "_")
		zcat $dir/${samp} | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}'  > ${shortname}.bed
		gzip ${shortname}.bed	
		mv ${shortname}.bed.gz ../bed_wgbs_data
	done < "sample_list"
done
