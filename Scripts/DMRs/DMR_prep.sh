#!/bin/bash -l

#A script to to turn .bed-formatted stranded CpG reports from Bismark into the Bsseq format to call DMRs using the BSmooth method

dir="../bed_wgbs_data"



ls "$dir" > sample_list_temp

while IFS= read -r sample
do


zcat $dir/$sample | awk '{if($7 == "CG" && prevPos == $2-1){CpG++} else{CpG=0}; if($7 == "CG"  && prevPos == $2-1 && $1 == prevChr && CpG == 2 &&  ($5+$6+prevM+prevC) < 200 ){CpG=0; print $1 "\t" $3 "\t" $5+prevM "\t" $5+$6+prevM+prevC}; if($7 == "CG"){prevM=$5; prevC=$6; prevChr=$1; prevPos=$2; CpG++}}' > bsseq_formatted_data/${sample}_data_formatted_temp &

p=$(($p+1))
remainder=$(( p  % 5 ))
if [[ $remainder -eq 0 ]] ; then
wait
fi


done <"sample_list_temp"

wait
