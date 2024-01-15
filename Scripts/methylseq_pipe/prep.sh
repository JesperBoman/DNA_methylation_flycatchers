#!/bin/bash -l


declare -a spec=("COL" "HYB" "PIE")
declare -a tiss=("Heart" "Kidney" "Liver" "Brain" "Testis")
dir="/PATH/TO/FASTQ"

for species in "${spec[@]}"
do

for tissue in "${tiss[@]}"
do

for i in 1 2 3 4 5
do

cat $dir/${species}0${i}-${tissue}*R1*fastq.gz > tmp_fastq/tmp_${species}0${i}-${tissue}_R1.fastq.gz &
cat $dir/${species}0${i}-${tissue}*R2*fastq.gz > tmp_fastq/tmp_${species}0${i}-${tissue}_R2.fastq.gz &

#echo $dir/${species}0${i}-${tissue}*R1*fastq.gz 
#echo $dir/${species}0${i}-${tissue}*R2*fastq.gz
done
		
if [ $(jobs -p | wc -w) -ge 20 ]; then
wait
fi

done
done
wait
