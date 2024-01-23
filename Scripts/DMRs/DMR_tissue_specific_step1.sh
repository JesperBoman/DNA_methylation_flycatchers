#!/bin/bash -l

ml bioinfo-tools BEDTools/2.29.2




#Obtaining tissue-specific DMRs

#Pipeline 1: This basically answers the question: what DMRs are unique for a certain comparison, compared to other comparisons where the focal tissue (Tissue1) and compared tissue (Tissue2) is not included.

dir="DMRs/collared"

#This is the order the comparisons been made. Keep this intact. Heart is the highest in a descending order of tissue comparisons.
declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")
declare -a arr2=("Heart" "Kidney" "Liver" "Brain" "Testis")


p=0
for foc in "${arr[@]}"
do

#Remove the focal tissue from this array
unset 'arr2[$p]'
p=$(($p+1))
k=0

for i in "${arr2[@]}"
do

#Enter focal tissue here
Tissue1=$foc
Tissue2=$i

awk 'NR==1 {print $0} NR>1{for(i=1; i<=NF; i++){ if(i==1){printf $i} else if(i==2){printf "\t" $2 - 1}  else{printf "\t" $i}} printf "\n" }' $dir/${Tissue1}v${Tissue2}/dmrs.txt > tmpA

echo "${Tissue1}v${Tissue2} is being compared to: " 
for f in $(seq 0 4)
do

if [[ "$Tissue1" != "${arr[$f]}" ]] && [[ "$Tissue2" != "${arr[$f]}" ]]; then

Tissue3=${arr[$f]}

Tissue4=${arr[$(($f+1))]}
if [[ ! -z $Tissue4 ]] && [[ $Tissue1 != $Tissue4 ]] && [[ $Tissue2 != $Tissue4 ]] ; then
echo "${Tissue3}v${Tissue4}"
awk 'NR==1 {print $0} NR>1{for(i=1; i<=NF; i++){ if(i==1){printf $i} else if(i==2){printf "\t" $2 - 1}  else{printf "\t" $i}} printf "\n" }' $dir/${Tissue3}v${Tissue4}/dmrs.txt > tmpB${k}
k=$(($k+1))
fi

Tissue4=${arr[$(($f+2))]}
if [[ ! -z $Tissue4 ]] && [[ $Tissue1 != $Tissue4 ]] && [[ $Tissue2 != $Tissue4 ]] ; then
echo "${Tissue3}v${Tissue4}"
awk 'NR==1 {print $0} NR>1{for(i=1; i<=NF; i++){ if(i==1){printf $i} else if(i==2){printf "\t" $2 - 1}  else{printf "\t" $i}} printf "\n" }' $dir/${Tissue3}v${Tissue4}/dmrs.txt > tmpB${k}
k=$(($k+1))
fi

Tissue4=${arr[$(($f+3))]}
if [[ ! -z $Tissue4 ]] && [[ $Tissue1 != $Tissue4 ]] && [[ $Tissue2 != $Tissue4 ]] ; then
echo "${Tissue3}v${Tissue4}"
awk 'NR==1 {print $0} NR>1{for(i=1; i<=NF; i++){ if(i==1){printf $i} else if(i==2){printf "\t" $2 - 1}  else{printf "\t" $i}} printf "\n" }' $dir/${Tissue3}v${Tissue4}/dmrs.txt > tmpB${k}
k=$(($k+1))
fi

Tissue4=${arr[$(($f+4))]}
if [[ ! -z $Tissue4 ]] && [[ $Tissue1 != $Tissue4 ]] && [[ $Tissue2 != $Tissue4 ]] ; then
echo "${Tissue3}v${Tissue4}"
awk 'NR==1 {print $0} NR>1{for(i=1; i<=NF; i++){ if(i==1){printf $i} else if(i==2){printf "\t" $2 - 1}  else{printf "\t" $i}} printf "\n" }' $dir/${Tissue3}v${Tissue4}/dmrs.txt > tmpB${k}
k=$(($k+1))
fi


fi
done        

bedtools intersect -header -f 0.25 -v -a tmpA -b tmpB*  > $dir/tissue_specific_results/${Tissue1}v${Tissue2}_uniq0.25f

rm tmpB*
k=0
done
done
