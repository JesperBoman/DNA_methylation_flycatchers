#!/bin/bash -l
#SBATCH -J DMR_ts_s2
#SBATCH -o DMR_ts_s2.output
#SBATCH -e DMR_ts_s2.error
#SBATCH --mail-user jesper.boman@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 00-00:30:00
#SBATCH -A p2018002
#SBATCH -p core
#SBATCH -n 1

ml bioinfo-tools BEDTools/2.29.2
dir="./collared/tissue_specific_results"

cd $dir

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for i in "${arr[@]}"
do

A=$(ls *$i* | head -n1)
B=$(ls *$i* | tail -n3)

#Check for 25 % reciprocal overlap between the comparisons including the focal tissue
bedtools intersect -header -f 0.25 -r -wo -C -a $A -b $B  > $dir/${i}_uniq0.25fr
bedtools intersect -header -f 0.25 -r -wao -a $A -b $B -filenames > $dir/${i}_uniq0.25fr_A_B

#Print the DMRs which are significant in all comparisons including the focal tissue
awk 'NR==1{print $0} NR!=1{if($18 > 0){a[$1,$2,$3]++; b[$1,$2,$3]=$0}} END{for (dmr in a){ if(a[dmr]>=3)print b[dmr]}}' $dir/${i}_uniq0.25fr > $dir/${i}_uniq0.25fr_all4comp 


#What's hyper and hypo means for a ts-DMR is dependent on whether it was Tissue1 or Tissue2


awk -v Acomp=$A -v focal=$i 'BEGIN{A_t1t2=match(Acomp, "^" focal)} NR==1{print $0} 
NR!=1{if($17 == "."){next} else{B_t1t2=match($17, "^" focal); if( (!($1,$2,$3,$17) in c) && ((A_t1t2 == B_t1t2 && $16 == $33) || (A_t1t2 != B_t1t2 && $16 != $33))){a[$1,$2,$3]++; b[$1,$2,$3]=$0; c[$1,$2,$3,$17]}}} 
END{for (dmr in a){if(a[dmr]>=3){ print b[dmr]}}}' $dir/${i}_uniq0.25fr_A_B  > $dir/${i}_uniq0.25fr_all4comp_sameDir

awk -v Acomp=$A -v focal=$i 'BEGIN{A_t1t2=match(Acomp, "^" focal)}
NR!=1{if($17 == "."){next} else{B_t1t2=match($17, "^" focal); if( (!($1,$2,$3,$17) in c) && ((A_t1t2 == B_t1t2 && $16 == $33) || (A_t1t2 != B_t1t2 && $16 != $33))){a[$1,$2,$3]++; b[$1,$2,$3]=Acomp "\t" $0; c[$1,$2,$3,$17]}}} 
END{for (dmr in a){if(a[dmr]>=3){ print b[dmr]}}}' $dir/${i}_uniq0.25fr_A_B  | awk -v Acomp=$A -v focal=$i '{if(match(Acomp, "^" focal)==1){if($17=="hyper"){hyper++} else{hypo++}} else{if($34=="hyper"){hyper++} else{hypo++}}} END{print focal "\t" hyper/NR}' >> freq_stats


if [[ $i != "Testis" ]] ; then

A=$(ls *$i* | grep -v "Testis" | grep -v "0.25fr" | head -n1 )
B=$(ls *$i* | grep -v "Testis" | grep -v "0.25fr" | tail -n2)

bedtools intersect -header -f 0.25 -r -wao -a $A -b $B -filenames > $dir/${i}_uniq0.25fr_A_B_TExcl

awk -v Acomp=$A -v focal=$i 'BEGIN{A_t1t2=match(Acomp, "^" focal)} NR==1{print $0} 
NR!=1{if($17 == "."){next} else{B_t1t2=match($17, "^" focal); if( (!($1,$2,$3,$17) in c) && ((A_t1t2 == B_t1t2 && $16 == $33) || (A_t1t2 != B_t1t2 && $16 != $33))){a[$1,$2,$3]++; b[$1,$2,$3]=$0; c[$1,$2,$3,$17]}}} 
END{for (dmr in a){if(a[dmr]>=2){ print b[dmr]}}}'  $dir/${i}_uniq0.25fr_A_B_TExcl  > $dir/${i}_uniq0.25fr_all4comp_sameDir_TExcl
fi




done
