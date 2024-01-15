#!/bin/bash -l




ml bioinfo-tools BEDTools/2.29.2

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

dir="../bed_wgbs_data" #Directory where the methylation calls were stored in .bed-format

for tissue in "${arr[@]}"
do

ls "$dir" | grep "$tissue" > sample_list

while IFS= read -r sample
do
	ind=$(echo $sample | cut -d- -f1)
	tissue=$(echo $sample | cut -d. -f1 | cut -d- -f2)

bedtools intersect -a <(zcat $dir/$sample) -b fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_plus_minus_100bp.bed -wa -wb > temp_${sample}.bed

if [ "$sample" != "HYB04-${tissue}.bed.gz" ] ; then
	spec=$(echo $sample | cut -d 0 -f1)

else
	spec="COL"
fi

awk -v ind="$ind" -v tissue="$tissue" -v parental="$spec" 'BEGIN {OFS="\t"; skip="FALSE"} {annot=$9 "\t" $10 "\t" $11; 
if(pos==$2 && $7=="CG" && prev=="CG" && NR!=1 && skip=="FALSE"){ 
        if(($5+$6+prev_total_reads)>2 && ($5+$6+prev_total_reads)<200){skip="TRUE";
        a[annot]+=1;  b[annot]+=($5+prev_meth_reads); c[annot]+=($6+prev_unmeth_reads); rat[annot]+=(($5+prev_meth_reads)/($5+$6+prev_total_reads)) }; 
    chr[annot]=$9; start[annot]=$10; end[annot]=$11}
    else{skip="FALSE"}; 
prev=$7; if(skip=="FALSE"){pos=$3} ; 
if($7=="CG"){
        totCpG[annot]+=0.5; prev_meth_reads=$5; prev_unmeth_reads=$6; prev_total_reads=$5+$6}
}
    
    END{for (wind in chr){if (a[wind]==""){print wind, "NA", "NA", "NA", "NA", "NA", "NA", ind, tissue, "Parental-" parental}
    else{print wind, totCpG[wind], a[wind], b[wind], c[wind], rat[wind], rat[wind]/a[wind], ind, tissue, "Parental-"  parental}}}' "temp_${sample}.bed"  >> meth_asm_snp_par_results

rm temp_${sample}.bed

done < "sample_list"
done
