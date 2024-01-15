#!/bin/bash -l





ml bioinfo-tools BEDTools/2.29.2

#awk '{print $1 "\t" $2-101 "\t" $2+100 "\t" $3 "\t" $4}' fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel > fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_plus_minus_100bp.bed

dir="Bismark_output"

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for tissue in "${arr[@]}"
do

for i in 1 2 5
do 

ls "$dir" | grep "HYB0${i}-${tissue}" | grep "CX_report" > list_tmp

while IFS= read -r CX_report
do
spec=$(echo "$CX_report" | cut -d \. -f3 | sed 's/genome1/COL/' | sed 's/genome2/PIE/')

#Filtering per cov 3 version
#zcat $dir/$CX_report | awk '{if($3+$4 < 3){next} else{OFS="\t"; if($6 == "CG"){print $1, $2-1, $2, $3, $4, $5, $6, $7 }}}' > "${dir}/HYB0${i}-${tissue}_${spec}_cov3_CpG.bed"

zcat $dir/$CX_report | awk '{OFS="\t"; if($6 == "CG"){print $1, $2-1, $2, $3, $4, $5, $6, $7 }}' > "${dir}/HYB0${i}-${tissue}_${spec}_CpG.bed"



bedtools intersect -a "${dir}/HYB0${i}-${tissue}_${spec}_CpG.bed" -b fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_plus_minus_100bp.bed -wa -wb > "${dir}/HYB0${i}-${tissue}_${spec}_CpG_100bpPMfd.bed"


awk -v ind="HYB0${i}" -v tissue="$tissue" -v parental="$spec" 'BEGIN {OFS="\t"; skip="FALSE"} {annot=$9 "\t" $10 "\t" $11; 
if(pos==$2 && $7=="CG" && prev=="CG" && NR!=1 && skip=="FALSE"){ 
        if(($5+$6+prev_total_reads)>2 && ($5+$6+prev_total_reads)<200){skip="TRUE";
        a[annot]+=1;  b[annot]+=($5+prev_meth_reads); c[annot]+=($6+prev_unmeth_reads); rat[annot]+=(($5+prev_meth_reads)/($5+$6+prev_total_reads)) }; 
    chr[annot]=$9; start[annot]=$10; end[annot]=$11}
    else{skip="FALSE"}; 
prev=$7; if(skip=="FALSE"){pos=$3} ; 
if($7=="CG"){
        totCpG[annot]+=0.5; prev_meth_reads=$5; prev_unmeth_reads=$6; prev_total_reads=$5+$6}
}
    
    END{for (wind in chr){if (a[wind]==""){print wind, "NA", "NA", "NA", "NA", "NA", "NA", ind, tissue, "Hybrid-" parental}
    else{print wind, totCpG[wind], a[wind], b[wind], c[wind], rat[wind], rat[wind]/a[wind], ind, tissue, "Hybrid-"  parental}}}' "${dir}/HYB0${i}-${tissue}_${spec}_CpG_100bpPMfd.bed"  >> meth_asm_snp_hyb_results


done < "list_tmp"

done
done
