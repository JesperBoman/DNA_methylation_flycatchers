#!/bin/bash -l

#Create a .bed-file with genes and exons

awk '{if($3 == "gene"){match($9, /ID=(.*);Name/, a); gene[a[1]]= $1 "\t" $4-1 "\t" $5 "\t" a[1] "\t" "gene" "\t" "NA" "\t" $7 };
	  if($3 == "exon"){match($9, /ID=(.*-gene-.*)-mRNA.*:exon:/, b);  exon_count[b[1]]++;
	  	if($7 == "+"){plus_exon[b[1], exon_count[b[1]]]= $1 "\t"$4-1 "\t" $5 "\t" b[1] "\t" "exon" "\t" "exon_" exon_count[b[1]] "\t" $7 }
 	    else{minus_exon[b[1] "|" exon_count[b[1]]]=$1 "\t" $4-1 "\t" $5 "\t" b[1] "\t" "exon" "\t" "exon_" }
	  };
	  }
END{ for(i in minus_exon){split(i, c, "|"); print minus_exon[i] exon_count[c[1]] - c[2] + 1 "\t" "-"  };
	 for(i in plus_exon){print plus_exon[i]};
	 for(i in gene){print gene[i]}
	 }' <(  sort -k1,1 -k4,4n -k5,5n fAlb15.chrom_run4_featuresOnly.gff) | sort -k1,1 -k2,2n -k3,3n > genes_and_exons.bed
	  
	 
	 
#Filter for genes that have both 5' and 3' UTRs 

#Warning this code seems to work but it creates several thousand whitespace rows
# - This does not seem to be caused by the input file directly since it has no whitespace rows according to: sed '/^[[:space:]]*$/d'


awk '{if($3 == "gene" && $7 == "+"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+).*Name/, a); plus_gene[a[1]]= $1 "\t" $4-1 "\t" $5 "\t" $7 "\t" a[1]} 
  if($3 == "gene" && $7 == "-"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+).*Name/, a); minus_gene[a[1]]= $1 "\t" $4-1 "\t" $5 "\t" $7 "\t" a[1]; minus_gene_end[$1,$4]} 
  if($3 == "five_prime_UTR" && $7 == "+"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+)-mRNA-[0-9]+:/, a); plus_5utr[a[1]]} 
  if($3 == "five_prime_UTR" && $7 == "-"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+)-mRNA-[0-9]+:/, a); minus_5utr[a[1]]} 
  
  if($3 == "three_prime_UTR" && $7 == "+"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+)-mRNA-[0-9]+:/, a); plus_3utr[a[1]]} 
  if($3 == "three_prime_UTR" && $7 == "-"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+)-mRNA-[0-9]+:/, a); minus_3utr[a[1]]} 

  }
  END{for (gene in plus_gene){if( (gene in plus_5utr) && (gene in plus_3utr) ){print plus_gene[gene]}} 
  for (gene in minus_gene){if( (gene in minus_5utr) && (gene in minus_3utr) ){print minus_gene[gene]}}}' fAlb15.chrom_run4_featuresOnly.gff > fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs.bed
  
  
  
  
  
  
  
#Consider filter away genes with multiple mRNAs (not implemented yet)


##

sort -k 4,4 genes_and_exons.bed > genes_and_exons_sort.bed
sort -k 5,5  fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs.bed | sed '/^[[:space:]]*$/d' > fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs_sort.bed



wc -l fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs.bed
#Result: 8563 

wc -l fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs.bed
#Result: 8563 genes with both five and three prime UTRs in the chromosome assembly

awk 'FNR==NR{a[$5]} FNR!=NR{if($4 in a)print $0 }' fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs.bed genes_and_exons.bed | cut -f4 | sort | uniq | wc -l
#Result: 8563 unique gene IDs so code seem to work


awk 'FNR==NR{a[$5]} FNR!=NR{if($4 in a)print $0 }' fAlb15.chrom_run4_genes_with_both_five_and_three_prime_UTRs_IDs.bed genes_and_exons.bed >  genes_and_exons_5n3p.bed


#Filter problematic genes?
#Shortest is 310 bp. No need to filter now.








#Gene profile - split up intervals in 100 bp segment 

#Chop up into 100 bp pieces
ml bioinfo-tools BEDOPS/2.4.39

#Obtain 5kb flanks (obs: does not check if it is beyond the upper limit of the chromosome)

awk '{if($5 == "gene"){if(($2-5000)<0){start=0}else{start=$2-5000}; print $1 "\t" start "\t" $3+5000 "\t" $4 "\t" $5 "\t" $6 "\t" $7}}' genes_and_exons_5n3p.bed > genes_5n3p_5kb_flanks.bed



while IFS=$'\t' read -r gene_entry
do
	gene=$(echo "$gene_entry" | cut -f4)
	strand=$(echo "$gene_entry" | cut -f7)
	bedops --chop 100 <(echo "$gene_entry") > temp_gene_split
	awk -v gene="$gene" -v strand="$strand" '{print $0 "\t" gene "\t" strand "\t" "seg_" NR}' temp_gene_split >> genes_5n3p_5kb_flanks_100bp_seg.bed
done < "genes_5n3p_5kb_flanks.bed" 


#Check
cat genes_5n3p_5kb_flanks_100bp_seg.bed | cut -f4 | sort | uniq | wc -l
#8563




###INTERSECTION WITH WGBS DATA


#create a sample_list

module load bioinfo-tools BEDTools/2.29.2

bed_wgbs_data_dir="../bed_wgbs_data"

while IFS= read -r sample
do
    ind=$(basename -s ".bed.gz" $sample)
	group=$(echo $ind | cut -f1 -d "-")

#Intersection of gene features
        bedtools intersect -a <(zcat $bed_wgbs_data_dir/$sample) \
        -b genes_5n3p_5kb_flanks_100bp_seg.bed -wa -wb > temp_prof_sample.bed
        

#Output: Chromosome Start End Gene_name Strand Total_CpG Post_filter_CpG Methylated_reads_per_feature Unmethylated_reads_per_feature  Sum_of_proportions_of_methylated_over_total_reads_per_CpG Mean_per_dinuc (i.e. Methylation level) Individual Tissue

awk -v ind="$ind" -v group="$group" 'BEGIN {OFS="\t"} {if(!(annot in pos)){skip[annot]=="FALSE"}; annot=$9 "\t" $10 "\t" $11 "\t" $12; 
if(pos[annot]==$2 && $7=="CG" && prev=="CG" && NR!=1 && skip[annot]=="FALSE"){ 
        if(($5+$6+prev_total_reads[annot])>5 && ($5+$6+prev_total_reads[annot])<200){skip[annot]="TRUE";
        a[annot]+=1;  b[annot]+=($5+prev_meth_reads[annot]); c[annot]+=($6+prev_unmeth_reads[annot]); rat[annot]+=(($5+prev_meth_reads[annot])/($5+$6+prev_total_reads[annot])) }; 
        ;  strand[annot]=$13; seg_num[annot]=$14}
    else{skip[annot]="FALSE"}; 
prev=$7; if(skip[annot]=="FALSE"){pos[annot]=$3} ; 
if($7=="CG"){
        totCpG[annot]+=0.5; prev_meth_reads[annot]=$5; prev_unmeth_reads[annot]=$6; prev_total_reads[annot]=$5+$6; strand[annot]=$13; seg_num[annot]=$14}
}
    
    END{for (wind in totCpG){if (a[wind]==""){print wind, seg_num[wind], strand[wind], totCpG[wind], "NA", "NA", "NA", "NA", "NA", ind, group}
    else{print wind, seg_num[wind], strand[wind], totCpG[wind], a[wind], b[wind], c[wind], rat[wind], rat[wind]/a[wind], ind, group}}}' temp_prof_sample.bed >> meth_gene_5kb_flanks_100bp_seg_data
#Cleaning up
        rm temp_prof_sample.bed
done < "sample_list"




#Check
cat meth_gene_5kb_flanks_100bp_seg_data | cut -f4 | sort | uniq | wc -l
#8563





#PROFILE COORDINATES

#Script for positions along the gene and ordering of upstream and downstream segments

grep -v "NA" meth_gene_5kb_flanks_100bp_seg_data > meth_gene_5kb_flanks_100bp_seg_data_noNA

awk '{if($5 == "gene")print}'  genes_and_exons_5n3p.bed >  genes_5n3p.bed


#Output: Chromosome Start End Gene_name Segment_# Strand Total_CpG Post_filter_CpG Methylated_reads_per_feature Unmethylated_reads_per_feature  Sum_of_proportions_of_methylated_over_total_reads_per_CpG Mean_per_dinuc (i.e. Methylation level) Individual Tissue Segment_type Segment_number (Coordinate from 0 to 1 if a genic element)


awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)}; FNR==NR{gl[$4]=$3-$2} FNR!=NR{split($5,seg_num,"_"); if(seg_num[2]<=50){if($6 == "+"){print $0 "\t" "Upstream" "\t" seg_num[2] "\t" "NA"} else{print $0 "\t" "Downstream" "\t" 51-seg_num[2] "\t" "NA"}}
         else if(50 < seg_num[2] && (seg_num[2] < (ceil(gl[$4]/100)+51))){print $0 "\t" "Genic" "\t" seg_num[2]-50 "\t" ((seg_num[2]-50.5)*100)/gl[$4]}
         else if((ceil(gl[$4]/100)+50) < seg_num[2]){if($6 == "+"){print $0 "\t" "Downstream" "\t" seg_num[2]-(ceil(gl[$4]/100)+50) "\t" "NA" } else{print $0 "\t" "Upstream" "\t" 50-(seg_num[2]-(ceil(gl[$4]/100)+51)) "\t" "NA" }}}' genes_5n3p.bed meth_gene_5kb_flanks_100bp_seg_data_noNA > meth_gene_5kb_flanks_100bp_seg_data_noNA_coords

awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)}; FNR==NR{gl[$4]=$3-$2} FNR!=NR{split($5,seg_num,"_"); if(seg_num[2]<=50){if($6 == "+"){print $0 "\t" "Upstream" "\t" seg_num[2] "\t" "NA"} else{print $0 "\t" "Downstream" "\t" 51-seg_num[2] "\t" "NA"}}
         else if(50 < seg_num[2] && (seg_num[2] < (ceil(gl[$4]/100)+51))){if($6 == "+"){print $0 "\t" "Genic" "\t" seg_num[2]-50 "\t" ((seg_num[2]-50.5)*100)/gl[$4]} else{print $0 "\t" "Genic" "\t" seg_num[2]-50 "\t" 1-(((seg_num[2]-50.5)*100)/gl[$4])} }
         else if((ceil(gl[$4]/100)+50) < seg_num[2]){if($6 == "+"){print $0 "\t" "Downstream" "\t" seg_num[2]-(ceil(gl[$4]/100)+50) "\t" "NA" } else{print $0 "\t" "Upstream" "\t" 50-(seg_num[2]-(ceil(gl[$4]/100)+51)) "\t" "NA" }}}' genes_5n3p.bed meth_gene_5kb_flanks_100bp_seg_data_noNA > meth_gene_5kb_flanks_100bp_seg_data_noNA_coords





#CGI ANNOTATION

module load bioinfo-tools BEDTools/2.29.2



bedtools intersect -a genes_5n3p_5kb_flanks.bed -b ../Meth_annot/CGI_Carina_liftedTo_chr_sort.bed -wa -wb > genes_5n3p_5kb_flanks_CGI_overlap.bed
bedtools intersect -a genes_5n3p_5kb_flanks_100bp_seg.bed -b ../Meth_annot/CGI_Carina_liftedTo_chr_sort.bed -wa -wb > genes_5n3p_5kb_flanks_100bp_seg_CGI_overlap.bed


awk 'FNR==NR{num[$1 "\t" $2 "\t" $3]++; oe[$1 "\t" $2 "\t" $3]+= ($12*$10); len[$1 "\t" $2 "\t" $3]+=$10 } FNR!=NR{prom=$1 "\t" $2 "\t" $3; if(prom in num){print $0 "\t" num[prom] "\t" len[prom] "\t" oe[prom]/len[prom] "\t" "Y" } else{print $0 "\t" "0" "\t" "NA" "\t" "NA" "\t" "N"}}' genes_5n3p_5kb_flanks_100bp_seg_CGI_overlap.bed meth_gene_5kb_flanks_100bp_seg_data_coords > meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot





#### ADD GENE EXPRESSION ESTIMATES ####

tissues="Heart Kidney Brain Testis Liver"

mkdir Expression_results_annotated

for tissue in $tissues
do
awk '{if(FNR==NR){a[$4]= $1 "_" $2 "_" $3}; if(FNR==1 && NR>FNR){print $0 "\t" "region_ID"}; if(FNR!=1 && NR>FNR && ($1 in a)){print $0 "\t" a[$1]}}' genes.all.run4.geneNames.bed /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/rnaseq_pipe/${tissue}/star_salmon/salmon.merged.gene_tpm_length_scaled.tsv > Expression_results_annotated/${tissue}_salmon.merged.gene_tpm_length_scaled.regionID.tsv


done
wait


#Add region ID (i.e. gene coordinates here to methylation data file)
awk 'FNR==NR{a[$4]=$1 "_" $2 "_" $3} FNR!=NR{if($4 in a){print $0 "\t" a[$4]}}' genes_5n3p.bed meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot > meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID

#


meth_table="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID"
tissues="Heart Kidney Brain Testis Liver"

for tissue in $tissues
do

grep $tissue $meth_table > temp_${tissue}


awk ' FNR==NR{if(NR!=1){COL01[$16]=$2; COL02[$16]=$3; COL03[$16]=$4; COL04[$16]=$5; COL05[$16]=$6; HYB04[$16]=$7; 
	HYB01[$16]=$8; HYB02[$16]=$9; HYB05[$16]=$10; 
	PIE01[$16]=$11; PIE02[$16]=$12; PIE03[$16]=$13; PIE04[$16]=$14;  PIE05[$16]=$15;}} FNR!=NR{ if($13 == "COL01"){print $0 "\t" COL01[$22]}
	if($13 == "COL02"){print $0 "\t" COL02[$22]}
	if($13 == "COL03"){print $0 "\t" COL03[$22]}
	if($13 == "COL04"){print $0 "\t" COL04[$22]}
	if($13 == "COL05"){print $0 "\t" COL05[$22]}
	if($13 == "HYB04"){print $0 "\t" HYB04[$22]}
	if($13 == "HYB01"){print $0 "\t" HYB01[$22]}
	if($13 == "HYB02"){print $0 "\t" HYB02[$22]}
	if($13 == "HYB05"){print $0 "\t" HYB05[$22]}
	if($13 == "PIE01"){print $0 "\t" PIE01[$22]}
	if($13 == "PIE02"){print $0 "\t" PIE02[$22]}
	if($13 == "PIE03"){print $0 "\t" PIE03[$22]}
	if($13 == "PIE04"){print $0 "\t" PIE04[$22]}
	if($13 == "PIE05"){print $0 "\t" PIE05[$22]}}' Expression_results_annotated/${tissue}_salmon.merged.gene_tpm_length_scaled.regionID.tsv temp_${tissue} >> meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr

done




#### Expression inheritance pattern DNA methylation gene profile code ####

Read in conk.df - use promoter IDs and Tissue as key

Reconstruct promoter IDs from gene IDs in meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr  


awk 'FNR == NR && FNR>1 {inpatExpr[$1 "\t" $2]=$7 } FNR!=NR {split($22, a, "_"); 
	if($6 == "+"){regID=a[1] "_" a[2]-2000 "_" a[2]+1 "\t" $14 };
	if($6 == "-"){regID=a[1] "_" a[3]-1 "_" a[3]+2000 "\t" $14 }; 
	if(regID in inpatExpr){print $0 "\t" inpatExpr[regID]} else{print $0 "\t" "NA" }}' conk.df.txt meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr > meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr






#### DMR presence/absence gene profile ####

module load bioinfo-tools BEDTools/2.29.2

dir="../DMRs"

awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvP_brain_dmrs.txt > CvP_brain_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvP_heart_dmrs.txt > CvP_heart_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $1,$2-1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' CvP_kidney_dmrs.txt > CvP_kidney_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvP_liver_dmrs.txt > CvP_liver_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvP_testis_dmrs.txt > CvP_testis_dmrs.bed

bedtools intersect -C -a genes_5n3p_5kb_flanks_100bp_seg.bed -b $dir/brain/CvP_brain_dmrs.bed $dir/heart/CvP_heart_dmrs.bed $dir/kidney/CvP_kidney_dmrs.bed $dir/liver/CvP_liver_dmrs.bed $dir/testis/CvP_testis_dmrs.bed  > genes_5n3p_5kb_flanks_100bp_seg_CvP_DMR_overlap.bed



awk 'BEGIN{tissue[1]="Brain"; tissue[2]="Heart"; tissue[3]="Kidney"; tissue[4]="Liver"; tissue[5]="Testis"}
 FNR==NR{if($8==0){dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=0} else{dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=1} } 
 FNR!=NR{segtiss=$1 "\t" $2 "\t" $3 "\t" $14; if(segtiss in dmr){print $0 "\t" dmr[segtiss]} else{print $0 "\t" "0"}}' genes_5n3p_5kb_flanks_100bp_seg_CvP_DMR_overlap.bed meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_HYB01_Liver >meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_HYB01_Liver

awk 'BEGIN{tissue[1]="Brain"; tissue[2]="Heart"; tissue[3]="Kidney"; tissue[4]="Liver"; tissue[5]="Testis"}
 FNR==NR{if($8==0){dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=0} else{dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=1} } 
 FNR!=NR{segtiss=$1 "\t" $2 "\t" $3 "\t" $14; if(segtiss in dmr){print $0 "\t" dmr[segtiss]} else{print $0 "\t" "0"}}' genes_5n3p_5kb_flanks_100bp_seg_CvP_DMR_overlap.bed meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr >meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR




#### DEseq CvP annot ####
module load bioinfo-tools R_packages/4.0.4

#Part 1 - Generate region-annotated expression tables

tissues="Heart Kidney Brain Testis Liver"

dir="."
#mkdir Expression_results_annotated

skip="N"

if [[ $skip == "Y" ]]; then
        echo "Skipping part 1 - annotating DeSeq2 output"

else

        for tissue in $tissues
        do
        inputDir="../rnaseq_pipe/${tissue}/star_salmon/deseq2_qc"
        cd $inputDir
        Rscript --vanilla $dir/DESeq_code_COLvPIE.R $tissue
        cd $dir

        awk '{if(FNR==NR){a[$4]= $1 "_" $2 "_" $3}; if(FNR!=1 && NR>FNR && ($1 in a)){print $3 "\t" $7 "\t" a[$1]}}' genes.all.run4.geneNames.bed $inputDir/${tissue}.diff.expr.res.COLvPIE.df.txt > Expression_results_annotated/${tissue}.diff.expr.res.COLvPIE.df.regionID.txt


done
fi

#Part 2 - Annotate meth table with region ID

#awk 'FNR==NR{a[$4]=$1 "_" $2 "_" $3} FNR!=NR{if($4 in a){print $0 "\t" a[$4]}}' genes_5n3p.bed meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot > meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID

#Part 3 - Annotate meth table with expression data

meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR"
tissues="Heart Kidney Brain Testis Liver"

for tissue in $tissues
do

grep $tissue $meth_table > temp_${tissue}


awk ' FNR==NR{stats[$3]=$1 "\t" $2} FNR!=NR{if($22 in stats){print $0 "\t" stats[$22]} else{print $0 "\t" "NA" "\t" "NA"}}' Expression_results_annotated/${tissue}.diff.expr.res.COLvPIE.df.regionID.txt temp_${tissue} >> ${meth_table}_l2fcCOLvPIE
rm temp_${tissue}
done



#### MethCOL mean, MethPIE mean, MethHYB mean ####

meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE"
tissues="Heart Kidney Brain Testis Liver"

for tissue in $tissues
do

grep $tissue $meth_table > temp_${tissue}


awk 'FNR==NR{if($13 ~ /COL/ || $13 == "HYB04"){methCOL[$1 "\t" $2 "\t" $3]+= $12; iCOL[$1 "\t" $2 "\t" $3]++} else if($13 ~ /PIE/ ){methPIE[$1 "\t" $2 "\t" $3]+= $12; iPIE[$1 "\t" $2 "\t" $3]++} else if($13 ~ /HYB/ ){methHYB[$1 "\t" $2 "\t" $3]+= $12; iHYB[$1 "\t" $2 "\t" $3]++}} 
	FNR!=NR{if(iCOL[$1 "\t" $2 "\t" $3] != 0 && iPIE[$1 "\t" $2 "\t" $3] != 0 && iHYB[$1 "\t" $2 "\t" $3] != 0){print $0 "\t" methCOL[$1 "\t" $2 "\t" $3]/iCOL[$1 "\t" $2 "\t" $3] "\t" methPIE[$1 "\t" $2 "\t" $3]/iPIE[$1 "\t" $2 "\t" $3] "\t" methHYB[$1 "\t" $2 "\t" $3]/iHYB[$1 "\t" $2 "\t" $3]}
	else{ print $0 "\t" "NA" "\t" "NA" "\t" "NA"}}' temp_${tissue} temp_${tissue} >> meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME


done















#### Fst gene profile ####

module load bioinfo-tools BEDTools/2.29.2

dir="../popgen"


#I'm was here using the weighted mean Fst, now I changed it to unweighted mean Fst, see Jackson et al. 2015 Heredity for a discussion of which to use
bedtools intersect -f 0.51 -r -wa -wb -a genes_5n3p_5kb_flanks_100bp_seg.bed -b $dir/CpG.windowed.weir.fst.chrom.strand.corr.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $12 }' > genes_5n3p_5kb_flanks_100bp_seg_CpGFst_overlap.bed
bedtools intersect -f 0.51 -r -wa -wb -a genes_5n3p_5kb_flanks_100bp_seg.bed -b $dir/non-CpG.windowed.weir.fst.chrom.strand.corr.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $12 }' > genes_5n3p_5kb_flanks_100bp_seg_non-CpGFst_overlap.bed

meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME"


awk 'FNR==NR{fst[$1 "\t" $2 "\t" $3]= $4 "\t" $5} 
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in fst){print $0 "\t" fst[seg]} else{print $0 "\t" NA" "\t" "NA"}}' genes_5n3p_5kb_flanks_100bp_seg_CpGFst_overlap.bed $meth_table  > ${meth_table}_fst_annot_tmp

awk 'FNR==NR{fst[$1 "\t" $2 "\t" $3]= $4 "\t" $5} 
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in fst){print $0 "\t" fst[seg]} else{print $0 "\t" "NA" "\t" "NA"}}' genes_5n3p_5kb_flanks_100bp_seg_non-CpGFst_overlap.bed  ${meth_table}_fst_annot_tmp > ${meth_table}_fst_annot
 
 rm ${meth_table}_fst_annot_tmp
 


 
 
 ## GC and CpG OE gene profile ###
 
module load bioinfo-tools BEDTools/2.29.2

genome="../methylseq_pipe/fAlb15.chrom_Olandfix_polyCorG100bpthres_hardmasked.fa"

bedtools nuc -bed genes_5n3p_5kb_flanks_100bp_seg.bed -fi $genome | awk 'NR>1{if(($9+$10+$11+$12) != 0){print $1 "\t" $2 "\t" $3 "\t"  ($10+$11)/($9+$10+$11+$12)} else{print $1 "\t" $2 "\t" $3 "\t" "NA"}}' > genes_5n3p_5kb_flanks_100bp_seg_GC_content.bed


bedtools nuc -bed genes_5n3p_5kb_flanks_100bp_seg.bed -fi $genome -pattern CG | awk 'NR>1{totbases=$9+$10+$11+$12; if(totbases != 0 && ($10 < 95) && ($11 < 95) ){print $1 "\t" $2 "\t" $3 "\t"  ($10+$11)/totbases "\t" ($10*$11)/(totbases) "\t" $16} else{print $1 "\t" $2 "\t" $3 "\t" "NA" "\t" "NA" "\t" "NA"}}' > genes_5n3p_5kb_flanks_100bp_seg_GC_content_CpGOE_.bed




meth_table="meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot"


awk 'FNR==NR{GCetc[$1 "\t" $2 "\t" $3]= $4 "\t" $5 "\t" $6} 
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in GCetc){print $0 "\t" GCetc[seg]} else{print $0 "\t" "NA" "\t" "NA" "\t" "NA"}}' genes_5n3p_5kb_flanks_100bp_seg_GC_content_CpGOE.bed $meth_table  > ${meth_table}_GCetc

 


awk 'FNR==NR{GCetc[$1 "\t" $2 "\t" $3]= $4 "\t" $5 "\t" $6} 
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in GCetc){print $0 "\t" GCetc[seg]} else{print $0 "\t" "NA" "\t" "NA" "\t" "NA"}}' genes_5n3p_5kb_flanks_100bp_seg_GC_content_CpGOE.bed meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_HYB02_Heart_Fst  > meth_gene_5kb_flanks_100bp_seg_data_noNA_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_HYB02_Heart_Fst_GCetc




# Exon/Intron annotation #


module load bioinfo-tools BEDTools/2.29.2

awk '{if($5 == "exon")print $0}' genes_and_exons_5n3p.bed > exons_5n3p.bed

bedtools subtract -a genes_5n3p.bed -b exons_5n3p.bed > introns_raw_5n3p.bed


sort -k1,1 -k2,2n -k3,3n introns_raw_5n3p.bed > introns_raw_5n3p_sort.bed

awk 'NR==FNR{gene_strand[$4]=$7; intron_num[$4]++} NR!=FNR{if ($7 == "+"){ a[$4]+=1; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "intron" "\t" "intron_"a[$4] "\t" $7}
	else{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "intron" "\t" "intron_"intron_num[$4]-a[$4] "\t" $7; a[$4]+=1}}' introns_raw_5n3p_sort.bed introns_raw_5n3p_sort.bed > intron_intronNumbers_5n3p.bed



segment="genes_5n3p_5kb_flanks_100bp_seg.bed"

bedtools intersect -wo -f 0.9 -a $segment -b intron_intronNumbers_5n3p.bed  > genes_5n3p_5kb_flanks_100bp_seg_intronOverlap.bed
bedtools intersect -wo -f 0.9 -a $segment -b exons_5n3p.bed > genes_5n3p_5kb_flanks_100bp_seg_exonOverlap.bed


meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot"

awk 'FNR==NR{if($11 == "intron"){intSeg[$1 "\t" $2 "\t" $3]= $12} else{exSeg[$1 "\t" $2 "\t" $3]= $12}}
 FNR!=NR{seg=$1 "\t" $2 "\t" $3; if(seg in intSeg){print $0 "\t" intSeg[seg]} else if(seg in exSeg){print $0 "\t" exSeg[seg]} else{print $0 "\t" "NA"}}' <(cat genes_5n3p_5kb_flanks_100bp_seg_intronOverlap.bed genes_5n3p_5kb_flanks_100bp_seg_exonOverlap.bed)  $meth_table > ${meth_table}_ExIn



#### Hybrid DMR presence/absence gene profile ####


module load bioinfo-tools BEDTools/2.29.2

dir="../DMRs"

awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvH_brain_dmrs.txt > CvH_brain_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvH_heart_dmrs.txt > CvH_heart_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $1,$2-1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' CvH_kidney_dmrs.txt > CvH_kidney_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvH_liver_dmrs.txt > CvH_liver_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' CvH_testis_dmrs.txt > CvH_testis_dmrs.bed

awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' PvH_brain_dmrs.txt > PvH_brain_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' PvH_heart_dmrs.txt > PvH_heart_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $1,$2-1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' PvH_kidney_dmrs.txt > PvH_kidney_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' PvH_liver_dmrs.txt > PvH_liver_dmrs.bed
awk 'NR==1{print $0} NR>1{OFS="\t"; print $2,$3-1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' PvH_testis_dmrs.txt > PvH_testis_dmrs.bed

bedtools intersect -C -a genes_5n3p_5kb_flanks_100bp_seg.bed -b $dir/brain/CvH_brain_dmrs.bed $dir/heart/CvH_heart_dmrs.bed $dir/kidney/CvH_kidney_dmrs.bed $dir/liver/CvH_liver_dmrs.bed $dir/testis/CvH_testis_dmrs.bed  > genes_5n3p_5kb_flanks_100bp_seg_CvH_DMR_overlap.bed
bedtools intersect -C -a genes_5n3p_5kb_flanks_100bp_seg.bed -b $dir/brain/PvH_brain_dmrs.bed $dir/heart/PvH_heart_dmrs.bed $dir/kidney/PvH_kidney_dmrs.bed $dir/liver/PvH_liver_dmrs.bed $dir/testis/PvH_testis_dmrs.bed  > genes_5n3p_5kb_flanks_100bp_seg_PvH_DMR_overlap.bed

meth_table="meth_gene_5kb_flanks_100bp_seg_data_coords_CGI_annot_regionID_Expr_InpatExpr_CvPDMR_l2fcCOLvPIE_GEnME_fst_annot_ExIn"

awk 'BEGIN{tissue[1]="Brain"; tissue[2]="Heart"; tissue[3]="Kidney"; tissue[4]="Liver"; tissue[5]="Testis"}
 FNR==NR{if($8==0){dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=0} else{dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=1} } 
 FNR!=NR{segtiss=$1 "\t" $2 "\t" $3 "\t" $14; if(segtiss in dmr){print $0 "\t" dmr[segtiss]} else{print $0 "\t" "0"}}' genes_5n3p_5kb_flanks_100bp_seg_CvH_DMR_overlap.bed $meth_table >$meth_table.tmp

awk 'BEGIN{tissue[1]="Brain"; tissue[2]="Heart"; tissue[3]="Kidney"; tissue[4]="Liver"; tissue[5]="Testis"}
 FNR==NR{if($8==0){dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=0} else{dmr[$1 "\t" $2 "\t" $3 "\t" tissue[$7]]=1} } 
 FNR!=NR{segtiss=$1 "\t" $2 "\t" $3 "\t" $14; if(segtiss in dmr){print $0 "\t" dmr[segtiss]} else{print $0 "\t" "0"}}' genes_5n3p_5kb_flanks_100bp_seg_PvH_DMR_overlap.bed $meth_table.tmp >${meth_table}_HDMR
 
 rm $meth_table.tmp
