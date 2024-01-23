#!/bin/bash -l

#translateChromToScaff.pl and translateScaffToChrom.pl are scripts made by Linnéa Smeds

module load bioinfo-tools vcftools/0.1.16 bcftools/1.12


Oland_collared="/crex/proj/sllstore2017033/nobackup/work/alexn/vcf/gt_OC.vcf.gz"
Oland_pied="/crex/proj/sllstore2017033/nobackup/work/alexn/vcf/gt_OP.vcf.gz"


bcftools merge -O z $Oland_collared $Oland_pied > OC.OP.comb.vcf.gz

bcftools query -l $Oland_collared > COL_samples
bcftools query -l $Oland_pied > PIE_samples

####
awk '{if($3 == "Other")print $1 "\t" $2 }' OC.OP.SNP_list.chrom_CpG_annot > OC.OP.non-CpG.SNPs
awk '{if($3 == "CpG_ans")print $1 "\t" $2 }' OC.OP.SNP_list.chrom_CpG_annot > OC.OP.CpG.SNPs

### 
perl translateChromToScaff.pl -in=OC.OP.non-CpG.SNPs -version=20130221 -level=strict -chrCol=1 -posCol=2 -out=OC.OP.non-CpG.SNPs.scaf &
perl translateChromToScaff.pl -in=OC.OP.CpG.SNPs -version=20130221 -level=strict -chrCol=1 -posCol=2 -out=OC.OP.CpG.SNPs.scaf &

wait


pop1="COL_samples"
pop2="PIE_samples"



awk '{print $2 "\t" $3}' OC.OP.CpG.SNPs.scaf > OC.OP.CpG.SNPs.scaf.list
awk '{print $2 "\t" $3}' OC.OP.non-CpG.SNPs.scaf > OC.OP.non-CpG.SNPs.scaf.list


vcftools --gzvcf OC.OP.comb.vcf.gz --positions OC.OP.CpG.SNPs.scaf.list --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size 100 --out CpG
vcftools --gzvcf OC.OP.comb.vcf.gz --positions OC.OP.non-CpG.SNPs.scaf.list --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size 100 --out non-CpG




##
perl translateScaffToChrom.pl -in=CpG.windowed.weir.fst -version=20130221 -level=strict -scafCol=1 -posCol=2,3 -out=CpG.windowed.weir.fst.chrom
perl translateScaffToChrom.pl -in=non-CpG.windowed.weir.fst -version=20130221 -level=strict -scafCol=1 -posCol=2,3 -out=non-CpG.windowed.weir.fst.chrom


awk '{if($2 < $3){print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6} else{print $1 "\t" $3-1 "\t" $2 "\t" $4 "\t" $5 "\t" $6}}'  CpG.windowed.weir.fst.chrom > CpG.windowed.weir.fst.chrom.strand.corr.bed
awk '{if($2 < $3){print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6} else{print $1 "\t" $3-1 "\t" $2 "\t" $4 "\t" $5 "\t" $6}}' non-CpG.windowed.weir.fst.chrom > non-CpG.windowed.weir.fst.chrom.strand.corr.bed


