#!/bin/bash -l


####
awk '{if($3 == "Other")print $1 "\t" $2 }' OC.OP.SNP_list.chrom_CpG_annot > OC.OP.non-CpG.SNPs
awk '{if($3 == "CpG_ans")print $1 "\t" $2 }' OC.OP.SNP_list.chrom_CpG_annot > OC.OP.CpG.SNPs

### 
perl /crex/proj/sllstore2017033/repos/scripts/translateChromToScaff.pl -in=OC.OP.non-CpG.SNPs -version=20130221 -level=strict -scafCol=1 -posCol=2 -out=OC.OP.non-CpG.SNPs.scaf &
perl /crex/proj/sllstore2017033/repos/scripts/translateChromToScaff.pl -in=OC.OP.CpG.SNPs -version=20130221 -level=strict -scafCol=1 -posCol=2 -out=OC.OP.CpG.SNPs.scaf &

wait



module load bioinfo-tools vcftools/0.1.16 

pop1="COL_samples"
pop2="PIE_samples"


vcftools --gzvcf OC.OP.comb.vcf.gz --positions <(awk '{print $2 "\t" $3}' OC.OP.CpG.SNPs.scaf) --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out CpG
vcftools --gzvcf OC.OP.comb.vcf.gz --positions <(awk '{print $2 "\t" $3}' OC.OP.non-CpG.SNPs.scaf) --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out non-CpG

awk '{print $2 "\t" $3}' OC.OP.CpG.SNPs.scaf > OC.OP.CpG.SNPs.scaf.list
awk '{print $2 "\t" $3}' OC.OP.non-CpG.SNPs.scaf > OC.OP.non-CpG.SNPs.scaf.list


vcftools --gzvcf OC.OP.comb.vcf.gz --positions OC.OP.CpG.SNPs.scaf.list --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size 2000 --out CpG_2kb
vcftools --gzvcf OC.OP.comb.vcf.gz --positions OC.OP.non-CpG.SNPs.scaf.list --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size 2000 --out non-CpG_2kb
