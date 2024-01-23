#!/bin/bash -l


module load bioinfo-tools vcftools/0.1.16 

hyperythra="HYP.gatk.allsites.vcf.bgz"
parva="parv.gatk.allsites.vcf.bgz"


vcftools --gzvcf $hyperythra --counts --out "hyperythra" &
vcftools --gzvcf $parva --counts --out "parva" &

wait

awk '{if($3 == 1 && $4 == 2){split($5, a, ":"); print $1 "\t" $2 "\t" a[1] }}' hyperythra.frq.count > hyperythra.AA &
awk '{if($3 == 1 && $4 == 2){split($5, a, ":"); print $1 "\t" $2 "\t" a[1] }}' parva.frq.count > parva.AA &

wait 
awk 'NR==FNR{a[$1 "\t" $2]=$3} NR!=FNR{if(a[$1 "\t" $2] == $3) print $1 "\t" $2 "\t" $3}' hyperythra.AA parva.AA > shared.ancestral.states

perl translateScaffToChrom.pl -in=shared.ancestral.states -version=20130221 -level=strict -scafCol=1 -posCol=2 -out=shared.ancestral.states.chrom
