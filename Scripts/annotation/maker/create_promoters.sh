#!/bin/bash -l

awk '{if($3 == "gene" && $7 == "+"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+).*Name/, a); plus_gene[$1,$4]= $1 "\t" $4 "\t" $5 "\t" $7 "\t" a[1]} 
  if($3 == "gene" && $7 == "-"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+).*Name/, a); minus_gene[$1,$5]= $1 "\t" $4 "\t" $5 "\t" $7 "\t" a[1]} 
  if($3 == "five_prime_UTR" && $7 == "+"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+)-mRNA-[0-9]+:/, a); plus_utr[$1,$4]= $1 "\t" $4 "\t" $5 "\t" $7 "\t" a[1]} 
  if($3 == "five_prime_UTR" && $7 == "-"){match($9, /ID=(.*gene-[0-9]+\.[0-9]+)-mRNA-[0-9]+:/, a); minus_utr[$1,$5]=  $1 "\t" $4 "\t" $5 "\t" $7 "\t" a[1]}} 
  END{for (beg in plus_gene){if(beg in plus_utr){print plus_gene[beg]}} for (beg in minus_gene){if(beg in minus_utr){print minus_gene[beg]}}}' fAlb15_run4_featuresOnly.gff > fAlb15_run4_genes_with_conc_five_prime_UTR_IDs.bed

perl translateScaffToChrom.pl -in=fAlb15_run4_featuresOnly.gff -version=20130221 -level=strict -scafCol=1 -posCol=4,5 -format=gtf -out=fAlb15.chrom_run4_featuresOnly.gff 

awk 'FNR==NR{a[$5]} FNR!=NR{if($3 == "gene"){ match($9, /ID=(.*gene-[0-9]+\.[0-9]+).*Name/, b); if(b[1] in a){if($7 =="+"){  print $1 "\t" $4-2001 "\t" $4} else if($7 == "-") {  print $1 "\t" $5-1 "\t" $5+2000}}}}' fAlb15_run4_genes_with_conc_five_prime_UTR_IDs.bed fAlb15.chrom_run4_featuresOnly.gff > promoters_2kb_run4_1.1.bed

