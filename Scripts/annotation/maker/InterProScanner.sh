#!/bin/bash -l


ml bioinfo-tools InterProScan/5.30-69.0

inputDir="/Scripts/annotation/maker/run4/fAlb15.maker.output"

sed 's/\s/_/g' $inputDir/fAlb15_run4.all.maker.proteins.fasta > $inputDir/fAlb15_run4.all.maker.proteins_noSpaceHead.fasta

protFasta="fAlb15_run4.all.maker.proteins_noSpaceHead.fasta"


/sw/apps/bioinfo/InterProScan/5.30-69.0/rackham/interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i $inputDir/$protFasta -o ${protFasta}.iprscan

prot_with_domain_n=$(awk '{a[$1]} END{for (prot in a){sum++} print sum }' ${protFasta}.iprscan)
prot_tot_n=$(grep ">" $inputDir/$protFasta | wc -l)

echo $prot_tot_n $prot_with_domain_n

echo "scale=3 ;  $prot_with_domain_n / $prot_tot_n" | bc


#ml maker/2.31.10

#ipr_update_gff $inputDir/fAlb15.all.gff ${protFasta}.iprscan > fAlb15.all.run4.putative_function.domain_added.gff
