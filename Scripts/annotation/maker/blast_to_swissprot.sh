#!/bin/bash -l


#Roughly following: http://www.yandell-lab.org/publications/pdf/maker_current_protocols.pdf - Support protocol 3

#wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz #Done 2021-11-01
#gunzip uniprot_sprot.fasta.gz
#mv uniprot_sprot.fasta blast_to_swissprot

ml bioinfo-tools blast/2.11.0+
#ml maker/2.31.10

#mkdir blast_to_swissprot

cd blast_to_swissprot

inputDir="/fAlb15.maker.output"
swissprot="uniprot_sprot.fasta"
makerprot="$inputDir/fAlb15_run4.all.maker.proteins_noSpaceHead.fasta"

makeblastdb -in $swissprot -input_type fasta -dbtype prot

blastp -db $swissprot -query $makerprot -out maker2sp.blastp -evalue .000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads 4

#maker_functional_gff $swissprot maker2sp.blastp.spaceHead $inputDir/fAlb15.all.gff >fAlb15.all.run4_blast2sp.gff
