#!/bin/bash -l


module load bioinfo-tools maker/2.31.10

cd fAlb15.maker.output

gff3_merge -n -d fAlb15_master_datastore_index.log -o fAlb15_run4.gff &

fasta_merge -d fAlb15_master_datastore_index.log -o fAlb15_run4 &

wait
