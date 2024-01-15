#!/bin/bash -l

ml bioinfo-tools Nextflow/20.10.0
#Prerequisites: install methylseq pipeline in current directory 
export NXF_HOME="Path to methylseq pipeline directory"


NXF_OPTS='-Xms1g -Xmx4g'


nextflow run nf-core/methylseq -resume "Liver" --reads 'tmp_fastq/*Liver*_R{1,2}.fastq.gz' -profile uppmax --project p2018002 --max_cpus 20 --max_memory 128GB --fasta fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa --out_dir ./results_inclReps --save_reference --comprehensive --cytosine_report --epignome
