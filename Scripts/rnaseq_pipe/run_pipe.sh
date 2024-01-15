#!/bin/bash -l


ml bioinfo-tools Nextflow/21.02.0-edge
export NXF_HOME="Path to methylseq pipeline directory"

NXF_OPTS='-Xms1g -Xmx4g'



#gffdir="/crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/annotation/maker/run4/fAlb15.maker.output/fAlb15.chrom_run4_featuresOnly.gff"


nextflow run nf-core/rnaseq -name "Testis" --input testis_rna_sample_sheet.csv -profile uppmax --project "" --max_cpus 20 --max_memory 128GB --fasta ../methylseq_pipe/fAlb15.chrom_Olandfix_hardmask_polyCorG100bpthres_softmask.fa --gtf fAlb15.chrom_run4_featuresOnly.gtf --outdir ./Testis --skip_preseq
