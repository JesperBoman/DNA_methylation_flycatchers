#!/bin/bash -l

source $HOME/Oyster_River_Protocol/software/anaconda/install/bin/activate orp

#cd $HOME/Oyster_River_Protocol/sampledata

#$HOME/Oyster_River_Protocol/oyster.mk \
#STRAND=RF \
#TPM_FILT=0.2 \
#MEM=15 \
#CPU=8 \
#READ1=test.1.fq.gz \
#READ2=test.2.fq.gz \
#RUNOUT=test


$HOME/Oyster_River_Protocol/oyster.mk \
TPM_FILT=1 \
STRAND=RF \
MEM=512 \
CPU=16 \
READ1=all_COL_R1.fq.gz \
READ2=all_COL_R2.fq.gz \
RUNOUT=COL
