#!/bin/bash -l

#The BUSCO installation that came with the Oyster-River Protocol installation did not work for me so I ran separately

module load bioinfo-tools BUSCO/3.0.2b

source $BUSCO_SETUP

run_BUSCO.py -i /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/annotation/ORP/assemblies/COL.ORP.fasta -o COL.ORP -l $BUSCO_LINEAGE_SETS/aves_odb9 -m transcriptome -c 20
run_BUSCO.py -i /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/annotation/ORP/reports/transrate_COL/COL.ORP/good.COL.ORP.fasta -o transrate.good.COL.ORP -l $BUSCO_LINEAGE_SETS/aves_odb9 -m transcriptome -c 20
