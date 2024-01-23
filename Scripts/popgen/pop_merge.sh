#Pop merge
#!/bin/bash -l

module load bioinfo-tools bcftools/1.12

Oland_collared="gt_OC.vcf.gz"
Oland_pied="gt_OP.vcf.gz"


bcftools merge -O z $Oland_collared $Oland_pied > OC.OP.comb.vcf.gz

bcftools query -l $Oland_collared > COL_samples
bcftools query -l $Oland_pied > PIE_samples
