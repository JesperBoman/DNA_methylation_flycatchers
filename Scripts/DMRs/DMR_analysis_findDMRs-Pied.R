#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs-Pied.R 
#### Jesper Boman - 2020-09-30

##### DMRs BETWEEN TISSUES WITHIN SPECIES #####

#### Pied ####

library("bsseq")
library("BiocParallel")

args = commandArgs(trailingOnly = TRUE)

#load("../../BS_data_P_KidneyvTestis.fit.rda")
#load("../../BS_data_P_KidneyvTestis.rda")

load(args[1])
load(args[2])


BS.cov <- getCoverage(BS_data.fit)

keepLoci<- which(rowSums(BS.cov[, 1:5] >= 2) >= 2 &
                            rowSums(BS.cov[, 6:10] >= 2) >= 2)
length(keepLoci)

comp.fit.filt <- BS_data.fit[keepLoci,]


comp.fit.filt.tstat <- BSmooth.tstat(comp.fit.filt, 
                                          group1 = colnames(BS_data)[1:5],
                                          group2 = colnames(BS_data)[6:10], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "t_plot.png")
plot(comp.fit.filt.tstat )
dev.off()


dmrs0 <- dmrFinder(comp.fit.filt.tstat , qcutoff = c(0.01, 0.99))
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(comp.fit.filt, file="comp.fit.filt.rda")
save(dmrs, file="dmrs.rda")
write.table(dmrs, file="dmrs.txt", sep="\t", quote =F, row.names=F)
dmrs$start <- dmrs$start - 1
write.table(dmrs, file="dmrs.bed", sep="\t", quote =F, row.names=F)
