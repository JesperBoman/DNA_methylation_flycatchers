#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs-Hybrid.R 
#### Jesper Boman - 2021-11-01

##### DMRs BETWEEN TISSUES WITHIN SPECIES #####

#### Hybrid ####

library("bsseq")
library("BiocParallel")


#load("../../BS_data_H_KidneyvBrain.fit.rda")
#load("../../BS_data_H_KidneyvBrain.rda")

args = commandArgs(trailingOnly = TRUE)

load(args[1])
load(args[2])

tissue1=args[3]
tissue2=args[4]

if(tissue1 == "Kidney" | tissue2 == "Kidney"){

if(tissue1 == "Kidney"){
colpatgroup1=c(1,2)
colpatgroup2=c(3,4,5)
}
if(tissue2 == "Kidney"){
colpatgroup1=c(1,2,3)
colpatgroup2=c(4,5)
}


BS.cov <- getCoverage(BS_data.fit)

keepLoci<- which(rowSums(BS.cov[, colpatgroup1] >= 2) >= 2 &
                            rowSums(BS.cov[, colpatgroup2] >= 2) >= 2)
length(keepLoci)

comp.fit.filt <- BS_data.fit[keepLoci,]


comp.fit.filt.tstat <- BSmooth.tstat(comp.fit.filt, 
                                          group1 = colnames(BS_data)[colpatgroup1],
                                          group2 = colnames(BS_data)[colpatgroup2], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

} else{
BS.cov <- getCoverage(BS_data.fit)

keepLoci<- which(rowSums(BS.cov[, 1:3] >= 2) >= 2 &
                            rowSums(BS.cov[, 4:6] >= 2) >= 2)
length(keepLoci)

comp.fit.filt <- BS_data.fit[keepLoci,]


comp.fit.filt.tstat <- BSmooth.tstat(comp.fit.filt, 
                                          group1 = colnames(BS_data)[1:3],
                                          group2 = colnames(BS_data)[4:6], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

}


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
