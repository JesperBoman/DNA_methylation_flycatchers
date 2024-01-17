#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2020-09-24

##### BRAIN ####

library("bsseq")
library("BiocParallel")


load("../bsseq_fit_rda/BS_data_Brain.fit.rda")
load("../bsseq_read_rda/BS_data_Brain.rda")

BS_data_brain <- BS_data
BS_data_brain.fit <- BS_data.fit

colnames(BS_data_brain)[colnames(BS_data_brain) == "HYB04-Brain"] <- "COL06-Brain" # Rename column
colnames(BS_data_brain.fit)[colnames(BS_data_brain.fit) == "HYB04-Brain"] <- "COL06-Brain" # Rename column


BS.cov <- getCoverage(BS_data_brain.fit)

#Pied vs collared
keepLociBrainCvP <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                             rowSums(BS.cov[, 10:14] >= 2) >= 2)
length(keepLociBrainCvP)

CvP_brain.fit.filt <- BS_data_brain.fit[keepLociBrainCvP,]

CvP_brain.fit.filt.tstat <- BSmooth.tstat(CvP_brain.fit.filt, 
                                           group1 = colnames(BS_data_brain)[c(1:5,8)],
                                           group2 = colnames(BS_data_brain)[10:14], 
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)

png(filename = "CvP_brain_t_plot.png")
plot(CvP_brain.fit.filt.tstat)
dev.off()

CvP_brain_dmrs0 <- dmrFinder(CvP_brain.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvP_brain_dmrs <- subset(CvP_brain_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvP_brain.fit.filt, file="CvP_brain.fit.filt.rda")
save(CvP_brain_dmrs, file="CvP_brain_dmrs.rda")
write.table(CvP_brain_dmrs, file="CvP_brain_dmrs.txt", sep="\t", quote =F)


#Collared vs hybrid

keepLociBrainCvH <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                             rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociBrainCvH)

CvH_brain.fit.filt <- BS_data_brain.fit[keepLociBrainCvH,]


CvH_brain.fit.filt.tstat <- BSmooth.tstat(CvH_brain.fit.filt, 
                                           group1 = colnames(BS_data_brain)[c(1:5,8)],
                                           group2 = colnames(BS_data_brain)[c(6:7, 9)], 
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)

png(filename = "CvH_brain_t_plot.png")
plot(CvH_brain.fit.filt.tstat)
dev.off()


CvH_brain_dmrs0 <- dmrFinder(CvH_brain.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvH_brain_dmrs <- subset(CvH_brain_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvH_brain.fit.filt, file="CvH_brain.fit.filt.rda")
save(CvH_brain_dmrs, file="CvH_brain_dmrs.rda")
write.table(CvH_brain_dmrs, file="CvH_brain_dmrs.txt", sep="\t", quote =F)





#Pied vs hybrid
keepLociBrainPvH <- which(rowSums(BS.cov[, 10:14] >= 2) >= 2 &
                             rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociBrainPvH)

PvH_brain.fit.filt <- BS_data_brain.fit[keepLociBrainPvH,]


PvH_brain.fit.filt.tstat <- BSmooth.tstat(PvH_brain.fit.filt, 
                                           group1 = colnames(BS_data_brain)[10:14],
                                           group2 = colnames(BS_data_brain)[c(6:7, 9)], 
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)

png(filename = "PvH_brain_t_plot.png")
plot(PvH_brain.fit.filt.tstat)
dev.off()

PvH_brain_dmrs0 <- dmrFinder(PvH_brain.fit.filt.tstat, qcutoff = c(0.01, 0.99))
PvH_brain_dmrs <- subset(PvH_brain_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(PvH_brain.fit.filt, file="PvH__brain.fit.filt.rda")
save(PvH_brain_dmrs, file="PvH_brain_dmrs.rda")
write.table(PvH_brain_dmrs, file="PvH_brain_dmrs.txt", sep="\t", quote =F)

