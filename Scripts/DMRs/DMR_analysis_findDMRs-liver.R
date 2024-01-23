#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2020-09-24

##### LIVER ####

library("bsseq")
library("BiocParallel")


load("../bsseq_fit_rda/BS_data_Liver.fit.rda")
load("../bsseq_read_rda/BS_data_Liver.rda")

BS_data_liver <- BS_data
BS_data_liver.fit <- BS_data.fit

colnames(BS_data_liver)[colnames(BS_data_liver) == "HYB04-Liver"] <- "COL06-Liver" # Rename column
colnames(BS_data_liver.fit)[colnames(BS_data_liver.fit) == "HYB04-Liver"] <- "COL06-Liver" # Rename column


BS.cov <- getCoverage(BS_data_liver.fit)

#Pied vs collared
keepLociLiverCvP <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                             rowSums(BS.cov[, 10:14] >= 2) >= 2)
length(keepLociLiverCvP)

CvP_liver.fit.filt <- BS_data_liver.fit[keepLociLiverCvP,]

CvP_liver.fit.filt.tstat <- BSmooth.tstat(CvP_liver.fit.filt, 
                                           group1 = colnames(BS_data_liver)[c(1:5,8)],
                                           group2 = colnames(BS_data_liver)[10:14], 
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)

png(filename = "CvP_liver_t_plot.png")
plot(CvP_liver.fit.filt.tstat)
dev.off()

CvP_liver_dmrs0 <- dmrFinder(CvP_liver.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvP_liver_dmrs <- subset(CvP_liver_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvP_liver.fit.filt, file="CvP_liver.fit.filt.rda")
save(CvP_liver_dmrs, file="CvP_liver_dmrs.rda")
write.table(CvP_liver_dmrs, file="CvP_liver_dmrs.txt", sep="\t", quote =F)


#Collared vs hybrid

keepLociLiverCvH <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                             rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociLiverCvH)

CvH_liver.fit.filt <- BS_data_liver.fit[keepLociLiverCvH,]


CvH_liver.fit.filt.tstat <- BSmooth.tstat(CvH_liver.fit.filt, 
                                           group1 = colnames(BS_data_liver)[c(1:5,8)],
                                           group2 = colnames(BS_data_liver)[c(6:7, 9)], 
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)

png(filename = "CvH_liver_t_plot.png")
plot(CvH_liver.fit.filt.tstat)
dev.off()


CvH_liver_dmrs0 <- dmrFinder(CvH_liver.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvH_liver_dmrs <- subset(CvH_liver_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvH_liver.fit.filt, file="CvH_liver.fit.filt.rda")
save(CvH_liver_dmrs, file="CvH_liver_dmrs.rda")
write.table(CvH_liver_dmrs, file="CvH_liver_dmrs.txt", sep="\t", quote =F)





#Pied vs hybrid
keepLociLiverPvH <- which(rowSums(BS.cov[, 10:14] >= 2) >= 2 &
                             rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociLiverPvH)

PvH_liver.fit.filt <- BS_data_liver.fit[keepLociLiverPvH,]


PvH_liver.fit.filt.tstat <- BSmooth.tstat(PvH_liver.fit.filt, 
                                           group1 = colnames(BS_data_liver)[10:14],
                                           group2 = colnames(BS_data_liver)[c(6:7, 9)], 
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)

png(filename = "PvH_liver_t_plot.png")
plot(PvH_liver.fit.filt.tstat)
dev.off()

PvH_liver_dmrs0 <- dmrFinder(PvH_liver.fit.filt.tstat, qcutoff = c(0.01, 0.99))
PvH_liver_dmrs <- subset(PvH_liver_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(PvH_liver.fit.filt, file="PvH__liver.fit.filt.rda")
save(PvH_liver_dmrs, file="PvH_liver_dmrs.rda")
write.table(PvH_liver_dmrs, file="PvH_liver_dmrs.txt", sep="\t", quote =F)
