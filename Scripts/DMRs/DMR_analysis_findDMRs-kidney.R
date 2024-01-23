#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2020-09-24


##### KIDNEY ####

library("bsseq")
library("BiocParallel")


load("../bsseq_fit_rda/BS_data_Kidney.fit.rda")
load("../bsseq_read_rda/BS_data_Kidney.rda")

BS_data_kidney <- BS_data
BS_data_kidney.fit <- BS_data.fit

colnames(BS_data_kidney)[colnames(BS_data_kidney) == "HYB04-Kidney"] <- "COL06-Kidney" # Rename column
colnames(BS_data_kidney.fit)[colnames(BS_data_kidney.fit) == "HYB04-Kidney"] <- "COL06-Kidney" # Rename column


BS.cov <- getCoverage(BS_data_kidney.fit)

#Pied vs collared
keepLociKidneyCvP <- which(rowSums(BS.cov[, c(1:5,7)] >= 2) >= 2 &
                            rowSums(BS.cov[, 9:13] >= 2) >= 2)
length(keepLociKidneyCvP)

CvP_kidney.fit.filt <- BS_data_kidney.fit[keepLociKidneyCvP,]

CvP_kidney.fit.filt.tstat <- BSmooth.tstat(CvP_kidney.fit.filt, 
                                          group1 = colnames(BS_data_kidney)[c(1:5,7)],
                                          group2 = colnames(BS_data_kidney)[9:13], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "CvP_kidney_t_plot.png")
plot(CvP_kidney.fit.filt.tstat)
dev.off()

CvP_kidney_dmrs0 <- dmrFinder(CvP_kidney.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvP_kidney_dmrs <- subset(CvP_kidney_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvP_kidney.fit.filt, file="CvP_kidney.fit.filt.rda")
save(CvP_kidney_dmrs, file="CvP_kidney_dmrs.rda")
write.table(CvP_kidney_dmrs, file="CvP_kidney_dmrs.txt", sep="\t", quote =F, row.names=F)


#Collared vs hybrid

keepLociKidneyCvH <- which(rowSums(BS.cov[, c(1:5,7)] >= 2) >= 2 &
                            rowSums(BS.cov[, c(6, 8)] >= 2) >= 2)
length(keepLociKidneyCvH)

CvH_kidney.fit.filt <- BS_data_kidney.fit[keepLociKidneyCvH,]


CvH_kidney.fit.filt.tstat <- BSmooth.tstat(CvH_kidney.fit.filt, 
                                          group1 = colnames(BS_data_kidney)[c(1:5,7)],
                                          group2 = colnames(BS_data_kidney)[c(6,8)], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "CvH_kidney_t_plot.png")
plot(CvH_kidney.fit.filt.tstat)
dev.off()


CvH_kidney_dmrs0 <- dmrFinder(CvH_kidney.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvH_kidney_dmrs <- subset(CvH_kidney_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvH_kidney.fit.filt, file="CvH_kidney.fit.filt.rda")
save(CvH_kidney_dmrs, file="CvH_kidney_dmrs.rda")
write.table(CvH_kidney_dmrs, file="CvH_kidney_dmrs.txt", sep="\t", quote =F, row.names =F)





#Pied vs hybrid
keepLociKidneyPvH <- which(rowSums(BS.cov[, 9:13] >= 2) >= 2 &
                            rowSums(BS.cov[, c(6,8)] >= 2) >= 2)
length(keepLociKidneyPvH)

PvH_kidney.fit.filt <- BS_data_kidney.fit[keepLociKidneyPvH,]


PvH_kidney.fit.filt.tstat <- BSmooth.tstat(PvH_kidney.fit.filt, 
                                          group1 = colnames(BS_data_kidney)[9:13],
                                          group2 = colnames(BS_data_kidney)[c(6,8)], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "PvH_kidney_t_plot.png")
plot(PvH_kidney.fit.filt.tstat)
dev.off()

PvH_kidney_dmrs0 <- dmrFinder(PvH_kidney.fit.filt.tstat, qcutoff = c(0.01, 0.99))
PvH_kidney_dmrs <- subset(PvH_kidney_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(PvH_kidney.fit.filt, file="PvH_kidney.fit.filt.rda")
save(PvH_kidney_dmrs, file="PvH_kidney_dmrs.rda")
write.table(PvH_kidney_dmrs, file="PvH_kidney_dmrs.txt", sep="\t", quote =F, row.names =F)
