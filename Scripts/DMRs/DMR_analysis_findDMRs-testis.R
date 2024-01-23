#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2020-09-24

library("bsseq")
library("BiocParallel")

load("../bsseq_fit_rda/BS_data_Testis.fit.rda")
load("../bsseq_read_rda/BS_data_Testis.rda")

BS_data_testis <- BS_data
BS_data_testis.fit <- BS_data.fit

colnames(BS_data_testis)[colnames(BS_data_testis) == "HYB04-Testis"] <- "COL06-Testis" # Rename column
colnames(BS_data_testis.fit)[colnames(BS_data_testis.fit) == "HYB04-Testis"] <- "COL06-Testis" # Rename column


BS.cov <- getCoverage(BS_data_testis.fit)

#Pied vs collared
keepLociTestisCvP <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                            rowSums(BS.cov[, 10:14] >= 2) >= 2)
length(keepLociTestisCvP)

CvP_testis.fit.filt <- BS_data_testis.fit[keepLociTestisCvP,]

CvP_testis.fit.filt.tstat <- BSmooth.tstat(CvP_testis.fit.filt, 
                                          group1 = colnames(BS_data_testis)[c(1:5,8)],
                                          group2 = colnames(BS_data_testis)[10:14], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "CvP_testis_t_plot.png")
plot(CvP_testis.fit.filt.tstat)
dev.off()

CvP_testis_dmrs0 <- dmrFinder(CvP_testis.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvP_testis_dmrs <- subset(CvP_testis_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(CvP_testis.fit.filt, file="CvP_testis.fit.filt.rda")
save(CvP_testis_dmrs, file="CvP_testis_dmrs.rda")
write.table(CvP_testis_dmrs, file="CvP_testis_dmrs.txt", sep="\t", quote =F)



#Collared vs hybrid

keepLociTestisCvH <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                            rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociTestisCvH)

CvH_testis.fit.filt <- BS_data_testis.fit[keepLociTestisCvH,]


CvH_testis.fit.filt.tstat <- BSmooth.tstat(CvH_testis.fit.filt, 
                                          group1 = colnames(BS_data_testis)[c(1:5,8)],
                                          group2 = colnames(BS_data_testis)[c(6:7, 9)], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)


CvH_testis_dmrs0 <- dmrFinder(CvH_testis.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvH_testis_dmrs <- subset(CvH_testis_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)


png(filename = "CvH_testis_t_plot.png")
plot(CvH_testis.fit.filt.tstat)
dev.off()



save(CvH_testis.fit.filt, file="CvH_testis.fit.filt.rda")
save(CvH_testis_dmrs, file="CvH_testis_dmrs.rda")
write.table(CvH_testis_dmrs, file="CvH_testis_dmrs.txt", , sep="\t", quote =F)



#Pied vs hybrid
keepLociTestisPvH <- which(rowSums(BS.cov[, 10:14] >= 2) >= 2 &
                            rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociTestisPvH)

PvH_testis.fit.filt <- BS_data_testis.fit[keepLociTestisPvH,]


PvH_testis.fit.filt.tstat <- BSmooth.tstat(PvH_testis.fit.filt, 
                                          group1 = colnames(BS_data_testis)[10:14],
                                          group2 = colnames(BS_data_testis)[c(6:7, 9)], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "PvH_testis_t_plot.png")
plot(PvH_testis.fit.filt.tstat)
dev.off()

PvH_testis_dmrs0 <- dmrFinder(PvH_testis.fit.filt.tstat, qcutoff = c(0.01, 0.99))
PvH_testis_dmrs <- subset(PvH_testis_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(PvH_testis.fit.filt, file="PvH__testis.fit.filt.rda")
save(PvH_testis_dmrs, file="PvH_testis_dmrs.rda")
write.table(PvH_testis_dmrs, file="PvH_testis_dmrs.txt", sep="\t", quote =F)
