#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2020-09-24

library("bsseq")
library("BiocParallel")

#### HEART ####
load("../bsseq_fit_rda/BS_data_Heart.fit.rda")
load("../bsseq_read_rda/BS_data_Heart.rda")


BS_data_heart <- BS_data
BS_data_heart.fit <- BS_data.fit

colnames(BS_data_heart)[colnames(BS_data_heart) == "HYB04-Heart"] <- "COL06-Heart" # Rename column
colnames(BS_data_heart.fit)[colnames(BS_data_heart.fit) == "HYB04-Heart"] <- "COL06-Heart" # Rename column


BS.cov <- getCoverage(BS_data_heart.fit)

#Pied vs collared
keepLociHeartCvP <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                    rowSums(BS.cov[, 10:14] >= 2) >= 2)
length(keepLociHeartCvP)

CvP_heart.fit.filt <- BS_data_heart.fit[keepLociHeartCvP,]

CvP_heart.fit.filt.tstat <- BSmooth.tstat(CvP_heart.fit.filt, 
                                        group1 = colnames(BS_data_heart)[c(1:5,8)],
                                        group2 = colnames(BS_data_heart)[10:14], 
                                        estimate.var = "same",
                                        local.correct = TRUE,
                                        verbose = TRUE)

png(filename = "CvP_heart_t_plot.png")
plot(CvP_heart.fit.filt.tstat)
dev.off()

CvP_heart_dmrs0 <- dmrFinder(CvP_heart.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvP_heart_dmrs <- subset(CvP_heart_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

pData <- pData(CvP_heart.fit.filt)
pData$col <- c(rep("red",times = 5), rep("darkgoldenrod3",times = 2), "red", "darkgoldenrod3", rep("blue",times = 5) )
pData(CvP_heart.fit.filt) <- pData

#png(filename = "CvP_heart_top4DMRs.png")
#par(mfrow=c(2,2))


#for(i in 1:4) {
 # plotRegion(CvP_heart.fit.filt, CvP_heart_dmrs[i,], extend = 10000, addRegions = CvP_heart_dmrs)
#}
#dev.off()

save(CvP_heart.fit.filt, file="CvP_heart.fit.filt.rda")
save(CvP_heart_dmrs, file="CvP_heart_dmrs.rda")

write.table(CvP_heart_dmrs, file="CvP_heart_dmrs.txt", sep="\t", quote =F)

#Collared vs hybrid

keepLociHeartCvH <- which(rowSums(BS.cov[, c(1:5,8)] >= 2) >= 2 &
                            rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociHeartCvH)

CvH_heart.fit.filt <- BS_data_heart.fit[keepLociHeartCvH,]


CvH_heart.fit.filt.tstat <- BSmooth.tstat(CvH_heart.fit.filt, 
                                          group1 = colnames(BS_data_heart)[c(1:5,8)],
                                          group2 = colnames(BS_data_heart)[c(6:7, 9)], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "CvH_heart_t_plot.png")
plot(CvH_heart.fit.filt.tstat)
dev.off()

CvH_heart_dmrs0 <- dmrFinder(CvH_heart.fit.filt.tstat, qcutoff = c(0.01, 0.99))
CvH_heart_dmrs <- subset(CvH_heart_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

pData <- pData(CvH_heart.fit.filt)
pData$col <- c(rep("red",times = 5), rep("darkgoldenrod3",times = 2), "red", "darkgoldenrod3", rep("blue",times = 5) )
pData(CvH_heart.fit.filt) <- pData

#png(filename = "CvH_heart_top4DMRs.png")
#par(mfrow=c(2,2))


#for(i in 1:4) {
#	plotRegion(CvH_heart.fit.filt, CvH_heart_dmrs[i,], extend = 10000, addRegions = CvH_heart_dmrs)
#}
#dev.off()

save(CvH_heart.fit.filt, file="CvH_heart.fit.filt.rda")
save(CvH_heart_dmrs, file="CvH_heart_dmrs.rda")

write.table(CvH_heart_dmrs, file="CvH_heart_dmrs.txt", sep="\t", quote =F)


#Pied vs hybrid
keepLociHeartPvH <- which(rowSums(BS.cov[, 10:14] >= 2) >= 2 &
                            rowSums(BS.cov[, c(6:7, 9)] >= 2) >= 2)
length(keepLociHeartPvH)

PvH_heart.fit.filt <- BS_data_heart.fit[keepLociHeartPvH,]


PvH_heart.fit.filt.tstat <- BSmooth.tstat(PvH_heart.fit.filt, 
                                          group1 = colnames(BS_data_heart)[10:14],
                                          group2 = colnames(BS_data_heart)[c(6:7, 9)], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "PvH_heart_t_plot.png")
plot(PvH_heart.fit.filt.tstat)
dev.off()

PvH_heart_dmrs0 <- dmrFinder(PvH_heart.fit.filt.tstat, qcutoff = c(0.01, 0.99))
PvH_heart_dmrs <- subset(PvH_heart_dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

pData <- pData(PvH_heart.fit.filt)
pData$col <- c(rep("red",times = 5), rep("darkgoldenrod3",times = 2), "red", "darkgoldenrod3", rep("blue",times = 5) )
pData(PvH_heart.fit.filt) <- pData

#png(filename = "PvH_heart_top4DMRs.png")
#par(mfrow=c(2,2))

#for(i in 1:4)  {
#  plotRegion(PvH_heart.fit.filt, PvH_heart_dmrs[i,], extend = 10000, addRegions = PvH_heart_dmrs)
#}
#dev.off()

save(PvH_heart.fit.filt, file="PvH_heart.fit.filt.rda")
save(PvH_heart_dmrs, file="PvH_heart_dmrs.rda")

write.table(PvH_heart_dmrs, file="PvH_heart_dmrs.txt", sep="\t", quote =F)
