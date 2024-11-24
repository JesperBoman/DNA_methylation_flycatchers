if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bsseq")

browseVignettes("bsseq")

library("bsseq")

#Example of plotting DMRs

load("CvP_heart.fit.filt.rda")
load("CvP_heart_dmrs.rda")

pData <- pData(CvP_heart.fit.filt)
pData$col <- c(rep("red",times = 5), rep("darkgoldenrod3",times = 2), "red", "darkgoldenrod3", rep("blue",times = 5) )
pData(CvP_heart.fit.filt) <- pData

png(filename = "CvP_heart_top4DMRs.png")
par(mfrow=c(2,2))

for(i in 1:4) {

  plotRegion(CvP_heart.fit.filt, CvP_heart_dmrs[i,], extend = 10000, addRegions = CvP_heart_dmrs)
}
dev.off()
