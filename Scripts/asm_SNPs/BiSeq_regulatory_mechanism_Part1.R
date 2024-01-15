#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

args


#### PART 1 #### 
library("BiSeq")
library("BiocParallel")

options(mc.cores=20)

bed_file <- read.table(args[1])

genomic_ranges <- GRanges(seqnames = bed_file$V1, ranges =IRanges(start = bed_file$V2, end = bed_file$V3, names = paste(bed_file$V1, bed_file$V2, bed_file$V3, sep="_")), acc_no=paste(bed_file$V1, bed_file$V2, bed_file$V3, sep="_"))


#for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
for (tissue in args[2:length(args)]  ){
  dir <- paste("/crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/methylseq_pipe/", tissue, "/bismark_methylation_calls/methylation_coverage",sep="")
  setwd(dir)
  files <- system("ls", intern=T)
  
  new_names <-  gsub("HYB04", "COL06", files)
  new_names <- gsub("-.*", "", new_names)
  species <- gsub("\\d", "",  new_names)
  colData <- DataFrame(group = species, row.names = new_names)
  
  bis_data<-readBismark(files=files , colData=colData)
  bis_data.red <- subsetByOverlaps(bis_data, genomic_ranges)
  bis_data.red <- bis_data.red[,c(1:5, 8, 10:14)]
  
  COLsamp <- if(tissue == "Kidney"){sample(c(1:6), size=2)} else{ sample(c(1:6), size=3)}
  paste("COL sample", colnames(bis_data.red)[COLsamp], "is used for", tissue)
  PIEsamp <- if(tissue == "Kidney"){ sample(c(7:11), size=2)} else{ sample(c(7:11), size=3)}
  paste("PIE samples", colnames(bis_data.red)[PIEsamp], "is used for", tissue)
  
  bis_data.red <- bis_data.red[, c(COLsamp, PIEsamp)]
  
  ov <- findOverlaps(bis_data.red, genomic_ranges)
  
  rowRanges(bis_data.red)$cluster.id[queryHits(ov)] <- genomic_ranges$acc_no[subjectHits(ov)] 
  
  bis_data.red.par <- bis_data.red
  
  
  ind.cov <- totalReads(bis_data.red) > 0
  quant <- quantile(totalReads(bis_data.red)[ind.cov], 0.45)
  bis_data.red.lim <- limitCov(bis_data.red, maxCov = quant)
  
  
  ## To shorten the run time set mc.cores, if possible!
  predictedMeth <- predictMeth(object = bis_data.red.lim, mc.cores = 20)
  predictedMeth_COLvPIE.par <- predictedMeth
  
  ## To shorten the run time set mc.cores, if possible! 
  #Takes ca 15 minutes with 20 cores on 1000 promoters
  betaResults_COLvPIE.par <- betaRegression(formula = ~group,
                                            link = "probit", 
                                            mc.cores = 20,
                                            object = predictedMeth_COLvPIE.par,
                                            type = "BR")
  
  
  predictedMeth_COLvPIE.par.null <- predictedMeth
  colData(predictedMeth_COLvPIE.par.null)$group.null <- if(tissue == "Kidney"){c(1,1,2,2)} else{c(1,1,1,2,2,2)}
  
  
  
  #Takes ca 20 minutes with 18 cores on 1000 promoters
  betaResultsNull <- betaRegression(formula = ~group.null,
                                    link = "probit",
                                    object = predictedMeth_COLvPIE.par.null,
                                    mc.cores = 20,
                                    type="BR")
  
  vario <- makeVariogram(betaResultsNull)
  #plot(vario$variogram$v)
  vario.sm <- smoothVariogram(vario, sill = 0.9)
  #lines(vario.sm$variogram[,c("h", "v.sm")],
  #     col = "red", lwd = 1.5)
  vario.aux <- makeVariogram(betaResults_COLvPIE.par, make.variogram=FALSE)
  vario.sm$pValsList <- vario.aux$pValsList
  locCor_COLvPIE.par <- estLocCor(vario.sm)
  clusters.rej.COLvPIE.par <- testClusters(locCor_COLvPIE.par, FDR.cluster = 0.9)
  
  
  
  
  
  
  
  
  ##### COLvPIE - in F1 hybrids #####
  
  dir <- "/crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/asm_SNPs/Bismark_output_CpG_only"
  setwd(dir)
  files <- system("ls", intern=T)
  files <- files[grep(tissue, files)]
  files <- files[grep("bismark.cov.gz", files)]
  
  
  new_names <- paste(gsub("-.*", "", files), gsub(".*cated\\.(.*)\\.bismark.*", "\\1", files), sep=".")
  new_names <- gsub("genome1", "COL", new_names)
  new_names <- gsub("genome2", "PIE", new_names)
  species <- gsub("HYB.*\\.", "",  new_names)
  
  colData <- DataFrame(group = species, row.names = new_names)
  
  bis_data<-readBismark(files=files , colData=colData)
  
  bis_data.red <- subsetByOverlaps(bis_data, genomic_ranges)
  ov <- findOverlaps(bis_data.red, genomic_ranges)
  
  rowRanges(bis_data.red)$cluster.id[queryHits(ov)] <- genomic_ranges$acc_no[subjectHits(ov)] 
  
  bis_data.red.hyb <- bis_data.red
  
  ind.cov <- totalReads(bis_data.red) > 0
  quant <- quantile(totalReads(bis_data.red)[ind.cov], 0.9)
  bis_data.red.lim <- limitCov(bis_data.red, maxCov = quant)
  
  
  ## To shorten the run time set mc.cores, if possible!
  predictedMeth <- predictMeth(object = bis_data.red.lim, mc.cores = 20)
  predictedMeth_COLvPIE.hyb <- predictedMeth
  
  ## To shorten the run time set mc.cores, if possible! 
  #Takes ca 15 minutes with 20 cores on 1000 promoters
  betaResults_COLvPIE.hyb <- betaRegression(formula = ~group,
                                            link = "probit", 
                                            mc.cores = 20,
                                            object = predictedMeth_COLvPIE.hyb,
                                            type = "BR")
  
  
  
  predictedMeth_COLvPIE.hyb.null <- predictedMeth
  
  colData(predictedMeth_COLvPIE.hyb.null)$group.null <- if(tissue == "Kidney"){c(1,2,1,2)} else{ c(1,2,1,2,1,2) }
  
  #Takes ca 20 minutes with 18 cores on 1000 promoters
  betaResultsNull <- betaRegression(formula = ~group.null,
                                    link = "probit",
                                    object = predictedMeth_COLvPIE.hyb.null,
                                    mc.cores = 20,
                                    type="BR")
  
  vario <- makeVariogram(betaResultsNull)
  #plot(vario$variogram$v)
  vario.sm <- smoothVariogram(vario, sill = 0.9)
  #lines(vario.sm$variogram[,c("h", "v.sm")],
  #     col = "red", lwd = 1.5)
  vario.aux <- makeVariogram(betaResults_COLvPIE.hyb, make.variogram=FALSE)
  vario.sm$pValsList <- vario.aux$pValsList
  locCor_COLvPIE.hyb <- estLocCor(vario.sm)
  clusters.rej.COLvPIE.hyb <- testClusters(locCor_COLvPIE.hyb, FDR.cluster = 0.9)
  

  save(bis_data.red.par, bis_data.red.hyb, predictedMeth_COLvPIE.hyb, predictedMeth_COLvPIE.par, file=paste(tissue, "BIGfiles_BiSeq_regulatory_mechanism_data.rda", sep="_"))
  
  
  
  
  
  ##### COL(par) v COL(hyb) ####
  
  pred.COL.par <- predictedMeth_COLvPIE.par[ ,grep("COL", colnames(predictedMeth_COLvPIE.par))]
  colnames(pred.COL.par) <- paste("PAR", colnames(pred.COL.par), sep="-")
  colData(pred.COL.par)$group <- rep("PAR", 3)
  
  
  pred.COL.hyb <- predictedMeth_COLvPIE.hyb[ ,grep("COL", colnames(predictedMeth_COLvPIE.hyb))]
  colnames(pred.COL.hyb) <- paste("HYB", colnames(pred.COL.hyb), sep="-")
  colData(pred.COL.hyb)$group <- rep("HYB", 3)
  
  pred.COL.comb <- combine(pred.COL.par, pred.COL.hyb)
  
  
  betaResults_COL.comb <- betaRegression(formula = ~group,
                                         link = "probit", 
                                         mc.cores = 20,
                                         object = pred.COL.comb,
                                         type = "BR")
  
  
  pred.COL.comb.null <- pred.COL.comb
  colData(pred.COL.comb.null)$group.null <- if(tissue == "Kidney"){c(1,1,2,2)} else{c(1,1,1,2,2,2)}
  
  
  
  #Takes ca 20 minutes with 18 cores on 1000 promoters
  betaResultsNull <- betaRegression(formula = ~group.null,
                                    link = "probit",
                                    object = pred.COL.comb.null,
                                    mc.cores = 20,
                                    type="BR")
  
  vario <- makeVariogram(betaResultsNull)
  #plot(vario$variogram$v)
  vario.sm <- smoothVariogram(vario, sill = 0.9)
  #lines(vario.sm$variogram[,c("h", "v.sm")],
  #    col = "red", lwd = 1.5)
  vario.aux <- makeVariogram(betaResults_COL.comb, make.variogram=FALSE)
  vario.sm$pValsList <- vario.aux$pValsList
  locCor_COL.comb <- estLocCor(vario.sm)
  try(clusters.rej.COL.comb <- testClusters(locCor_COL.comb, FDR.cluster = 0.9))
  
  
  
  #### PIE(par) v PIE(hyb) #####
  
  
  pred.PIE.par <- predictedMeth_COLvPIE.par[ ,grep("PIE", colnames(predictedMeth_COLvPIE.par))]
  colnames(pred.PIE.par) <- paste("PAR", colnames(pred.PIE.par), sep="-")
  colData(pred.PIE.par)$group <- rep("PAR", 3)
  
  
  pred.PIE.hyb <- predictedMeth_COLvPIE.hyb[ ,grep("PIE", colnames(predictedMeth_COLvPIE.hyb))]
  colnames(pred.PIE.hyb) <- paste("HYB", colnames(pred.PIE.hyb), sep="-")
  colData(pred.PIE.hyb)$group <- rep("HYB", 3)
  
  pred.PIE.comb <- combine(pred.PIE.par, pred.PIE.hyb)
  
  betaResults_PIE.comb <- betaRegression(formula = ~group,
                                         link = "probit", 
                                         mc.cores = 20,
                                         object = pred.PIE.comb,
                                         type = "BR")
  
  
  pred.PIE.comb.null <- pred.PIE.comb
  
  colData(pred.PIE.comb.null)$group.null <- if(tissue == "Kidney"){c(1,1,2,2)} else{c(1,1,1,2,2,2)}
  
  
  
  #Takes ca 20 minutes with 18 cores on 1000 promoters
  betaResultsNull <- betaRegression(formula = ~group.null,
                                    link = "probit",
                                    object = pred.PIE.comb.null,
                                    mc.cores = 20,
                                    type="BR")
  
  vario <- makeVariogram(betaResultsNull)
  #plot(vario$variogram$v)
  vario.sm <- smoothVariogram(vario, sill = 0.9)
  #lines(vario.sm$variogram[,c("h", "v.sm")],
  #     PIE = "red", lwd = 1.5)
  vario.aux <- makeVariogram(betaResults_PIE.comb, make.variogram=FALSE)
  vario.sm$pValsList <- vario.aux$pValsList
  locCor_PIE.comb <- estLocCor(vario.sm)
  try(clusters.rej.PIE.comb <- testClusters(locCor_PIE.comb, FDR.cluster = 0.9))
  
  
  
  
  
  
  
  setwd("/PATH/to/asm_SNPs")
  
  save(locCor_COLvPIE.par, locCor_COLvPIE.hyb, clusters.rej.COLvPIE.par, clusters.rej.COLvPIE.hyb,
       locCor_COL.comb, locCor_PIE.comb,   file=paste(tissue, "BiSeq_regulatory_mechanism_data.rda", sep="_"))
  
  
}
