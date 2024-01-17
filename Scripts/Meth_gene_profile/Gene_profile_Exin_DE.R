#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_exin.R input_data "tissue" "group"
#### Jesper Boman - 2022-02-08

library(ggplot2)
library(reshape2)
library(plyr)



# Exon/Intron annotation ####

args = commandArgs(trailingOnly=T)

args

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")
colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Individual", "Tissue", "Segment_type", "Segment_number", "Coordinate", 
                    "Number_of_CGIs", "Summed_length_CGI", "Weighted_mean_obs_over_exp_CpG_for_CGI", "CGI", "region_ID", "Expression_level", "InpatExpr", "CvP_DMR", "log2FoldChange.COLvPIE", 
                    "padj.COLvPIE", "Meth_mean_COL", "Meth_mean_PIE", "Meth_mean_HYB", "CpG_variants", "CpG_Fst", "non_CpG_variants", "non_CpG_Fst", "Exin_annotation")

promData <- data[data$Segment_type == "Upstream" & data$Segment_number > 30 ,]
CGIgenes <- as.data.frame(table(promData$CGI, promData$Gene_name))
CGIgenes <- dcast(CGIgenes, Var2~Var1)
colnames(CGIgenes) <- c("Gene", "N", "Y")
CGIgenes$CGI <- ifelse(CGIgenes$Y > 0, "Y", "N")


data$CGI_prom <- NA
for (i in 1:length(CGIgenes$Gene)){
  matches <- CGIgenes$Gene[i] == data$Gene_name
  m_list <- (1:length(matches))[matches]
  data$CGI_prom[m_list] <- CGIgenes$CGI[i]
  
  if( i %% 1000 == 0){print(i)}
}

data <- data[!is.na(data$CGI_prom),]

data$seg_ID <- paste(data$Chromosome, data$Start, data$End, sep="_")
data$Individual <- gsub("HYB04", "COL06", data$Individual)
data$Species <- gsub("\\d", "",  data$Individual)
data$Exin_num <- gsub(".*_", "", data$Exin_annotation)
data$Exin_annotation <- gsub("_.*", "", data$Exin_annotation)

#Filtering for callable windows
data$CpG_Fst <- ifelse(data$CpG_Fst == 0 & data$CpG_variants == 0 , NA, data$CpG_Fst)
data$non_CpG_Fst <- ifelse(data$non_CpG_Fst == 0 & data$non_CpG_variants == 0, NA, data$non_CpG_Fst)
#

#Setting negative Fst as 0
data$CpG_Fst <- ifelse(data$CpG_Fst < 0, 0, data$CpG_Fst)
data$non_CpG_Fst <- ifelse(data$non_CpG_Fst < 0, 0, data$non_CpG_Fst)
#

data$Meth.CvP <- abs(data$Meth_mean_COL-data$Meth_mean_PIE) 

data$DE_status <- ifelse(data$padj.COLvPIE < 0.1 & !is.na(data$padj.COLvPIE), "DE", "non-DE")

data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149


res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "Individual", "CGI_prom", "DE_status"), function(x) c(mean(x$CvP_DMR), mean(x$Meth.CvP, na.rm=T), mean(x$Mean_per_dinuc, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$CpG_Fst, na.rm=T), cor(x$Mean_per_dinuc, log10(x$Expression_level+1), method="spearman", use = "complete.obs" )))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
subData <- subData[!is.na(subData$Exin_annotation),]
res_body <- ddply(subData, c("Segment_number", "Segment_type","Individual", "Exin_annotation", "CGI_prom", "DE_status"), function(x) c(mean(x$CvP_DMR), mean(x$Meth.CvP, na.rm=T), mean(x$Mean_per_dinuc, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$CpG_Fst, na.rm=T), cor(x$Mean_per_dinuc, log10(x$Expression_level+1), method="spearman", use = "complete.obs" )))


res_up_and_down$Exin_annotation <- NA
res_up_and_down[res_up_and_down$Segment_type == "Upstream", ]$Exin_annotation <- "upstream"
res_up_and_down1 <- res_up_and_down[res_up_and_down$Segment_type == "Upstream", ]

res_up_and_down[res_up_and_down$Segment_type == "Downstream", ]$Exin_annotation <- "downstream"
res_up_and_down2 <- res_up_and_down[res_up_and_down$Segment_type == "Downstream", ]

all.res <- rbind(res_up_and_down1, res_up_and_down2, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type","Individual", "CGI_prom", "DE_status", "DMR_frequency", "Meth.CvP","Meth", "non_CpG_Fst", "CpG_Fst", "Cor", "Exin_annotation")
all.res$Segment_number <- as.integer(all.res$Segment_number)

save(all.res, file = paste(args[2], args[3],  "Exin_CvP_DE.rda", sep="_" ))
