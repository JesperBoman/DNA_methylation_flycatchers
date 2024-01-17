#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_HDMR.R input_data "tissue"
#### Jesper Boman - 2022-02-24

library(reshape2)
library(plyr)

args = commandArgs(trailingOnly=T)

args

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")

colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Individual", "Tissue", "Segment_type", "Segment_number", "Coordinate", 
                    "Number_of_CGIs", "Summed_length_CGI", "Weighted_mean_obs_over_exp_CpG_for_CGI", "CGI", "region_ID", "Expression_level", "InpatExpr", "CvP_DMR", "log2FoldChange.COLvPIE", 
                    "padj.COLvPIE", "Meth_mean_COL", "Meth_mean_PIE", "Meth_mean_HYB","CpG_variants", "CpG_Fst", "non_CpG_variants", "non_CpG_Fst", "Exin_annotation", "HvC_DMR", "HvP_DMR")


data$seg_ID <- paste(data$Chromosome, data$Start, data$End, sep="_")
data$Individual <- gsub("HYB04", "COL06", data$Individual)
data$Species <- gsub("\\d", "",  data$Individual)
data$Exin_num <- gsub(".*_", "", data$Exin_annotation)
data$Exin_annotation <- gsub("_.*", "", data$Exin_annotation)


data<-data[!is.na(data$InpatExpr),]
data$DE_status_HvC <- ifelse(data$InpatExpr == "Pied-dominant" | data$InpatExpr == "Additive" | data$InpatExpr == "Overdominant" | data$InpatExpr == "Underdominant", "DE", "non-DE")
data$DE_status_HvP <- ifelse(data$InpatExpr == "Collared-dominant" | data$InpatExpr == "Additive" | data$InpatExpr == "Overdominant" | data$InpatExpr == "Underdominant","DE", "non-DE")


data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status_HvC"), function(x) c(mean(x$HvC_DMR), sd(x[x$Individual == unique(x$Individual)[1],]$HvC_DMR, na.rm=T)/sqrt(length(x[x$Individual == unique(x$Individual)[1],]$HvC_DMR))))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status_HvC"), function (x) c(mean(x$HvC_DMR), sd(x[x$Individual == unique(x$Individual)[1],]$HvC_DMR, na.rm=T)/sqrt(length(x[x$Individual == unique(x$Individual)[1],]$HvC_DMR))))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status_HvC", "DMR_freq", "sem")
all.res$Segment_number <- as.integer(all.res$Segment_number)

save(all.res, file = paste(args[2],  "HvC_DMR_DE_stat.rda", sep="_" ))


res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status_HvP"), function(x) c(mean(x$HvP_DMR), sd(x[x$Individual == unique(x$Individual)[1],]$HvP_DMR, na.rm=T)/sqrt(length(x[x$Individual == unique(x$Individual)[1],]$HvP_DMR))))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status_HvP"), function (x) c(mean(x$HvP_DMR), sd(x[x$Individual == unique(x$Individual)[1],]$HvP_DMR, na.rm=T)/sqrt(length(x[x$Individual == unique(x$Individual)[1],]$HvP_DMR))))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status_HvP", "DMR_freq", "sem")
all.res$Segment_number <- as.integer(all.res$Segment_number)

save(all.res, file = paste(args[2],  "HvP_DMR_DE_stat.rda", sep="_" ))
