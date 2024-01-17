#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_HybFst.R input_data "tissue"
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
data$DE_status_HvP <- ifelse(data$InpatExpr == "Collared-dominant" | data$InpatExpr == "Additive" | data$InpatExpr == "Overdominant" | data$InpatExpr == "Underdominant", "DE", "non-DE")


#Filtering for callable windows
data$CpG_Fst <- ifelse(data$CpG_Fst == 0 & data$CpG_variants == 0 , NA, data$CpG_Fst)
data$non_CpG_Fst <- ifelse(data$non_CpG_Fst == 0 & data$non_CpG_variants == 0, NA, data$non_CpG_Fst)
#

#Setting negative Fst as 0
data$CpG_Fst <- ifelse(data$CpG_Fst < 0, 0, data$CpG_Fst)
data$non_CpG_Fst <- ifelse(data$non_CpG_Fst < 0, 0, data$non_CpG_Fst)
#


data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149


#### DE_status ####

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status_HvC", "DE_status_HvP"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status_HvC", "DE_status_HvP"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status_HvC", "DE_status_HvP", "CpG_Fst", "non_CpG_Fst")
all.res$Segment_number <- as.integer(all.res$Segment_number)

save(all.res, file= paste(args[2], "hybridDE_status_Fst.rda", sep="_"))
