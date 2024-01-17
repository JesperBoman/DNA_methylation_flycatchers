#!/usr/bin/env Rscript
#### A script to obtain R data that can be used to plot gene profile features
#### Usage: Rscript Gene_profile_acrossTissCor.R input_data "group"
#### Jesper Boman - 2022-02-21

library(reshape2)
library(plyr)

args = commandArgs(trailingOnly=T)

args

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")

colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Individual", "Tissue", "Segment_type", "Segment_number", "Coordinate", 
                    "Number_of_CGIs", "Summed_length_CGI", "Weighted_mean_obs_over_exp_CpG_for_CGI", "CGI", "region_ID", "Expression_level", "InpatExpr", "CvP_DMR", "log2FoldChange.COLvPIE", 
                    "padj.COLvPIE", "Meth_mean_COL", "Meth_mean_PIE", "CpG_variants", "CpG_Fst", "non_CpG_variants", "non_CpG_Fst", "Exin_annotation")


data$seg_ID <- paste(data$Chromosome, data$Start, data$End, sep="_")
data$Individual <- gsub("HYB04", "COL06", data$Individual)
data$Species <- gsub("\\d", "",  data$Individual)
data$Exin_num <- gsub(".*_", "", data$Exin_annotation)
data$Exin_annotation <- gsub("_.*", "", data$Exin_annotation)



data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "Gene_name"), function(x) cor(x$Mean_per_dinuc, x$Expression_level, method="spearman", use = "complete.obs" ))


subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData,  c("Segment_number", "Segment_type", "Gene_name"), function(x) cor(x$Mean_per_dinuc, x$Expression_level, method="spearman", use = "complete.obs" ))


all.res <- rbind(res_up_and_down, res_body)

colnames(all.res) <- c("Segment_number", "Segment_type", "Gene_name", "Cor")
all.res$Segment_number <- as.integer(all.res$Segment_number)

promData <- data[data$Segment_type == "Upstream" & data$Segment_number > 30 ,]
CGIgenes <- as.data.frame(table(promData$CGI, promData$Gene_name))
CGIgenes <- dcast(CGIgenes, Var2~Var1)
colnames(CGIgenes) <- c("Gene", "N", "Y")
CGIgenes$CGI <- ifelse(CGIgenes$Y > 0, "Y", "N")

all.res$CGI_prom <- NA
for (i in 1:length(CGIgenes$Gene)){
  matches <- CGIgenes$Gene[i] == all.res$Gene_name
  m_list <- (1:length(matches))[matches]
  all.res$CGI_prom[m_list] <- CGIgenes$CGI[i]
  
  if( i %% 1000 == 0){print(i)}
}


save(all.res, file = paste(args[2],  "acrossTissCor_noBrain.rda", sep="_" ))
