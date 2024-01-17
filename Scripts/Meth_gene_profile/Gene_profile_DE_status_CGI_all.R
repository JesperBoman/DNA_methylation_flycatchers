#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_plots_comb.R input_data "tissue"
#### Jesper Boman - 2022-03-24

#### DE stat profiles per promoter type ####
library(plyr)
library(reshape2)
library(ggplot2)


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

data$DE_status <- ifelse(data$padj.COLvPIE < 0.1 & !is.na(data$padj.COLvPIE), "DE", "non-DE")


#Filtering for callable windows
data$CpG_Fst <- ifelse(data$CpG_Fst == 0 & data$CpG_variants == 0 , NA, data$CpG_Fst)
data$non_CpG_Fst <- ifelse(data$non_CpG_Fst == 0 & data$non_CpG_variants == 0, NA, data$non_CpG_Fst)
#

#Setting negative Fst as 0
data$CpG_Fst <- ifelse(data$CpG_Fst < 0, 0, data$CpG_Fst)
data$non_CpG_Fst <- ifelse(data$non_CpG_Fst < 0, 0, data$non_CpG_Fst)
#

data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149

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


#### CvP DE_status ####

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$CvP_DMR, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$CvP_DMR, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status", "CGI_prom", "CpG_Fst", "non_CpG_Fst", "CvP_DMR_freq")
all.res$Segment_number <- as.integer(all.res$Segment_number)


save(all.res, file= paste(args[2], "CvP_DE_status_CGI_all.rda", sep="_"))

print(args[2])
table(data[match(unique(data$Gene_name), data$Gene_name), ]$DE_status, data[match(unique(data$Gene_name), data$Gene_name), ]$CGI_prom )


#### HvC DE_status ####

data$DE_status_HvC <- ifelse(data$InpatExpr == "Pied-dominant" | data$InpatExpr == "Additive" | data$InpatExpr == "Overdominant" | data$InpatExpr == "Underdominant", "DE", "non-DE")

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status_HvC", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$HvC_DMR, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status_HvC", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$HvC_DMR, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status_HvC", "CGI_prom", "CpG_Fst", "non_CpG_Fst", "HvC_DMR_freq")
all.res$Segment_number <- as.integer(all.res$Segment_number)


save(all.res, file= paste(args[2], "HvC_DE_status_CGI_all.rda", sep="_"))

print(args[2])
table(data[match(unique(data$Gene_name), data$Gene_name), ]$DE_status_HvC, data[match(unique(data$Gene_name), data$Gene_name), ]$CGI_prom )



#### HvP DE_status ####

data$DE_status_HvP <- ifelse(data$InpatExpr == "Collared-dominant" | data$InpatExpr == "Additive" | data$InpatExpr == "Overdominant" | data$InpatExpr == "Underdominant", "DE", "non-DE")

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status_HvP", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$HvP_DMR, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status_HvP", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$HvP_DMR, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status_HvP", "CGI_prom", "CpG_Fst", "non_CpG_Fst", "HvP_DMR_freq")
all.res$Segment_number <- as.integer(all.res$Segment_number)


save(all.res, file= paste(args[2], "HvP_DE_status_CGI_all.rda", sep="_"))

print(args[2])
table(data[match(unique(data$Gene_name), data$Gene_name), ]$DE_status_HvP, data[match(unique(data$Gene_name), data$Gene_name), ]$CGI_prom )





#### Transgressive DE_status ####

data$DE_status_Transgressive <- ifelse(data$InpatExpr == "Overdominant" | data$InpatExpr == "Underdominant", "DE", "non-DE")



res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status_Transgressive", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$HvC_DMR, na.rm=T), mean(x$HvP_DMR, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status_Transgressive", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T), mean(x$HvC_DMR, na.rm=T), mean(x$HvP_DMR, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status_Transgressive", "CGI_prom", "CpG_Fst", "non_CpG_Fst", "HvC_DMR_freq", "HvP_DMR_freq")
all.res$Segment_number <- as.integer(all.res$Segment_number)


save(all.res, file= paste(args[2], "Transgressive_DE_status_CGI_all.rda", sep="_"))

print(args[2])
table(data[match(unique(data$Gene_name), data$Gene_name), ]$DE_status_Transgressive, data[match(unique(data$Gene_name), data$Gene_name), ]$CGI_prom )
