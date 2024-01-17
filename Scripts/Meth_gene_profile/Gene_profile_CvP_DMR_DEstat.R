#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_plots_comb.R input_data "tissue"
#### Jesper Boman - 2021-10-25

#### DMR profile ####
library(plyr)
library(reshape2)
library(ggplot2)


args = commandArgs(trailingOnly=T)

args

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")

colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Individual", "Tissue", "Segment_type", "Segment_number", "Coordinate", 
                    "Number_of_CGIs", "Summed_length_CGI", "Weighted_mean_obs_over_exp_CpG_for_CGI", "CGI", "region_ID", "Expression_level", "InpatExpr", "CvP_DMR", "log2FoldChange.COLvPIE", "padj.COLvPIE", "Meth_mean_COL", "Meth_mean_PIE", "Meth_mean_HYB")

data$seg_ID <- paste(data$Chromosome, data$Start, data$End, sep="_")
data$Individual <- gsub("HYB04", "COL06", data$Individual)
data$Species <- gsub("\\d", "",  data$Individual)


data$DE_status <- ifelse(data$padj.COLvPIE < 0.1 & !is.na(data$padj.COLvPIE), "DE", "non-DE")


data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status"), function(x) c(mean(x$CvP_DMR), sd(x[x$Individual == unique(x$Individual)[1],]$CvP_DMR, na.rm=T)/sqrt(length(x[x$Individual == unique(x$Individual)[1],]$CvP_DMR))))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status"), function (x) c(mean(x$CvP_DMR), sd(x[x$Individual == unique(x$Individual)[1],]$CvP_DMR, na.rm=T)/sqrt(length(x[x$Individual == unique(x$Individual)[1],]$CvP_DMR))))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status", "DMR_freq", "sem")
all.res$Segment_number <- as.integer(all.res$Segment_number)

ggplot(data=all.res, aes(x=Segment_number, y=DMR_freq, col=DE_status, fill=DE_status, ymin=DMR_freq-sem, ymax=DMR_freq+sem)) +
  theme_classic()+
  geom_smooth(method = 'loess')+
  geom_ribbon(alpha=0.3, colour=NA)+
  geom_line() +
  ylab("DMR frequency")+
  ylim(0,0.1)+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

ggsave(filename = paste(args[2], "CvP_DMR_DEstat.png", sep="_"), plot = last_plot())
save(all.res, file = paste(args[2], "CvP_DMR_DEstat.rda", sep="_"))
table(data[match(unique(data$Gene_name), data$Gene_name), ]$DE_status )
