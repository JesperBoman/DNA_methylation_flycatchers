#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_plots_comb.R input_data "tissue"
#### Jesper Boman - 2021-11-05

#### Fst profile ####
library(plyr)
library(reshape2)
library(ggplot2)


args = commandArgs(trailingOnly=T)

args

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")


colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Individual", "Tissue", "Segment_type", "Segment_number", "Coordinate", 
                    "Number_of_CGIs", "Summed_length_CGI", "Weighted_mean_obs_over_exp_CpG_for_CGI", "CGI", "region_ID", "Expression_level", "InpatExpr", "CvP_DMR", "log2FoldChange.COLvPIE", 
                    "padj.COLvPIE", "Meth_mean_COL", "Meth_mean_PIE", "Meth_mean_HYB", "CpG_variants", "CpG_Fst", "non_CpG_variants", "non_CpG_Fst")


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


#### DE_status ####

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "DE_status"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "DE_status"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "DE_status", "CpG_Fst", "non_CpG_Fst")
all.res$Segment_number <- as.integer(all.res$Segment_number)

ggplot(data=all.res, aes(x=Segment_number, y=CpG_Fst, col=DE_status))+
  theme_classic()+
  geom_line() +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  #geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))
ggsave(filename = paste(args[2], "DE_status_CpG_Fst.png", sep="_"), plot = last_plot())

ggplot(data=all.res, aes(x=Segment_number, y=non_CpG_Fst, col=DE_status))+
  theme_classic()+
  geom_line() +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  #geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))
ggsave(filename = paste(args[2], "DE_status_Fst.png", sep="_"), plot = last_plot())
save(all.res, file= paste(args[2], "DE_status_Fst.rda", sep="_"))





#### CGI_prom ####

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type",  "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "CGI_prom"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "CGI_prom", "CpG_Fst", "non_CpG_Fst")
all.res$Segment_number <- as.integer(all.res$Segment_number)

ggplot(data=all.res, aes(x=Segment_number, y=CpG_Fst, col=CGI_prom))+
  theme_classic()+
  geom_line() +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  #geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))
ggsave(filename = paste(args[2], "CGIprom_CpG_Fst.png", sep="_"), plot = last_plot())

ggplot(data=all.res, aes(x=Segment_number, y=non_CpG_Fst, col=CGI_prom))+
  theme_classic()+
  geom_line() +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  #geom_hline(yintercept=0, lty=2)+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))
ggsave(filename = paste(args[2], "CGIprom_status_Fst.png", sep="_"), plot = last_plot())
save(all.res, file= paste(args[2], "CGIprom_status_Fst.rda", sep="_"))





#### Cor: Fst v abs(log2fold gene expression) ####

res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "CGI_prom"), 
	function(x) c(cor(x[x$Individual == unique(x$Individual)[1],]$CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$log2FoldChange.COLvPIE), method="spearman", use = "complete.obs"), cor(x[x$Individual == unique(x$Individual)[1],]$non_CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$log2FoldChange.COLvPIE), method="spearman", use = "complete.obs") ))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "CGI_prom"), 
	function(x) c(cor(x[x$Individual == unique(x$Individual)[1],]$CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$log2FoldChange.COLvPIE), method="spearman", use = "complete.obs"), cor(x[x$Individual == unique(x$Individual)[1],]$non_CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$log2FoldChange.COLvPIE), method="spearman", use = "complete.obs") ))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "CGI_prom", "cor_CpG_Fst", "cor_non_CpG_Fst")
all.res$Segment_number <- as.integer(all.res$Segment_number)

ggplot(data=all.res, aes(x=Segment_number, y=cor_CpG_Fst, lty=CGI_prom))+
  theme_classic()+
  geom_line(colour="red", alpha=0.15) +
  geom_hline(yintercept=0, lty=2)+
  geom_smooth(data=all.res, aes(x=Segment_number, y=cor_CpG_Fst), col="red", method = "loess", span = 0.15, method.args = list(degree=1))+
  geom_line( aes(y=cor_non_CpG_Fst, lty=CGI_prom), colour="blue", alpha=0.15)+
  geom_smooth(data=all.res, aes(x=Segment_number, y=cor_non_CpG_Fst),  method = "loess", span = 0.15, method.args = list(degree=1))+
  ylab("cor_Fst_gene_expression")+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

ggsave(filename = paste(args[2], "CGIprom_cor_Fst_GE.png", sep="_"), plot = last_plot())




#### Cor: Fst v abs(log2fold meth) ####


#data$Meth.CvP.log2 <- ifelse(data$Meth_mean_COL/data$Meth_mean_PIE == 0, 0, log((data$Meth_mean_COL+1e9)/(data$Meth_mean_PIE+1e9), base = 2))
data$Meth.CvP.log2 <- abs(data$Meth_mean_COL-data$Meth_mean_PIE) 


res_up_and_down <- ddply(data[data$Segment_type != "Genic",], c("Segment_number", "Segment_type", "CGI_prom"), 
	function(x) c(cor(x[x$Individual == unique(x$Individual)[1],]$CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$Meth.CvP.log2), method="spearman", use = "complete.obs"), cor(x[x$Individual == unique(x$Individual)[1],]$non_CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$Meth.CvP.log2), method="spearman", use = "complete.obs") ))

subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData, c("Segment_number", "Segment_type", "CGI_prom"), 
	function(x) c(cor(x[x$Individual == unique(x$Individual)[1],]$CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$Meth.CvP.log2), method="spearman", use = "complete.obs"), cor(x[x$Individual == unique(x$Individual)[1],]$non_CpG_Fst, abs(x[x$Individual == unique(x$Individual)[1],]$Meth.CvP.log2), method="spearman", use = "complete.obs") ))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "CGI_prom", "cor_CpG_Fst", "cor_non_CpG_Fst")
all.res$Segment_number <- as.integer(all.res$Segment_number)

ggplot(data=all.res, aes(x=Segment_number, y=cor_CpG_Fst, lty=CGI_prom))+
  theme_classic()+
  geom_line(colour="red", alpha=0.15) +
  geom_hline(yintercept=0, lty=2)+
  geom_smooth(data=all.res, aes(x=Segment_number, y=cor_CpG_Fst), col="red", method = "loess", span = 0.15, method.args = list(degree=1))+
  geom_line( aes(y=cor_non_CpG_Fst, lty=CGI_prom), colour="blue", alpha=0.15)+
  geom_smooth(data=all.res, aes(x=Segment_number, y=cor_non_CpG_Fst),  method = "loess", span = 0.15, method.args = list(degree=1))+
  ylab("cor_Fst_absolute_methylation_divergence")+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))

ggsave(filename = paste(args[2], "CGIprom_cor_Fst_Meth.png", sep="_"), plot = last_plot())
save(all.res, file= paste(args[2], "CGIprom_cor_Fst_Meth.rda", sep="_"))
