#!/usr/bin/env Rscript
#### A script to plot gene profile features
#### Usage: Rscript Gene_profile_plots_comb.R input_data "tissue" "group"
#### Jesper Boman - 2021-06-04

library(ggplot2)
library(reshape2)
library(plyr)


args = commandArgs(trailingOnly=T)

args

data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")

colnames(data) <- c("Chromosome", "Start", "End", "Gene_name", "Segment_num_text", "Strand", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_feature", "Unmethylated_reads_per_feature",  
                    "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_per_dinuc", "Individual", "Tissue", "Segment_type", "Segment_number", "Coordinate", 
                    "Number_of_CGIs", "Summed_length_CGI", "Weighted_mean_obs_over_exp_CpG_for_CGI", "CGI", "region_ID", "Expression_level")

data$Individual <- gsub("HYB04", "COL06", data$Individual)

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



data$Expression_cat <- NA
for (ind in unique(data$Individual)){
  
data[data$Expression_level <= quantile(data[data$Individual == ind,]$Expression_level, 1/5, na.rm=T) & data$Individual == ind,]$Expression_cat <- "L"
data[data$Expression_level > quantile(data[data$Individual == ind,]$Expression_level, 1/5, na.rm=T) & data$Expression_level <= quantile(data[data$Individual == ind,]$Expression_level, 4/5, na.rm=T)  & data$Individual == ind,]$Expression_cat <- "M"
data[data$Expression_level > quantile(data[data$Individual == ind,]$Expression_level, 4/5, na.rm=T) & data$Expression_level <= quantile(data[data$Individual == ind,]$Expression_level, 1, na.rm=T) & data$Individual == ind,]$Expression_cat <- "H"

}





data[data$Segment_type == "Downstream",]$Segment_number <- data[data$Segment_type == "Downstream",]$Segment_number + 149
res_up_and_down <- ddply(data[data$Segment_type != "Genic" & !is.na(data$CGI_prom),], c("Segment_number", "Segment_type", "Individual", "CGI_prom", "Expression_cat"), function(x) c(mean(x$Mean_per_dinuc, na.rm=T), sd(x$Mean_per_dinuc, na.rm=T)/sqrt(length(x$Mean_per_dinuc))))


subData <- data[-grep("stream", data$Segment_type),]
br <- quantile(subData$Coordinate, seq(0, 1, length=100), na.rm=TRUE)
mid <- br[-length(br)] + diff(br)/2
subData$Segment_number <- cut(subData$Coordinate, breaks=br, labels=51:149, include.lowest=T)
res_body <- ddply(subData[!is.na(subData$CGI_prom),], c("Segment_number", "Segment_type", "Individual", "CGI_prom", "Expression_cat"), function(x) c(mean(x$Mean_per_dinuc, na.rm=T), sd(x$Mean_per_dinuc, na.rm=T)/sqrt(length(x$Mean_per_dinuc))))

all.res <- rbind(res_up_and_down, res_body)
colnames(all.res) <- c("Segment_number", "Segment_type", "Individual", "CGI_prom", "Expression_cat", "Methylation_level", "sem")
all.res$Segment_number <- as.integer(all.res$Segment_number)
all.res$Expression_cat <- factor(all.res$Expression_cat, levels=c("L","M","H"))

for (i in c("Y", "N")){
  ggplot(data=all.res[all.res$CGI_prom == i,], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Individual), col=Expression_cat, fill=Expression_cat,  ymin=Methylation_level-sem, ymax=Methylation_level+sem)) +
    theme_classic()+
    geom_line() +
    geom_ribbon(alpha=0.3, colour = NA) +
    ylab("Methylation level")+
    labs(colour="Expression category")+
    ylim(0,1)+
    guides(fill=FALSE)+
    scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
    theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=12),  legend.title=element_text(size=14))
  ggsave(filename = paste(args[2], args[3], "CGI_prom", i,  "Meth_x_RNA_LMH.png", sep="_"), plot = last_plot())
}

save(all.res, file = paste(args[2], args[3],  "Meth_x_RNA_LMH.rda", sep="_" ))
