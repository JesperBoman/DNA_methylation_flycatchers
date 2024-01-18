#!/usr/bin/env Rscript

library(ggplot2)
library(viridis)
library(ade4)
library(lmodel2)
library(reshape2)
library(ggforce)
library(plyr)

args = commandArgs(trailingOnly=T)

args



#load("prom_meth_annot_PCA_data.rda")

meth_annot_data <- read.table("meth_promoters_2kb_run4_1.1.bed_results_perPROM_nm_SNPsIncl_nP_1.1_CGI_annot")
meth_annot_data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")
#colnames(meth_annot_data) <- c("Chromosome", "Start", "End", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_promoter", "Unmethylated_reads_per_promoter",  "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_of_ratios_per_dinuc",  "Individual", "Tissue")
#colnames(meth_annot_data) <- c("Chromosome", "Start", "End", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_promoter", "Unmethylated_reads_per_promoter",  "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_of_ratios_per_dinuc", "CGI_num_per_region", "CGI_sum_len", "Mean_obs_over_expected_CG",  "Individual", "Tissue")

#Promoter specific
colnames(meth_annot_data) <- c("Chromosome", "Start", "End", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_promoter", "Unmethylated_reads_per_promoter",  "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_of_ratios_per_dinuc", "Individual", "Tissue", "CGI_num_per_region", "CGI_sum_len", "Mean_obs_over_expected_CG")
#

meth_annot_data$Individual <- gsub("HYB04", "COL06", meth_annot_data$Individual)
meth_annot_data$Species <- gsub("\\d", "",  meth_annot_data$Individual)
meth_annot_data$region_ID <- paste(meth_annot_data$Chromosome, meth_annot_data$Start, meth_annot_data$End, sep="_") 

#Promoter specific
meth_annot_data$Promoter_type <- ifelse(meth_annot_data$CGI_num_per_region == 0, "Other", "CGI")

meth_annot_data_wide_CGI <- dcast(meth_annot_data[meth_annot_data$Promoter_type == "CGI",], Tissue+Individual+Species~region_ID, value.var = "Mean_of_ratios_per_dinuc")
meth_annot_data_wide_CGI_noNA <- meth_annot_data_wide_CGI[ , colSums(is.na(meth_annot_data_wide_CGI)) == 0]

meth_annot_data_wide_Other <- dcast(meth_annot_data[meth_annot_data$Promoter_type == "Other",], Tissue+Individual+Species~region_ID, value.var = "Mean_of_ratios_per_dinuc")
meth_annot_data_wide_Other_noNA <- meth_annot_data_wide_Other[ , colSums(is.na(meth_annot_data_wide_Other)) == 0]
#



# Principal component analysis ####

firstdatacol=4

if (args[3] == "skipPCA"){
  print("PCA skipped")
} else{
  
  
  
  #PCA per region
  
  for (promoter_type in c("Other", "CGI")){
    
    if(promoter_type == "CGI"){
      p<-dudi.pca(meth_annot_data_wide_CGI_noNA[, c(firstdatacol:ncol(meth_annot_data_wide_CGI_noNA))], scale=F,scann=F,n=10)
      meth_annot_PC_axes <- cbind(as.data.frame(p$li), Species=meth_annot_data_wide_CGI_noNA$Species, Individual=meth_annot_data_wide_CGI_noNA$Individual, Tissue=meth_annot_data_wide_CGI_noNA$Tissue) 
    }
    
    if(promoter_type == "Other"){
    p<-dudi.pca(meth_annot_data_wide_Other_noNA[, c(firstdatacol:ncol(meth_annot_data_wide_Other_noNA))], scale=F,scann=F,n=10)
    meth_annot_PC_axes <- cbind(as.data.frame(p$li), Species=meth_annot_data_wide_Other_noNA$Species, Individual=meth_annot_data_wide_Other_noNA$Individual, Tissue=meth_annot_data_wide_Other_noNA$Tissue) 
    }
  
  
  for (pc in c(2, 3, 4, 5)){

    eig=round(100*p$eig/sum(p$eig),2) 
    
    n1=paste("Axis", 1,sep="")
    n2=paste("Axis", pc ,sep="")
    
    ggplot(data=meth_annot_PC_axes, aes_string(x=n1, y=n2, col="Tissue", shape="Species"))+geom_point(size=4)+
      theme_classic()+
      scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
      #geom_text(aes(label=Individual), hjust=-0.2, vjust=0)+
      xlab(paste("PC1 (",eig[1],"% explained variance)")) + 
      ylab(paste("PC2 (",eig[pc],"% explained variance)")) + 
      theme(aspect.ratio=1, legend.position = "right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))
    
    ggsave(filename = paste(args[2], "allTissues", promoter_type, "PC1", pc, ".png", sep="_"), plot = last_plot())
  }
    

  }
  
  
  #TISSUE-specific PCA
  

  for (promoter_type in c("Other", "CGI")){
    
  
  for (tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
    
    if(promoter_type == "CGI"){
      meth_annot_tissue <- meth_annot_data_wide_CGI_noNA[meth_annot_data_wide_CGI_noNA$Tissue==tissue ,]
    }
    if(promoter_type == "Other"){
      meth_annot_tissue <- meth_annot_data_wide_Other_noNA[meth_annot_data_wide_Other_noNA$Tissue==tissue ,]
    }
    
    p <- dudi.pca(meth_annot_tissue[,c(firstdatacol:ncol(meth_annot_tissue))], scale=F,scann=F,n=10)
    tissue_PC_axes <- cbind(as.data.frame(p$li), Species=meth_annot_tissue$Species, Individual=meth_annot_tissue$Individual, Tissue=meth_annot_tissue$Tissue)
    

    for (pc in 2:5){
      eig=round(100*p$eig/sum(p$eig),2) 
      
      n1=paste("Axis", 1,sep="")
      n2=paste("Axis", pc ,sep="")
      
      ggplot(data=tissue_PC_axes, aes_string(x=n1, y=n2, col="Species"))+geom_point(size=9)+
        scale_color_manual(name="Species", values = c("red", "darkgoldenrod3", "blue")) +
        theme_classic()+
        #geom_text(aes(label=Individual), hjust=-0.2, vjust=0)+
        xlab(paste0("PC1 (",eig[1],"%)")) + 
        ylab(paste0("PC2 (",eig[pc],"%)")) + 
        theme(aspect.ratio=1, legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))
      
      #ggsave(filename = paste(args[2], tissue, promoter_type, "PC1", pc, ".png", sep="_"), plot = last_plot())
  }
  }
  }
}


#BCA analysis ####

#Based on original script by Carina Mugal

bca_tissPCA <- data.frame()
for (promoter_type in c("Other", "CGI")){
  
  
  for (tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
    
    if(promoter_type == "CGI"){
      meth_annot_tissue <- meth_annot_data_wide_CGI_noNA[meth_annot_data_wide_CGI_noNA$Tissue==tissue  & meth_annot_data_wide_CGI_noNA$Species != "HYB" ,]
    }
    if(promoter_type == "Other"){
      meth_annot_tissue <- meth_annot_data_wide_Other_noNA[meth_annot_data_wide_Other_noNA$Tissue==tissue & meth_annot_data_wide_Other_noNA$Species != "HYB",]
    }
    
    p <- dudi.pca(meth_annot_tissue[,c(firstdatacol:ncol(meth_annot_tissue))], scale=F,scann=F,n=10)
    tissue_PC_axes <- cbind(as.data.frame(p$li), Species=meth_annot_tissue$Species, Individual=meth_annot_tissue$Individual, Tissue=meth_annot_tissue$Tissue)
    
    #between species
    bp=bca(p,fac=factor(meth_annot_tissue$Species),nf=1,scannf=FALSE)
    randbp <- randtest(bca(p,fac=factor(meth_annot_tissue$Species),scannf=FALSE), 999)
    plot(randbp, main = "Monte-Carlo test")
    c(bp$ratio,randbp$pvalue)
    
    bca_tissPCA <- rbind(bca_tissPCA, cbind(ratio=bp$ratio, pval=randbp$pvalue, promoter_type, tissue))
}}    
    

# across all tissues: compute the proportion of variance associated with species identity and tissue (hybrids removed)
# -> The variance associated with species identity is very small and does not significantly differ from random expectations (0.65% of total variance explained by difference between species).
#    It is still small but becomes significant if you compute it after removing the part of the variance explained by difference between tissues (4.66%)

if(promoter_type == "CGI"){
  p<-dudi.pca(meth_annot_data_wide_CGI_noNA[meth_annot_data_wide_CGI_noNA$Species != "HYB", c(firstdatacol:ncol(meth_annot_data_wide_CGI_noNA)) ], scale=F,scann=F,n=10)
  dataset<-meth_annot_data_wide_CGI_noNA[meth_annot_data_wide_CGI_noNA$Species != "HYB",]
}

if(promoter_type == "Other"){
  p<-dudi.pca(meth_annot_data_wide_Other_noNA[meth_annot_data_wide_Other_noNA$Species != "HYB", c(firstdatacol:ncol(meth_annot_data_wide_Other_noNA)) ], scale=F,scann=F,n=10)
  dataset<-meth_annot_data_wide_Other_noNA[meth_annot_data_wide_Other_noNA$Species != "HYB",]
}

set.seed(1000)
#between species
bp=bca(p,fac=factor(dataset$Species),nf=1,scannf=FALSE)
randbp <- randtest(bca(p,fac=factor(dataset$Species),scannf=FALSE), 999)
plot(randbp, main = "Monte-Carlo test")
c(bp$ratio,randbp$pvalue)

loocv(bp)

# between tissues
bt=bca(p,fac=factor(dataset$Tissue),nf=1,scannf=FALSE)
randbt <- randtest(bca(p,fac=factor(dataset$Tissue),scannf=FALSE), 999)

print("Between tissues")
c(bt$ratio, randbt$pvalue)


# between species, after removing between tissues effects
wt=wca(p,fac=factor(dataset$Tissue),scannf=FALSE)
bpwt=bca(wt,fac=factor(dataset$Species),nf=1,scannf=FALSE)
randbpwt <- randtest(bca(wt,fac=factor(dataset$Species),scannf=FALSE), 999)
plot(randbpwt, main = "Monte-Carlo test")
# bp = wt * bpwt
wt$ratio
c(bpwt$ratio,randbpwt$pvalue)

