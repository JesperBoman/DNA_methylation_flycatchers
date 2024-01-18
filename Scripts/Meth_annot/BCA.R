#!/usr/bin/env Rscript

library(ade4)
library(lmodel2)
library(reshape2)

args = commandArgs(trailingOnly=T)

args
meth_annot_data <- read.table(file=args[1], na.strings=c(""," ", "NA"), sep="\t")
colnames(meth_annot_data) <- c("Chromosome", "Start", "End", "region_ID", "Total_CpG", "Post_filter_CpG", "Methylated_reads_per_promoter", "Unmethylated_reads_per_promoter",  "Sum_of_proportions_of_methylated_over_total_reads_per_CpG", "Mean_of_ratios_per_dinuc",  "Individual", "Tissue")


meth_annot_data$Individual <- gsub("HYB04", "COL06", meth_annot_data$Individual)
meth_annot_data$Species <- gsub("\\d", "",  meth_annot_data$Individual)


meth_annot_data_wide <- dcast(meth_annot_data, Tissue+Individual+Species~region_ID, value.var = "Mean_of_ratios_per_dinuc")
meth_annot_data_wide_noNA <- meth_annot_data_wide[ , colSums(is.na(meth_annot_data_wide)) == 0]



firstdatacol=4

#All tissues

p<-dudi.pca(meth_annot_data_wide_noNA[meth_annot_data_wide_noNA$Species != "HYB", c(firstdatacol:ncol(meth_annot_data_wide_noNA)) ], scale=F,scann=F,n=10)
dataset<-meth_annot_data_wide_noNA[meth_annot_data_wide_noNA$Species != "HYB",]

#between species
bp=bca(p,fac=factor(dataset$Species),nf=1,scannf=FALSE)
randbp <- randtest(bca(p,fac=factor(dataset$Species),scannf=FALSE), 999)
plot(randbp, main = "Monte-Carlo test")

print("Between species")
c(bp$ratio,randbp$pvalue)

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
print("Between species, after removing between tissues effects")
wt$ratio
c(bpwt$ratio,randbpwt$pvalue)




bca_tissPCA <- data.frame()

for (tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
    
    meth_annot_tissue <- meth_annot_data_wide_noNA[meth_annot_data_wide$Tissue==tissue  & meth_annot_data_wide$Species != "HYB" ,]

    p <- dudi.pca(meth_annot_tissue[,c(firstdatacol:ncol(meth_annot_tissue))], scale=F,scann=F,n=10)
    tissue_PC_axes <- cbind(as.data.frame(p$li), Species=meth_annot_tissue$Species, Individual=meth_annot_tissue$Individual, Tissue=meth_annot_tissue$Tissue)
    
    #between species
    bp=bca(p,fac=factor(meth_annot_tissue$Species),nf=1,scannf=FALSE)
    randbp <- randtest(bca(p,fac=factor(meth_annot_tissue$Species),scannf=FALSE), 999)
    plot(randbp, main = "Monte-Carlo test")
    c(bp$ratio,randbp$pvalue)
    
    bca_tissPCA <- rbind(bca_tissPCA, cbind(ratio=bp$ratio, pval=randbp$pvalue, tissue))
}  

print(bca_tissPCA)

