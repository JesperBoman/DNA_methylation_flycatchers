#!/usr/bin/env Rscript

library("reshape2")
library("plyr")

args = commandArgs(trailingOnly=T)


TS <- read.table("tmp_TSproms_PEM")
nonTS <- read.table("tmp_nonTSproms_PEM")
ref_test <- read.table("ref_test")

colnames(ref_test) <- c("Species", "Tissue", "Gene", "region_ID", "PEM", "Chromosome", "Chr_type", "Promoter_type", "Set")
colnames(TS) <- c("Species", "Tissue", "Gene", "region_ID", "PEM", "Chromosome", "Chr_type", "Promoter_type", "Set")
colnames(nonTS) <- c("Species", "Tissue", "Gene", "region_ID", "PEM", "Chromosome", "Chr_type", "Promoter_type", "Set")
df_com <- rbind(ref_test, TS)  


rank_table <- ddply(df_com, c("Gene"), function (x) rank(x$PEM))
colnames(rank_table) <-  c("Gene", unique(df_com$Tissue))
rank_table_long <- melt(rank_table)

tissue=args[1]
TS_ranks_frq <- as.data.frame(table(rank_table_long$variable, rank_table_long$value))
TS_ranks_frq <- TS_ranks_frq[TS_ranks_frq$Var1 == tissue,]

TS_ranks_frq$Exp  <- sum(TS_ranks_frq$Freq)/5
Chi2pval <- pchisq(sum(((TS_ranks_frq$Exp -TS_ranks_frq$Freq)^2 / TS_ranks_frq$Exp)),  df=4, lower.tail=F) 
fisherpval <-fisher.test(cbind(TS_ranks_frq$Freq, TS_ranks_frq$Exp))$p.value

print(paste(unique(TS$Species), unique(TS$Tissue), unique(TS$Promoter_type), sum(TS_ranks_frq$Freq), Chi2pval, fisherpval))


#

test <- try(t.test(TS[TS$PEM > 0,]$PEM, nonTS[nonTS$PEM >0,]$PEM))

pval <- try(wilcox.test(TS[TS$PEM > 0,]$PEM, nonTS[nonTS$PEM >0,]$PEM)$p.value)



print(paste("Overexpressed", unique(TS$Species), unique(TS$Tissue), unique(TS$Promoter_type),  test$conf.int[1],  test$conf.int[2],  test$p.value, pval[1]))
#
#

test <-try(t.test(TS[TS$PEM < 0,]$PEM, nonTS[nonTS$PEM <0,]$PEM))

pval <- try(wilcox.test(TS[TS$PEM < 0,]$PEM, nonTS[nonTS$PEM < 0,]$PEM)$p.value)
print(paste("Underexpressed", unique(TS$Species), unique(TS$Tissue), unique(TS$Promoter_type),  test$conf.int[1],  test$conf.int[2],  test$p.value, pval[1]))
#
