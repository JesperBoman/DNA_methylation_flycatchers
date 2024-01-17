#!/usr/bin/env Rscript

#I open this script in RStudio and run it from there

library(ggplot2)

setwd(paste(getwd(), "Downloads", sep="/"))
load(file="COL_Brain_Meth_x_RNA_Ind_cor.rda")

COL.Brain.all.res <- all.res
COL.Brain.all.res$Species <- "COL"

load(file="PIE_Brain_Meth_x_RNA_Ind_cor.rda")
PIE.Brain.all.res <- all.res
PIE.Brain.all.res$Species <- "PIE"

load(file="HYB_Brain_Meth_x_RNA_Ind_cor.rda")
HYB.Brain.all.res <- all.res
HYB.Brain.all.res$Species <- "HYB"

Brain.spec.all.res <- rbind(COL.Brain.all.res, PIE.Brain.all.res, HYB.Brain.all.res)

ggplot(data=Brain.spec.all.res[Brain.spec.all.res$CGI_prom == "Y",], aes(x=Segment_number, y=Cor, group = interaction(Species, Individual), col=Species)) +
  theme_classic()+
  scale_color_manual(name="Species", values = c("red", "darkgoldenrod3", "blue")) +
  geom_hline(yintercept=0, lty=2)+
  geom_line() +
  ylab("Spearman's r")+
  labs(fill="Individual")+
  ylim(-0.5,0.5)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=Brain.spec.all.res[Brain.spec.all.res$CGI_prom == "Y",], aes(x=Segment_number, y=Cor, col=Species)) +
  theme_classic()+
  scale_color_manual(name="Species", values = c("red", "darkgoldenrod3", "blue")) +
  geom_hline(yintercept=0, lty=2)+
  geom_line()+
  ylab("Spearman's r")+
  ylim(-0.5,0.5)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))




load(file="COL_Liver_Meth_x_RNA_Ind_cor.rda")
COL.Liver.all.res <- all.res
COL.Liver.all.res$Species <- "COL"

load(file="PIE_Liver_Meth_x_RNA_Ind_cor.rda")
PIE.Liver.all.res <- all.res
PIE.Liver.all.res$Species <- "PIE"

load(file="HYB_Liver_Meth_x_RNA_Ind_cor.rda")
HYB.Liver.all.res <- all.res
HYB.Liver.all.res$Species <- "HYB"

Liver.spec.all.res <- rbind(COL.Liver.all.res, PIE.Liver.all.res, HYB.Liver.all.res)

ggplot(data=Liver.spec.all.res[Liver.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Cor, group = interaction(Species, Individual), col=Species)) +
  theme_classic()+
  scale_color_manual(name="Species", values = c("red", "darkgoldenrod3", "blue")) +
  geom_hline(yintercept=0, lty=2)+
  geom_line() +
  ylab("Spearman's r")+
  labs(fill="Individual")+
  ylim(-0.5,0.5)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


ggplot(data=all.agg.cgi.full, aes(x=Segment_number, y=Prop, fill=Individual)) +
  theme_classic()+
  geom_line(size=1) +
  ylab("CGI proportion")+
  labs(fill="Individual")+
  ylim(0,1)+
  scale_x_continuous(name ="5 Kb Upstream                        Gene body                         5 Kb Downstream", labels=c(-5,0,0.5,0,5))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=20), legend.text=element_text(size=24),  legend.title=element_text(size=26))


load(file="COL_Heart_Meth_x_RNA_Ind_cor.rda")
COL.Heart.all.res <- all.res
COL.Heart.all.res$Species <- "COL"

load(file="PIE_Heart_Meth_x_RNA_Ind_cor.rda")
PIE.Heart.all.res <- all.res
PIE.Heart.all.res$Species <- "PIE"

load(file="HYB_Heart_Meth_x_RNA_Ind_cor.rda")
HYB.Heart.all.res <- all.res
HYB.Heart.all.res$Species <- "HYB"

Heart.spec.all.res <- rbind(COL.Heart.all.res, PIE.Heart.all.res, HYB.Heart.all.res)


load(file="COL_Kidney_Meth_x_RNA_Ind_cor.rda")
COL.Kidney.all.res <- all.res
COL.Kidney.all.res$Species <- "COL"

load(file="PIE_Kidney_Meth_x_RNA_Ind_cor.rda")
PIE.Kidney.all.res <- all.res
PIE.Kidney.all.res$Species <- "PIE"

load(file="HYB_Kidney_Meth_x_RNA_Ind_cor.rda")
HYB.Kidney.all.res <- all.res
HYB.Kidney.all.res$Species <- "HYB"

Kidney.spec.all.res <- rbind(COL.Kidney.all.res, PIE.Kidney.all.res, HYB.Kidney.all.res)


load(file="COL_Testis_Meth_x_RNA_Ind_cor.rda")
COL.Testis.all.res <- all.res
COL.Testis.all.res$Species <- "COL"

load(file="PIE_Testis_Meth_x_RNA_Ind_cor.rda")
PIE.Testis.all.res <- all.res
PIE.Testis.all.res$Species <- "PIE"

load(file="HYB_Testis_Meth_x_RNA_Ind_cor.rda")
HYB.Testis.all.res <- all.res
HYB.Testis.all.res$Species <- "HYB"

Testis.spec.all.res <- rbind(COL.Testis.all.res, PIE.Testis.all.res, HYB.Testis.all.res)



Brain.spec.all.res$Tissue <- "Brain"
Heart.spec.all.res$Tissue <- "Heart"
Kidney.spec.all.res$Tissue <- "Kidney"
Liver.spec.all.res$Tissue <- "Liver"
Testis.spec.all.res$Tissue <- "Testis"

all.tiss.spec.all.res <- rbind(Brain.spec.all.res, Heart.spec.all.res, Kidney.spec.all.res, Liver.spec.all.res, Testis.spec.all.res)

ggplot(data=all.tiss.spec.all.res[all.tiss.spec.all.res$CGI_prom == "Y",], aes(x=Segment_number, y=Cor, group = interaction(Tissue, Species, Individual), col=Tissue, fill=Tissue,lty=Species )) +
  theme_classic()+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) + 
  
  geom_line() +
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_hline(yintercept=0, lty=2)+
  ylab(expression("Spearman's" ~ rho))+
  #geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  ylim(-0.5,0.5)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))




ggplot(data=all.tiss.spec.all.res[all.tiss.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Cor, group = interaction(Tissue, Species, Individual), col=Tissue, fill=Tissue,lty=Species )) +
  theme_classic()+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) + 

  geom_line() +
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_hline(yintercept=0, lty=2)+
  ylab(expression("Spearman's" ~ rho))+
  #geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  ylim(-0.5,0.5)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="right", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))








load(file="COL_Brain_Meth_x_RNA_LMH.rda")
COL.Brain.all.res <- all.res
COL.Brain.all.res$Species <- "COL"

load(file="PIE_Brain_Meth_x_RNA_LMH.rda")
PIE.Brain.all.res <- all.res
PIE.Brain.all.res$Species <- "PIE"

load(file="HYB_Brain_Meth_x_RNA_LMH.rda")
HYB.Brain.all.res <- all.res
HYB.Brain.all.res$Species <- "HYB"

Brain.spec.all.res <- rbind(COL.Brain.all.res, PIE.Brain.all.res, HYB.Brain.all.res)



 ggplot(data=Brain.spec.all.res[Brain.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Individual), col=Expression_cat, fill=Expression_cat,  ymin=Methylation_level-sem, ymax=Methylation_level+sem, lty = Species )) +
  #scale_color_manual(name="Expression category", values = c("#004D40", "#FFC107", "#D81B60")) +
 # scale_fill_manual(name="Species", values = c("#004D40", "#FFC107", "#D81B60")) +
  theme_bw()+
  scale_color_viridis_d(option="viridis")+
  scale_fill_viridis_d(option="viridis")+
  geom_line() +
  geom_ribbon(alpha=0.3, color=NA) +
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("Methylation level")+
  labs(colour="Expression category")+
  ylim(0,1)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size = 6)))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


load(file="COL_Liver_Meth_x_RNA_LMH.rda")
COL.Liver.all.res <- all.res
COL.Liver.all.res$Species <- "COL"

load(file="PIE_Liver_Meth_x_RNA_LMH.rda")
PIE.Liver.all.res <- all.res
PIE.Liver.all.res$Species <- "PIE"

load(file="HYB_Liver_Meth_x_RNA_LMH.rda")
HYB.Liver.all.res <- all.res
HYB.Liver.all.res$Species <- "HYB"

Liver.spec.all.res <- rbind(COL.Liver.all.res, PIE.Liver.all.res, HYB.Liver.all.res)

LY <- ggplot(data=Liver.spec.all.res[Liver.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Individual), col=Expression_cat, lty=Species, fill=Expression_cat,  ymin=Methylation_level-sem, ymax=Methylation_level+sem )) +
  theme_bw()+
  scale_color_viridis_d(option="viridis")+
  scale_fill_viridis_d(option="viridis")+
  geom_line() +
  geom_ribbon(alpha=0.3, colour = NA) +
  ylab("Methylation level")+
  labs(colour="Expression category")+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylim(0,1)+
  guides(fill=FALSE)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


load(file="COL_Testis_Meth_x_RNA_LMH.rda")
COL.Testis.all.res <- all.res
COL.Testis.all.res$Species <- "COL"

load(file="PIE_Testis_Meth_x_RNA_LMH.rda")
PIE.Testis.all.res <- all.res
PIE.Testis.all.res$Species <- "PIE"

load(file="HYB_Testis_Meth_x_RNA_LMH.rda")
HYB.Testis.all.res <- all.res
HYB.Testis.all.res$Species <- "HYB"

Testis.spec.all.res <- rbind(COL.Testis.all.res, PIE.Testis.all.res, HYB.Testis.all.res)

TY <- ggplot(data=Testis.spec.all.res[Testis.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Individual), col=Expression_cat, lty=Species, fill=Expression_cat,  ymin=Methylation_level-sem, ymax=Methylation_level+sem )) +
  theme_bw()+
  scale_color_viridis_d(option="viridis")+
  scale_fill_viridis_d(option="viridis")+
  geom_line() +
  geom_ribbon(alpha=0.3, colour = NA) +
  ylab("Methylation level")+
  labs(colour="Expression category")+
  ylim(0,1)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  guides(fill=FALSE)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))



load(file="COL_Heart_Meth_x_RNA_LMH.rda")
COL.Heart.all.res <- all.res
COL.Heart.all.res$Species <- "COL"

load(file="PIE_Heart_Meth_x_RNA_LMH.rda")
PIE.Heart.all.res <- all.res
PIE.Heart.all.res$Species <- "PIE"

load(file="HYB_Heart_Meth_x_RNA_LMH.rda")
HYB.Heart.all.res <- all.res
HYB.Heart.all.res$Species <- "HYB"

Heart.spec.all.res <- rbind(COL.Heart.all.res, PIE.Heart.all.res, HYB.Heart.all.res)

HY <- ggplot(data=Heart.spec.all.res[Heart.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Individual), col=Expression_cat, lty=Species, fill=Expression_cat,  ymin=Methylation_level-sem, ymax=Methylation_level+sem )) +
  theme_bw()+
  scale_color_viridis_d(option="viridis")+
  scale_fill_viridis_d(option="viridis")+
  geom_line() +
  geom_ribbon(alpha=0.3, colour = NA) +
  #ylab("Methylation level")+
  labs(colour="Expression category")+
  ylim(0,1)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  guides(fill=FALSE)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


load(file="COL_Kidney_Meth_x_RNA_LMH.rda")
COL.Kidney.all.res <- all.res
COL.Kidney.all.res$Species <- "COL"

load(file="PIE_Kidney_Meth_x_RNA_LMH.rda")
PIE.Kidney.all.res <- all.res
PIE.Kidney.all.res$Species <- "PIE"

load(file="HYB_Kidney_Meth_x_RNA_LMH.rda")
HYB.Kidney.all.res <- all.res
HYB.Kidney.all.res$Species <- "HYB"

Kidney.spec.all.res <- rbind(COL.Kidney.all.res, PIE.Kidney.all.res, HYB.Kidney.all.res)

KY <- ggplot(data=Kidney.spec.all.res[Kidney.spec.all.res$CGI_prom == "N",], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Individual), col=Expression_cat, lty=Species, fill=Expression_cat,  ymin=Methylation_level-sem, ymax=Methylation_level+sem )) +
  theme_bw()+
  scale_color_viridis_d(option="viridis")+
  scale_fill_viridis_d(option="viridis")+
  geom_line() +
  geom_ribbon(alpha=0.3, colour = NA) +
  ylab("Methylation level")+
  labs(colour="Expression category")+
  ylim(0,1)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  guides(fill=FALSE)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

Brain.spec.all.res$Tissue <- "Brain"
Heart.spec.all.res$Tissue <- "Heart"
Kidney.spec.all.res$Tissue <- "Kidney"
Liver.spec.all.res$Tissue <- "Liver"
Testis.spec.all.res$Tissue <- "Testis"

all.tiss.spec.all.res <- rbind(Brain.spec.all.res, Heart.spec.all.res, Kidney.spec.all.res, Liver.spec.all.res, Testis.spec.all.res)

library("plyr")
all.tiss.spec.all.res.mean <- ddply(all.tiss.spec.all.res, c("Segment_number", "Segment_type", "CGI_prom", "Expression_cat", "Species", "Tissue"),
                                    function(x) c(mean(x$Methylation_level, na.rm=T), mean(x$sem, na.rm=T)))

colnames(all.tiss.spec.all.res.mean) <- c("Segment_number", "Segment_type", "CGI_prom", "Expression_cat","Species","Tissue", "Methylation_level", "sem")

ggplot(data=all.tiss.spec.all.res.mean[all.tiss.spec.all.res.mean$CGI_prom == "Y",], aes(x=Segment_number, y=Methylation_level, group = interaction(Expression_cat, Tissue, Species), col=Tissue, fill=Tissue,  ymin=Methylation_level-sem, ymax=Methylation_level+sem )) +
  theme_classic()+
  #scale_color_manual(name="Expression category", values = c("#004D40", "#FFC107", "#D81B60")) +
  #scale_fill_manual(name="Species", values = c("#004D40", "#FFC107", "#D81B60")) +
  geom_line() +
  geom_ribbon(alpha=0.3, colour = NA) +
  ylab("Methylation level")+
  labs(colour="Expression category")+
  ylim(0,1)+
  guides(fill=FALSE)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

library(gridExtra)
grid.arrange(HY,KY,LY,TY, nrow=1)





# CvP DE Fst ####

library(plyr)
data_DE_status_Fst <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
load(paste(tissue, "DE_status_nonCpG_Fst.rda", sep="_"))
all.res$Tissue <- tissue
data_DE_status_Fst <- rbind(data_DE_status_Fst, all.res)
}

data_DE_status_Fst_all.tiss <- ddply(data_DE_status_Fst, c("Segment_number", "Segment_type", "DE_status"), function(x) c(mean(x$CpG_Fst, na.rm=T), mean(x$non_CpG_Fst, na.rm=T)))
colnames(data_DE_status_Fst_all.tiss) <- c("Segment_number", "Segment_type", "DE_status", "CpG_Fst", "non_CpG_Fst")

ggplot(data=data_DE_status_Fst_all.tiss, aes(x=Segment_number, y=non_CpG_Fst, col=DE_status))+
  theme_classic()+
  geom_line() +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_DE_status_Fst_all.tiss, aes(x=Segment_number, y=CpG_Fst, col=DE_status))+
  theme_classic()+
  geom_line() +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_DE_status_Fst, aes(x=Segment_number, y=CpG_Fst, col=Tissue, fill=Tissue, group = interaction(DE_status, Tissue), lty=DE_status))+
  theme_classic()+
  ylim(0, 0.3)+
  #geom_line() +
  ylab(expression("CpG" ~ italic("F"["ST"])))+
  #ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_DE_status_Fst, aes(x=Segment_number, y=non_CpG_Fst, col=Tissue, fill=Tissue, group = interaction(DE_status, Tissue), lty=DE_status))+
  theme_classic()+
  ylim(0, 0.3)+
  #geom_line() +
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  #ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

# HYB DE Fst ####

library(plyr)
data_DE_status_Fst_Hyb <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "hybridDE_status_Fst.rda", sep="_"))
  all.res$Tissue <- tissue
  data_DE_status_Fst_Hyb <- rbind(data_DE_status_Fst_Hyb, all.res)
}

ggplot(data=data_DE_status_Fst_Hyb, aes(x=Segment_number, y=CpG_Fst, col=Tissue, fill=Tissue, group = interaction(DE_status_HvC, Tissue), lty=DE_status_HvC))+
  theme_classic()+

  #geom_line() +
  ylab(expression("CpG" ~ italic("F"["ST"])))+
  #ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_DE_status_Fst_Hyb, aes(x=Segment_number, y=non_CpG_Fst, col=Tissue, fill=Tissue, group = interaction(DE_status_HvC, Tissue), lty=DE_status_HvC))+
  theme_classic()+

  #geom_line() +
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  #ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_DE_status_Fst_Hyb, aes(x=Segment_number, y=CpG_Fst, col=Tissue, fill=Tissue, group = interaction(DE_status_HvP, Tissue), lty=DE_status_HvP))+
  theme_classic()+
  #geom_line() +
  ylab(expression("CpG" ~ italic("F"["ST"])))+
  #ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_DE_status_Fst_Hyb, aes(x=Segment_number, y=non_CpG_Fst, col=Tissue, fill=Tissue, group = interaction(DE_status_HvP, Tissue), lty=DE_status_HvP))+
  theme_classic()+
  ylim(0, 0.3)+
  #geom_line() +
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  #ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))





# Cor Fst CvP ####
data_cor_Fst_Meth <- data.frame()
for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "CGIprom_cor_Fst_Meth.rda", sep="_"))
  all.res$Tissue <- tissue
  data_cor_Fst_Meth<- rbind(data_cor_Fst_Meth, all.res)
}
#data_cor_Fst_Meth[data_cor_Fst_Meth$Tissue == tissue,]
#tissue="Testis"
#expression(italic("F"["ST"]))
#colour="#E66100"
#"#5D3A9B"
data_cor_Fst_Meth$CGI_prom <- ifelse(data_cor_Fst_Meth$CGI_prom == "Y", "CGI", "Other")


ggplot(data=data_cor_Fst_Meth, aes(x=Segment_number, y=cor_CpG_Fst, lty=CGI_prom, group=interaction(Tissue, CGI_prom)))+
  
  theme_classic()+
  #geom_line(aes(colour="CpG"), alpha=0.15) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_smooth(data=data_cor_Fst_Meth, aes(x=Segment_number, y=cor_CpG_Fst, col="CpG", fill="CpG"), fill= "#E66100",method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.1)+
  #geom_line( aes(y=cor_non_CpG_Fst, lty=CGI_prom, colour="non-CpG"), alpha=0.15)+
  geom_smooth(data=data_cor_Fst_Meth, aes(x=Segment_number, y=cor_non_CpG_Fst, col="non-CpG"), fill="#5D3A9B",  method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.1)+
  ylab(expression("Spearman's" ~ rho))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  scale_color_manual(name=expression(italic("F"["ST"])), breaks=c("CpG", "non-CpG"), values = c("CpG"="#E66100", "non-CpG"="#5D3A9B"))+
  scale_linetype(name='Promoter type')+
  theme(aspect.ratio=1, legend.position="right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


#Tissue colouration
ggplot(data=data_cor_Fst_Meth, aes(x=Segment_number, y=cor_CpG_Fst, lty=CGI_prom, group=interaction(Tissue, CGI_prom)))+
  theme_classic()+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_smooth(data=data_cor_Fst_Meth, aes(x=Segment_number, y=cor_CpG_Fst, col=Tissue, fill= Tissue),method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.1)+
  geom_smooth(data=data_cor_Fst_Meth, aes(x=Segment_number, y=cor_non_CpG_Fst, col=Tissue,fill=Tissue),  method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.1)+
  ylab(expression("Spearman's" ~ rho))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  #scale_color_manual(name=expression(italic("F"["ST"])), breaks=c("CpG", "non-CpG"), values = c("CpG"="#E66100", "non-CpG"="#5D3A9B"))+
  scale_linetype(name='Promoter type')+
  theme(aspect.ratio=1, legend.position="right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))






data_cor_Fst_Meth_all.tiss <- ddply(data_cor_Fst_Meth, c("Segment_number", "Segment_type", "CGI_prom "), function(x) c(mean(x$cor_CpG_Fst, na.rm=T), mean(x$cor_non_CpG_Fst, na.rm=T)))
colnames(data_cor_Fst_Meth_all.tiss) <- c("Segment_number", "Segment_type", "CGI_prom", "cor_CpG_Fst", "cor_non_CpG_Fst")

ggplot(data=data_cor_Fst_Meth_all.tiss, aes(x=Segment_number, y=cor_CpG_Fst, lty=CGI_prom))+
  theme_classic()+
  
  geom_line(colour="#E66100", alpha=0.15) +
  geom_hline(yintercept=0, lty=2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_smooth(data=data_cor_Fst_Meth_all.tiss, aes(x=Segment_number, y=cor_CpG_Fst), col="#E66100", method = "loess", span = 0.15, method.args = list(degree=1))+
  geom_line( aes(y=cor_non_CpG_Fst, lty=CGI_prom), colour="#5D3A9B", alpha=0.15)+
  geom_smooth(data=data_cor_Fst_Meth_all.tiss, aes(x=Segment_number, y=cor_non_CpG_Fst), colour="#5D3A9B",  method = "loess", span = 0.15, method.args = list(degree=1))+
  ylab(expression("Spearman's" ~ rho))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


#CvP DMR v DEstat ####

data_CvP_DMR_DEstat <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "CvP_DMR_DEstat.rda", sep="_"))
  all.res$Tissue <- tissue
  data_CvP_DMR_DEstat <- rbind(data_CvP_DMR_DEstat, all.res)
}

for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
print(tissue)
print(t.test(data_CvP_DMR_DEstat[data_CvP_DMR_DEstat$Segment_number>50 & data_CvP_DMR_DEstat$Segment_number<150 & data_CvP_DMR_DEstat$DE_status == "DE" & data_CvP_DMR_DEstat$Tissue == tissue, ]$DMR_freq, 
             data_CvP_DMR_DEstat[data_CvP_DMR_DEstat$Segment_number>50 & data_CvP_DMR_DEstat$Segment_number<150 & data_CvP_DMR_DEstat$DE_status == "non-DE" & data_CvP_DMR_DEstat$Tissue == tissue,, ]$DMR_freq,paired=T))
}



#ymin=DMR_freq-sem, ymax=DMR_freq+sem
ggplot(data=data_CvP_DMR_DEstat, aes(x=Segment_number, y=DMR_freq, col=Tissue, fill=Tissue, group = interaction(DE_status, Tissue), lty=DE_status))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("DMR frequency")+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_CvP_DMR_DEstat, aes(x=Segment_number, y=DMR_freq, col=Tissue, fill=Tissue, lty=DE_status))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("DMR frequency")+
  facet_grid(~Tissue)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(strip.text.x = element_text(size = 20), aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=14),  legend.title=element_text(size=14))




#HvC DMR v DEstat ####

data_HvC_DMR_DEstat <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "HvC_DMR_DE_stat.rda", sep="_"))
  all.res$Tissue <- tissue
  data_HvC_DMR_DEstat <- rbind(data_HvC_DMR_DEstat, all.res)
}


#ymin=DMR_freq-sem, ymax=DMR_freq+sem
ggplot(data=data_HvC_DMR_DEstat, aes(x=Segment_number, y=DMR_freq, col=Tissue, fill=Tissue, group = interaction(DE_status_HvC, Tissue), lty=DE_status_HvC))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("DMR frequency")+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="nonr", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


ggplot(data=data_HvC_DMR_DEstat, aes(x=Segment_number, y=DMR_freq, col=Tissue, fill=Tissue, lty=DE_status_HvC))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("DMR frequency")+
  facet_grid(~Tissue)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(strip.text.x = element_text(size = 20), aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=14),  legend.title=element_text(size=14))


#HvP DMR v DEstat ####

data_HvP_DMR_DEstat <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "HvP_DMR_DE_stat.rda", sep="_"))
  all.res$Tissue <- tissue
  data_HvP_DMR_DEstat <- rbind(data_HvP_DMR_DEstat, all.res)
}


#ymin=DMR_freq-sem, ymax=DMR_freq+sem
ggplot(data=data_HvP_DMR_DEstat, aes(x=Segment_number, y=DMR_freq, col=Tissue, fill=Tissue, group = interaction(DE_status_HvP, Tissue), lty=DE_status_HvP))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("DMR frequency")+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))



#Exin ####

data_Exin<- data.frame()

for (species in c("COL", "PIE", "HYB")){
for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(species, tissue, "Exin.rda", sep="_"))
  all.res$Species <- species
  all.res$Tissue <- tissue
  data_Exin <- rbind(data_Exin, all.res)
}
}

ggplot(data=data_Exin[data_Exin$CGI_prom == "Y" & data_Exin$Tissue == "Kidney",], aes(x=Segment_number, y=DMR_frequency, fill=Exin_annotation, col=Exin_annotation, lty=Species, group = interaction(Exin_annotation, Tissue, Species, Individual)))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  #geom_hline(yintercept=0, lty=2) +
  #ylim(-0.5, 0.5)+
  #ylim(0, 1)+
  scale_color_manual(name="Annotation", values = c("black", "red", "orange", "black")) +
  scale_fill_manual(name="Annotation", values = c("black", "red", "orange", "black")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab(expression("Spearman's" ~ rho))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_Exin[data_Exin$CGI_prom == "N" & data_Exin$Tissue == "Testis",], aes(x=Segment_number, y=Meth, lty=Species, col=Exin_annotation,fill=Exin_annotation, group = interaction(Exin_annotation, Tissue, Individual)))+
  theme_classic()+
  #theme_bw()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_hline(yintercept=0, lty=2) +
 # ylim(-0.5, 0.5)+
  ylim(0, 1)+
  scale_color_manual(name="Annotation", values = c("black", "red", "orange", "black")) +
  scale_fill_manual(name="Annotation", values = c("black", "red", "orange", "black")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
#ylim(0,1)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

#

ggplot(data=data_Exin[data_Exin$CGI_prom == "N" & data_Exin$Tissue == "Brain" & data_Exin$Individual == "COL01" ,], aes(x=Segment_number, y=Meth.CvP, lty=Species, col=Exin_annotation,fill=Exin_annotation, group = interaction(Exin_annotation, Tissue, Individual)))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Annotation", values = c("black", "red", "orange", "black")) +
  scale_fill_manual(name="Annotation", values = c("black", "red", "orange", "black")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("Relative methylation diff.")+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))




# Exin - Fst ####
colnames(data_Exin) <- c("Segment_number", "Segment_type","Individual", "CGI_prom","DMR_frequency", "Meth.CvP","Meth", "non_CpG_Fst", "CpG_Fst", "Cor", "Exin_annotation", "Species", "Tissue")


ggplot(data=data_Exin[data_Exin$CGI_prom == "Y" & data_Exin$Individual == "COL01" & data_Exin$Tissue == "Heart" ,], aes(x=Segment_number, y=non_CpG_Fst, col=Exin_annotation,fill=Exin_annotation, group = interaction(Exin_annotation), lty="CGI"))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  #geom_hline(yintercept=0, lty=2)+
  scale_color_manual(name="Annotation", values = c("black", "red", "orange", "black"), labels=c("Downstream", "Exon", "Intron", "Upstream")) +
  scale_fill_manual(name="Annotation", values = c("black", "red", "orange", "black"), labels=c("Downstream", "Exon", "Intron", "Upstream")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  geom_smooth(data=data_Exin[data_Exin$CGI_prom == "N" & data_Exin$Individual == "COL01" & data_Exin$Tissue == "Heart" ,], aes(x=Segment_number, y=non_CpG_Fst, lty="Other"),  method = "loess", span = 0.15, method.args = list(degree=1),  alpha=0.2)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  scale_linetype_manual(name="Promoter type", breaks=c("Other", "CGI"), values=c("Other"=2, "CGI"=1))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_Exin[data_Exin$CGI_prom == "Y" & data_Exin$Individual == "COL01" & data_Exin$Tissue == "Heart" ,], aes(x=Segment_number, y=CpG_Fst, col=Exin_annotation,fill=Exin_annotation, group = interaction(Exin_annotation), lty="CGI"))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  #geom_hline(yintercept=0, lty=2)+
  scale_color_manual(name="Annotation", values = c("black", "red", "orange", "black"), labels=c("Downstream", "Exon", "Intron", "Upstream")) +
  scale_fill_manual(name="Annotation", values = c("black", "red", "orange", "black"), labels=c("Downstream", "Exon", "Intron", "Upstream")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab(expression("CpG" ~ italic("F"["ST"])))+
  geom_smooth(data=data_Exin[data_Exin$CGI_prom == "N" & data_Exin$Individual == "COL01" & data_Exin$Tissue == "Heart" ,], aes(x=Segment_number, y=CpG_Fst, lty="Other"),  method = "loess", span = 0.15, method.args = list(degree=1),  alpha=0.2)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  scale_linetype_manual(name="Promoter type", breaks=c("Other", "CGI"), values=c("Other"=2, "CGI"=1))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


# Exin - Methylation differentiation ####
#Exin colouration
ggplot(data=data_Exin[data_Exin$CGI_prom == "Y" & data_Exin$Individual == "COL01" ,], aes(x=Segment_number, y=Meth.CvP, col=Exin_annotation,fill=Exin_annotation, group = interaction(Exin_annotation, Tissue), lty="CGI"))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  #geom_hline(yintercept=0, lty=2)+
  ylab(expression(italic("M"["diff"])))+
  scale_color_manual(name="Annotation", values = c("black", "red", "orange", "black"), labels=c("Downstream", "Exon", "Intron", "Upstream")) +
  scale_fill_manual(name="Annotation", values = c("black", "red", "orange", "black"), labels=c("Downstream", "Exon", "Intron", "Upstream")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_smooth(data=data_Exin[data_Exin$CGI_prom == "N" & data_Exin$Individual == "COL01" ,], aes(x=Segment_number, y=Meth.CvP, lty="Other"),  method = "loess", span = 0.15, method.args = list(degree=1),  alpha=0.2)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  scale_linetype_manual(name="Promoter type", breaks=c("Other", "CGI"), values=c("Other"=2, "CGI"=1))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))



#Tissue colouration
ggplot(data=data_Exin[data_Exin$CGI_prom == "Y" & data_Exin$Individual == "COL01" ,], aes(x=Segment_number, y=Meth.CvP, col=Tissue, fill=Tissue, group = interaction(Exin_annotation, Tissue), lty="CGI"))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  #geom_hline(yintercept=0, lty=2)+
  ylab(expression(italic("M"["diff"])))+
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
 geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_smooth(data=data_Exin[data_Exin$CGI_prom == "N" & data_Exin$Individual == "COL01" ,], aes(x=Segment_number, y=Meth.CvP, lty="Other"),  method = "loess", span = 0.15, method.args = list(degree=1),  alpha=0.2)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  scale_linetype_manual(name="Promoter type", breaks=c("Other", "CGI"), values=c("Other"=2, "CGI"=1))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))





# Across tissue cor ####

load("COL_acrossTissCor_noBrain.rda")
COL.acrossCor <- all.res
mean(COL.acrossCor[COL.acrossCor$Segment_number == "50" & COL.acrossCor$CGI_prom == "Y",]$Cor, na.rm=T)
COL.acrossCor.mean <- aggregate(Cor ~ CGI_prom+Segment_number+Segment_type, COL.acrossCor, mean)
COL.acrossCor.mean$Species <- "COL"

load("HYB_acrossTissCor_noBrain.rda")
HYB.acrossCor <- all.res
mean(HYB.acrossCor[HYB.acrossCor$Segment_number == "50" & HYB.acrossCor$CGI_prom == "Y",]$Cor, na.rm=T)
HYB.acrossCor.mean <- aggregate(Cor ~ CGI_prom+Segment_number+Segment_type, HYB.acrossCor, mean)
HYB.acrossCor.mean$Species <- "HYB"


load("PIE_acrossTissCor_noBrain.rda")
PIE.acrossCor <- all.res
mean(PIE.acrossCor[PIE.acrossCor$Segment_number == "50" & PIE.acrossCor$CGI_prom == "Y",]$Cor, na.rm=T)
PIE.acrossCor.mean <- aggregate(Cor ~ CGI_prom+Segment_number+Segment_type, PIE.acrossCor, mean)
PIE.acrossCor.mean$Species <- "PIE"

acrossCor <- rbind(COL.acrossCor.mean, PIE.acrossCor.mean, HYB.acrossCor.mean)
acrossCor$CGI_prom <- ifelse(acrossCor$CGI_prom == "Y", "CGI", "Other")

ggplot(data=acrossCor, aes(x=Segment_number, y=Cor, col=Species, lty=CGI_prom )) +
  theme_classic()+
  scale_color_manual(name="Species", values = c("red", "darkgoldenrod3", "blue")) +
  geom_line() +
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  geom_hline(yintercept=0, lty=2)+
  ylab(expression("Spearman's" ~ rho))+
  #geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1))+
  #ylim(-0.5,0.5)+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  scale_linetype(name="Promoter type")+
  
  theme(aspect.ratio=1, legend.position="right",  axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


