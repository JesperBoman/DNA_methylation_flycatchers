#### PART 2 ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiSeq")
library("BiSeq")
library("plyr")
library("lmodel2")
library("ggplot2")

FDR.level <- 0.1


reg_bed_file <- read.table("fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_plus_minus_100bp.bed")


regulatory_class_data_stat_all_tiss <- data.frame()
regulatory_class_data_stat_freq_all_tiss <- data.frame()


for (tissue in c("Heart", "Brain", "Testis", "Liver" )){ #,"Heart", "Liver", "Brain", "Testis"
  
  load(paste(tissue, "BiSeq_regulatory_mechanism_data.rda", sep="_"))
  
  #COLvPIE.par
  if (length(try(testClusters(locCor_COLvPIE.par, FDR.cluster = FDR.level), silent=T)) == 1){
    clusters.rej.COLvPIE.par<-NA
    not.rej.COLvPIE.par <-NA
    rej.COLvPIE.par  <- NA
  } else{
    clusters.rej.COLvPIE.par <- testClusters(locCor_COLvPIE.par, FDR.cluster = FDR.level)
    not.rej.COLvPIE.par <- as.data.frame(clusters.rej.COLvPIE.par$clusters.not.reject)
    rej.COLvPIE.par <- as.data.frame(clusters.rej.COLvPIE.par$clusters.reject)
  }
  
  #COLvPIE.hyb
  if (length(try(testClusters(locCor_COLvPIE.hyb, FDR.cluster = FDR.level), silent=T)) == 1){
    clusters.rej.COLvPIE.hyb<-NA
    not.rej.COLvPIE.hyb <- NA
    rej.COLvPIE.hyb <-NA
  } else{
    clusters.rej.COLvPIE.hyb <- testClusters(locCor_COLvPIE.hyb, FDR.cluster = FDR.level)
    not.rej.COLvPIE.hyb <- as.data.frame(clusters.rej.COLvPIE.hyb$clusters.not.reject)
    rej.COLvPIE.hyb <- as.data.frame(clusters.rej.COLvPIE.hyb$clusters.reject)
  }
  
  #COL(par) v COL(hyb)
  if (length(try(testClusters(locCor_COL.comb, FDR.cluster = FDR.level), silent=T)) == 1){
    clusters.rej.COL.comb<-NA
    not.rej.COL.comb <- NA
    rej.COL.comb <-NA
  } else{
    clusters.rej.COL.comb <- testClusters(locCor_COL.comb, FDR.cluster = FDR.level)
    not.rej.COL.comb <- as.data.frame(clusters.rej.COL.comb$clusters.not.reject)
    rej.COL.comb<- as.data.frame(clusters.rej.COL.comb$clusters.reject)
  }
  
  #PIE(par) v PIE(hyb)
  if (length(try(testClusters(locCor_PIE.comb, FDR.cluster = FDR.level), silent=T)) == 1){
    clusters.rej.PIE.comb<-NA
    not.rej.PIE.comb <- NA
    rej.PIE.comb <-NA
  } else{
    clusters.rej.PIE.comb <- testClusters(locCor_PIE.comb, FDR.cluster = FDR.level)
    not.rej.PIE.comb <- as.data.frame(clusters.rej.PIE.comb$clusters.not.reject)
    rej.PIE.comb<- as.data.frame(clusters.rej.PIE.comb$clusters.reject)
  }
  
  
  
  
  
  
  ambig.stat <- NA
  cons.stat <- NA
  cons.shifted.level.stat <- NA
  trans.stat <- NA
  cis.stat <- NA
  cis_plus_trans.stat <- NA
  compensatory.stat <- NA
  
  for(region in paste(reg_bed_file$V1, reg_bed_file$V2, reg_bed_file$V3, sep="_")){
    
    if( empty(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]) == T){COLvPIE.par.test=0} 
    else {COLvPIE.par.test=length(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]$X)}
    
    if( empty(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]) == T ){COLvPIE.hyb.test=0} 
    else {COLvPIE.hyb.test=length(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]$X)}
    
    if( empty(rej.COL.comb[rej.COL.comb$X ==  region,]) == T ){COL.comb.test=0} 
    else {COL.comb.test=length(rej.COL.comb[rej.COL.comb$X ==  region,]$X)}
    
    if( empty(rej.PIE.comb[rej.PIE.comb$X ==  region,]) == T ){PIE.comb.test=0} 
    else {PIE.comb.test=length(rej.PIE.comb[rej.PIE.comb$X ==  region,]$X)}

    
  if( empty(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]) == T & empty(not.rej.COLvPIE.par[not.rej.COLvPIE.par$X ==  region,]) == T | empty(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]) == T & empty(not.rej.COLvPIE.hyb[not.rej.COLvPIE.hyb$X ==  region,]) == T ){
      ambig.stat <- c(ambig.stat, region)

    }
    else{
    
    
    if(COLvPIE.par.test == 0 & COLvPIE.hyb.test == 0 ){
      if(COL.comb.test != 0 & PIE.comb.test != 0 ){
        cons.shifted.level.stat <- c(cons.shifted.level.stat, region)
      }
      else{
        cons.stat <- c(cons.stat, region)
      }

    }
    
    if(COLvPIE.par.test != 0 & COLvPIE.hyb.test == 0 & (COL.comb.test != 0 | PIE.comb.test != 0)){
      trans.stat <- c(trans.stat, region)

    }
    
    if(COLvPIE.par.test == 0 & COLvPIE.hyb.test != 0 & (COL.comb.test != 0 | PIE.comb.test != 0)){
      compensatory.stat <- c(compensatory.stat, region)
      
    }
    if(COLvPIE.par.test != 0 & COLvPIE.hyb.test != 0 & COL.comb.test == 0 & PIE.comb.test == 0 ){
      cis.stat <- c(cis.stat, region)
      
    }
    if(COLvPIE.par.test != 0 & COLvPIE.hyb.test != 0 & (COL.comb.test != 0 | PIE.comb.test != 0)){
        cis_plus_trans.stat <- c(cis_plus_trans.stat, region)
        
    }
    }
    
  }  

  ambig.stat <- as.data.frame(ambig.stat[-1])
  colnames(ambig.stat)<-"region_ID"
  ambig.stat$Regpat <- "Ambiguous"
  
  cons.stat <- as.data.frame(cons.stat[-1])
  colnames(cons.stat)<-"region_ID"
  cons.stat$Regpat <- "Conserved"
  
  cons.shifted.level.stat <- as.data.frame(cons.shifted.level.stat[-1])
  colnames(cons.shifted.level.stat)<-"region_ID"
  try(cons.shifted.level.stat$Regpat <- "Conserved_shifted_level", silent = T)
  
  trans.stat <- as.data.frame(trans.stat[-1])
  colnames(trans.stat)<-"region_ID"
  try(trans.stat$Regpat  <- "Trans", silent = T)
  
  compensatory.stat <- as.data.frame(compensatory.stat[-1])
  colnames(compensatory.stat)<-"region_ID"
  try(compensatory.stat$Regpat  <- "Compensatory", silent = T)
  
  cis.stat <- as.data.frame(cis.stat[-1])
  colnames(cis.stat)<-"region_ID"
  try(cis.stat$Regpat  <- "Cis", silent = T)
  
  cis_plus_trans.stat <- as.data.frame(cis_plus_trans.stat[-1])
  colnames(cis_plus_trans.stat)<-"region_ID"
  try(cis_plus_trans.stat$Regpat  <- "Cis + Trans", silent = T)
  
  
  regulatory_class_data_stat <- rbind(ambig.stat, cons.stat, cons.shifted.level.stat, trans.stat, compensatory.stat, cis.stat, cis_plus_trans.stat)
  regulatory_class_data_stat$Chromosome <- sub('(Chr.*)_.*_.*', '\\1', regulatory_class_data_stat$region_ID)
  regulatory_class_data_stat$Chr_type <- ifelse(regulatory_class_data_stat$Chromosome == "ChrZ", "Z", "A") 
  regulatory_class_data_stat$Tissue <- tissue
  
  print(paste("The number of cons.shifted.level for ", tissue, " :", sep=""))
  print(length(cons.shifted.level.stat$Regpat))

  #The mean of the differences per CpG for each for BiSeq
  regulatory_class_data_stat$P_C_v_P_biseq <- NA
  regulatory_class_data_stat$H_C_v_P_biseq <- NA
  regulatory_class_data_stat$COL_Par_v_Hyb_biseq <- NA
  regulatory_class_data_stat$PIE_Par_v_Hyb_biseq <- NA
  
  for(region in regulatory_class_data_stat[regulatory_class_data_stat$Regpat != "Ambiguous",]$region_ID){
    if( empty(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]) == T){COLvPIE.par.test=0} 
    else {COLvPIE.par.test=length(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]$X)}
    
    if( empty(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]) == T ){COLvPIE.hyb.test=0} 
    else {COLvPIE.hyb.test=length(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]$X)}
    
    if( empty(rej.COL.comb[rej.COL.comb$X ==  region,]) == T ){COL.comb.test=0} 
    else {COL.comb.test=length(rej.COL.comb[rej.COL.comb$X ==  region,]$X)}
    
    if( empty(rej.PIE.comb[rej.PIE.comb$X ==  region,]) == T ){PIE.comb.test=0} 
    else {PIE.comb.test=length(rej.PIE.comb[rej.PIE.comb$X ==  region,]$X)}
    
    
    
    if(COLvPIE.par.test!=0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$P_C_v_P_biseq <- mean(clusters.rej.COLvPIE.par$CpGs.clust.reject[region][[1]][4]$meth.group1 - clusters.rej.COLvPIE.par$CpGs.clust.reject[region][[1]][5]$meth.group2)}
    if(COLvPIE.hyb.test!=0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$H_C_v_P_biseq <- mean(clusters.rej.COLvPIE.hyb$CpGs.clust.reject[region][[1]][4]$meth.group1 - clusters.rej.COLvPIE.hyb$CpGs.clust.reject[region][[1]][5]$meth.group2)}
    if(COLvPIE.par.test==0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$P_C_v_P_biseq <- mean(clusters.rej.COLvPIE.par$CpGs.clust.not.reject[region][[1]][4]$meth.group1 - clusters.rej.COLvPIE.par$CpGs.clust.not.reject[region][[1]][5]$meth.group2)}
    if(COLvPIE.hyb.test==0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$H_C_v_P_biseq <- mean(clusters.rej.COLvPIE.hyb$CpGs.clust.not.reject[region][[1]][4]$meth.group1 - clusters.rej.COLvPIE.hyb$CpGs.clust.not.reject[region][[1]][5]$meth.group2)}
 
    if(COL.comb.test!=0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$COL_Par_v_Hyb_biseq <- mean(clusters.rej.COL.comb$CpGs.clust.reject[region][[1]][4]$meth.group1 - clusters.rej.COL.comb$CpGs.clust.reject[region][[1]][5]$meth.group2)}
    if(PIE.comb.test!=0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$PIE_Par_v_Hyb_biseq <- mean(clusters.rej.PIE.comb$CpGs.clust.reject[region][[1]][4]$meth.group1 - clusters.rej.PIE.comb$CpGs.clust.reject[region][[1]][5]$meth.group2)}
    if(COL.comb.test==0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$COL_Par_v_Hyb_biseq <- mean(clusters.rej.COL.comb$CpGs.clust.not.reject[region][[1]][4]$meth.group1 - clusters.rej.COL.comb$CpGs.clust.not.reject[region][[1]][5]$meth.group2)}
    if(PIE.comb.test==0){regulatory_class_data_stat[regulatory_class_data_stat$region_ID == region,]$PIE_Par_v_Hyb_biseq <- mean(clusters.rej.PIE.comb$CpGs.clust.not.reject[region][[1]][4]$meth.group1 - clusters.rej.PIE.comb$CpGs.clust.not.reject[region][[1]][5]$meth.group2)}
    
 }
  
  regulatory_class_data_stat$Regpat <- ifelse(regulatory_class_data_stat$Regpat == "Cis" & (regulatory_class_data_stat$P_C_v_P_biseq*regulatory_class_data_stat$H_C_v_P_biseq < 0), "Cis x Trans", regulatory_class_data_stat$Regpat)
  regulatory_class_data_stat$Regpat <- ifelse(regulatory_class_data_stat$Regpat == "Cis + Trans" & (regulatory_class_data_stat$P_C_v_P_biseq*regulatory_class_data_stat$H_C_v_P_biseq < 0), "Cis x Trans", regulatory_class_data_stat$Regpat)
  
  regulatory_class_data_stat_all_tiss<- regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue != tissue,]
  
  regulatory_class_data_stat_all_tiss <- rbind(regulatory_class_data_stat_all_tiss, regulatory_class_data_stat)
  
  
  regulatory_class_data_stat_freq <- as.data.frame(table(regulatory_class_data_stat$Regpat, regulatory_class_data_stat$Chr_type))
  colnames(regulatory_class_data_stat_freq) <- c("Regpat", "Chr_type", "Freq")
  
  regulatory_class_data_stat_freq$Tissue <- tissue
  
  regulatory_class_data_stat_freq_all_tiss<- regulatory_class_data_stat_freq_all_tiss[regulatory_class_data_stat_freq_all_tiss$Tissue != tissue,]
  
  regulatory_class_data_stat_freq_all_tiss <- rbind(regulatory_class_data_stat_freq_all_tiss, regulatory_class_data_stat_freq)
}

regulatory_class_data_stat$P_C_v_P_biseq <- ifelse(is.nan(regulatory_class_data_stat$P_C_v_P_biseq), NA, regulatory_class_data_stat$P_C_v_P_biseq)
regulatory_class_data_stat$H_C_v_P_biseq <- ifelse(is.nan(regulatory_class_data_stat$H_C_v_P_biseq), NA, regulatory_class_data_stat$H_C_v_P_biseq)
head(regulatory_class_data_stat)



regulatory_class_data_stat_all_tiss$Regpat <- factor(regulatory_class_data_stat_all_tiss$Regpat ,levels = c("Trans", "Cis", "Compensatory", "Cis x Trans", "Cis + Trans", "Conserved", "Conserved_shifted_level", "Ambiguous"))
head(regulatory_class_data_stat_all_tiss)

regulatory_class_data_stat_freq_all_tiss$Regpat <- factor(regulatory_class_data_stat_freq_all_tiss$Regpat ,levels = c("Trans", "Cis", "Compensatory", "Cis x Trans", "Cis + Trans", "Conserved", "Conserved_shifted_level", "Ambiguous"))



#COL v PIE plot
tissue="Heart"
ggplot(data=regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level", ], aes(y=H_C_v_P_biseq, x=P_C_v_P_biseq, col=Regpat))+geom_point(size=5)+
  scale_color_manual("Reg. pattern", values=c("#C3D898", "#FB4B4E", "#136F63", "turquoise", "yellow", "black", "grey"))+
  xlim(-1, 1)+
  ylim(-1, 1)+
  theme_classic()+
  ylab("Hybrid: COL - PIE")+
  xlab("Parental: COL - PIE")+
  geom_hline(yintercept=0, alpha = 0.7)+
  geom_vline(xintercept=0, alpha = 0.7)+
  theme(aspect.ratio=1, legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


#Frequency per regulatory pattern plot
regred <- regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level",]
regred_prop <- as.data.frame(prop.table(table(regred$Regpat, regred$Tissue, regred$Chr_type),c(2,3)))
colnames(regred_prop) <- c("Regpat", "Tissue", "Chr_type", "Freq")

ggplot(regred_prop[regred_prop$Regpat != "Conserved",], aes(fill=Regpat, y=Freq, x=Tissue)) + 
  theme_classic()+
  scale_fill_manual("Reg. pattern", values=c("#C3D898", "#FB4B4E", "#136F63", "turquoise", "yellow", "black", "grey"))+
  geom_bar(position = "stack", stat="identity", colour="black")+
  ylab("Proportion")+
  facet_wrap(~Chr_type)+
  theme(strip.text = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=16))



#Pie chart
library(scales)

all_tiss <- as.data.frame(table(regulatory_class_data_stat_all_tiss$Regpat))
colnames(all_tiss) <- c("Regpat", "Freq")

bp<- ggplot(all_tiss[all_tiss$Regpat != "Ambiguous"   & all_tiss$Regpat != "Conserved_shifted_level",], aes(x="", y=Freq, fill=Regpat))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  theme(legend.text=element_text(size=22),  legend.title=element_text(size=24))+

  #scale_fill_manual("Reg. pattern", values=c("red", "blue", "turquoise", "orange","black", "grey"))
  scale_fill_manual("Reg. pattern", values=c("#C3D898", "#FB4B4E", "#136F63", "turquoise", "yellow", "black", "pink", "grey"))






#Inheritance plot 
library(lmodel2)
tissue="Heart"

inmod <- lmodel2(COL_Par_v_Hyb_biseq~PIE_Par_v_Hyb_biseq, data=regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue  & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level",])


ggplot(data=regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level", ], aes(y=COL_Par_v_Hyb_biseq, x=PIE_Par_v_Hyb_biseq, col=Regpat))+geom_point(size=5)+
  scale_color_manual("Reg. pattern", values=c("#C3D898", "#FB4B4E", "#136F63", "turquoise", "yellow", "black", "grey"))+
  xlim(-1, 1)+
  ylim(-1, 1)+
  theme_classic()+
  ylab("COL: Parental - Hybrid")+
  xlab("PIE: Parental - Hybrid")+
  geom_hline(yintercept=0, alpha = 0.7)+
  geom_vline(xintercept=0, alpha = 0.7)+
  annotate("text",x=0.55, y= 0.8, size=10, label=bquote(italic(p) ~ "≈" ~ .(round(inmod$P.param,2))), col="black")+
  geom_abline(slope=inmod$regression.results[2,3], intercept=inmod$regression.results[2,2], col="black", size=1)+ #Major axis regression
 theme(aspect.ratio=1, legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))



#Follow-up stats ####

#Between-tissue
fisher.test(regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level", ]$Regpat, regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level", ]$Tissue, simulate.p.value=TRUE, B=1e6)

#Between-tissue: divergent
fisher.test(regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved", ]$Regpat, regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level"  & regulatory_class_data_stat_all_tiss$Regpat != "Conserved", ]$Tissue, simulate.p.value=TRUE, B=1e6)

tissue="Testis"
#Between chromosome types for a given tissue
fisher.test(regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level", ]$Regpat, regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level", ]$Chr_type)

#Between chromosome types for a given tissue: divergent
fisher.test(regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved", ]$Regpat, regulatory_class_data_stat_all_tiss[regulatory_class_data_stat_all_tiss$Tissue == tissue & regulatory_class_data_stat_all_tiss$Regpat != "Ambiguous" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved_shifted_level" & regulatory_class_data_stat_all_tiss$Regpat != "Conserved", ]$Chr_type)



#Cis vs trans
noCons <- regulatory_class_data_stat_freq_all_tiss[regulatory_class_data_stat_freq_all_tiss$Regpat != "Conserved" & regulatory_class_data_stat_freq_all_tiss$Regpat != "Conserved_shifted_level" & regulatory_class_data_stat_freq_all_tiss$Regpat != "Ambiguous",]
df <- data.frame()
mat2 <- matrix(nrow=0, ncol=3)

for(tissue in c("Brain", "Heart", "Liver", "Testis")){
cis_A_n <- sum(noCons[noCons$Tissue == tissue & noCons$Regpat == "Cis" & noCons$Chr_type == "A",3])
trans_A_n <- sum(noCons[noCons$Tissue == tissue & noCons$Regpat == "Trans" & noCons$Chr_type == "A",3])

print(tissue)
print(binom.test(cis_A_n, cis_A_n+trans_A_n))

cis_Z_n <- sum(noCons[noCons$Tissue == tissue & noCons$Regpat == "Cis" & noCons$Chr_type == "Z",3])
trans_Z_n <- sum(noCons[noCons$Tissue == tissue & noCons$Regpat == "Trans" & noCons$Chr_type == "Z",3])

print(tissue)
print(binom.test(cis_Z_n, cis_Z_n+trans_Z_n))


cis_n <- sum(noCons[noCons$Tissue == tissue & noCons$Regpat == "Cis",3])
trans_n <- sum(noCons[noCons$Tissue == tissue & noCons$Regpat == "Trans",3])

print(tissue)
bin<-(binom.test(cis_n, cis_n+trans_n))
print(bin)



mat <- matrix(c(trans_A_n, cis_A_n, trans_Z_n, cis_Z_n), nrow=2, ncol=2)

fp <- fisher.test(mat)
df <- rbind(df, cbind(Tissue=tissue, or=fp$estimate, pval=fp$p.value, totRat_low=0.5* 1/bin$conf.int[1], totRat=trans_n/cis_n, totRat_high=0.5* 1/bin$conf.int[2]) )
mat <- matrix(c(trans_A_n, cis_A_n, trans_Z_n, cis_Z_n, tissue, tissue), nrow=2, ncol=3)
mat2 <- rbind(mat2, mat)
}
