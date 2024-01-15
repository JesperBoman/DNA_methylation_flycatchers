#### PART 2 ####
install.packages("BiSeq")
library("BiSeq")
library("plyr")
tissue="Heart" #Fill in tissue here

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
    
    if( is.na(rej.COLvPIE.par) == T ){COLvPIE.par.test=0} 
    else {COLvPIE.par.test=length(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]$X)}
    
    if( is.na(rej.COLvPIE.hyb) == T ){COLvPIE.hyb.test=0} 
    else {COLvPIE.hyb.test=length(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]$X)}
    
    if( is.na(rej.COL.comb) == T ){COL.comb.test=0} 
    else {COL.comb.test=length(rej.COL.comb[rej.COL.comb$X ==  region,]$X)}
    
    if( is.na(rej.PIE.comb) == T ){PIE.comb.test=0} 
    else {PIE.comb.test=length(rej.PIE.comb[rej.PIE.comb$X ==  region,]$X)}

    
#    H_C_v_P <-sum(meth_annot_data_reg[meth_annot_data_reg$Type == "Hybrid-COL" & meth_annot_data_reg$region_ID == region & meth_annot_data_reg$Tissue == tissue, ]$Mean_of_ratios_per_dinuc, na.rm=T) - sum(meth_annot_data_reg[meth_annot_data_reg$Type == "Hybrid-PIE" & meth_annot_data_reg$region_ID == region & meth_annot_data_reg$Tissue == tissue, ]$Mean_of_ratios_per_dinuc, na.rm=T)
#    P_C_v_P <-sum(meth_annot_data_reg[meth_annot_data_reg$Type == "Parental-COL" & meth_annot_data_reg$region_ID == region & meth_annot_data_reg$Tissue == tissue, ]$Mean_of_ratios_per_dinuc, na.rm=T) - sum(meth_annot_data_reg[meth_annot_data_reg$Type == "Parental-PIE" & meth_annot_data_reg$region_ID == region & meth_annot_data_reg$Tissue == tissue, ]$Mean_of_ratios_per_dinuc, na.rm=T)
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
    if( is.na(rej.COLvPIE.par) == T){COLvPIE.par.test=0} 
    else {COLvPIE.par.test=length(rej.COLvPIE.par[rej.COLvPIE.par$X ==  region,]$X)}
    if( is.na(rej.COLvPIE.hyb) == T){COLvPIE.hyb.test=0} 
    else {COLvPIE.hyb.test=length(rej.COLvPIE.hyb[rej.COLvPIE.hyb$X ==  region,]$X)}
    if( is.na(rej.COL.comb) == T ){COL.comb.test=0} 
    else {COL.comb.test=length(rej.COL.comb[rej.COL.comb$X ==  region,]$X)}
    if( is.na(rej.PIE.comb) == T ){PIE.comb.test=0} 
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
