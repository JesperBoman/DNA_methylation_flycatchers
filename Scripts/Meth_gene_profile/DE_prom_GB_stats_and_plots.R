#CvP DEstat ####

data_CvP_all_DEstat <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "CvP_DE_status_CGI_all.rda", sep="_"))
  all.res$Tissue <- tissue
  data_CvP_all_DEstat <- rbind(data_CvP_all_DEstat, all.res)
}

#CvP - DMR_freq ####

#Promoter
CGI_prom="N"
for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t<-t.test(data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>30 & data_CvP_all_DEstat$Segment_number<51 & data_CvP_all_DEstat$DE_status == "DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$CvP_DMR_freq, 
               data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>30 & data_CvP_all_DEstat$Segment_number<51 & data_CvP_all_DEstat$DE_status == "non-DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$CvP_DMR_freq,paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=20))
  }



#Gene body
CGI_prom="Y"
for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>50 & data_CvP_all_DEstat$Segment_number<150 & data_CvP_all_DEstat$DE_status == "DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$CvP_DMR_freq, 
               data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>50 & data_CvP_all_DEstat$Segment_number<150 & data_CvP_all_DEstat$DE_status == "non-DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$CvP_DMR_freq,paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=20))
  
  }


#CvP - non-CpG Fst ####
#Promoter
CGI_prom="N"
for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>30 & data_CvP_all_DEstat$Segment_number<51 & data_CvP_all_DEstat$DE_status == "DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
               data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>30 & data_CvP_all_DEstat$Segment_number<51 & data_CvP_all_DEstat$DE_status == "non-DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=20))
  }


#Gene body
CGI_prom="N"
for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>50 & data_CvP_all_DEstat$Segment_number<150 & data_CvP_all_DEstat$DE_status == "DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
               data_CvP_all_DEstat[data_CvP_all_DEstat$Segment_number>50 & data_CvP_all_DEstat$Segment_number<150 & data_CvP_all_DEstat$DE_status == "non-DE" & data_CvP_all_DEstat$Tissue == tissue & data_CvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=20))
}
  
  

  
  #
  data_CvP_all_DEstat$Promoter_type <- ifelse(data_CvP_all_DEstat$CGI_prom == "Y", "CGI", "Other")
  
  ggplot(data=data_CvP_all_DEstat[data_CvP_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=CvP_DMR_freq  , col=Tissue, fill=Tissue, group = interaction(DE_status, Tissue), lty=DE_status))+
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
  
  
  #Control plot of CGI CvP DMR and CGI proportion per segment
  
  head(data_CvP_all_DEstat)
  
  ggplot(data=data_CvP_all_DEstat[data_CvP_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=CvP_DMR_freq  , col=Tissue, fill=Tissue, group = interaction(DE_status, Tissue), lty=DE_status))+
    theme_classic()+
    scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
    scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
    geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
    geom_vline(xintercept=30, alpha=0.3)+
    geom_vline(xintercept=50, alpha=0.3)+
    geom_vline(xintercept=150, alpha=0.3)+
    ylab("DMR frequency")+
    geom_line(data=all.agg.cgi.full, aes(x=Segment_number, y=Prop), col="black")+
    scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
    theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))
  

  ggplot(data=data_CvP_all_DEstat[data_CvP_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=CpG_Fst  , col=Tissue, fill=Tissue, group = interaction(DE_status, Tissue), lty=DE_status))+
    theme_classic()+
    #geom_ribbon(alpha=0.3, colour=NA)+
    #geom_line() +
    scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
    scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
    geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
    geom_vline(xintercept=50, alpha=0.3)+
    geom_vline(xintercept=150, alpha=0.3)+
    ylab(expression("CpG" ~ italic("F"["ST"])))+
    scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
    theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))
  
  
  
  
  
  data_Exin_DE<- data.frame()
  
  for (species in c("COL", "PIE", "HYB")){
    for (tissue in c("Heart", "Kidney", "Liver", "Testis", "Brain")){
      load(paste(species, tissue, "Exin_CvP_DE.rda", sep="_"))
      all.res$Species <- species
      all.res$Tissue <- tissue
      data_Exin_DE <- rbind(data_Exin_DE, all.res)
    }
  }
  
  
  ggplot(data=data_Exin_DE[data_Exin_DE$CGI_prom == "Y",], aes(x=Segment_number, y=non_CpG_Fst, fill=Exin_annotation, col=Exin_annotation, lty=DE_status, group = interaction(Exin_annotation, Tissue, Species, Individual, DE_status)))+
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
    facet_wrap(~Tissue)+
    scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
    theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))
  

  
  
  #HvC DEstat ####
  
  data_HvC_all_DEstat <- data.frame()
  
  for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
    load(paste(tissue, "HvC_DE_status_CGI_all.rda", sep="_"))
    all.res$Tissue <- tissue
    data_HvC_all_DEstat <- rbind(data_HvC_all_DEstat, all.res)
  }
  
  #HvC - DMR_freq ####
  data_HvC_all_DEstat<-data_HvC_all_DEstat[!is.na(data_HvC_all_DEstat$DE_status_HvC),]
  
  #Promoter
  CGI_prom="N"
  for ( tissue in c("Brain","Heart", "Kidney", "Liver", "Testis")){
    print(tissue)
    t<-t.test(data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>30 & data_HvC_all_DEstat$Segment_number<51 & data_HvC_all_DEstat$DE_status_HvC == "DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq, 
              data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>30 & data_HvC_all_DEstat$Segment_number<51 & data_HvC_all_DEstat$DE_status_HvC == "non-DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq,paired=T)
    print(t$estimate)
    print(t$p.value)
  }
  
  
  
  #Gene body
  CGI_prom="N"
  for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
    print(tissue)
    t <- t.test(data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>50 & data_HvC_all_DEstat$Segment_number<150 & data_HvC_all_DEstat$DE_status_HvC == "DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq, 
                data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>50 & data_HvC_all_DEstat$Segment_number<150 & data_HvC_all_DEstat$DE_status_HvC == "non-DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq,paired=T)
    print(t$estimate)
    print(t$p.value)
  }
  
  
  #HvC - non-CpG Fst ####
  #Promoter
  CGI_prom="N"
  for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
    print(tissue)
    t <- t.test(data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>30 & data_HvC_all_DEstat$Segment_number<51 & data_HvC_all_DEstat$DE_status_HvC == "DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
                data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>30 & data_HvC_all_DEstat$Segment_number<51 & data_HvC_all_DEstat$DE_status_HvC == "non-DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
    print(t$estimate)
    print(t$p.value)
  }
  
  
  #Gene body
  CGI_prom="N"
  for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
    print(tissue)
    t <- t.test(data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>50 & data_HvC_all_DEstat$Segment_number<150 & data_HvC_all_DEstat$DE_status_HvC == "DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
                data_HvC_all_DEstat[data_HvC_all_DEstat$Segment_number>50 & data_HvC_all_DEstat$Segment_number<150 & data_HvC_all_DEstat$DE_status_HvC == "non-DE" & data_HvC_all_DEstat$Tissue == tissue & data_HvC_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
    print(t$estimate)
    print(t$p.value)
  }
  }



#
data_HvC_all_DEstat$Promoter_type <- ifelse(data_HvC_all_DEstat$CGI_prom == "Y", "CGI", "Other")

ggplot(data=data_HvC_all_DEstat[data_HvC_all_DEstat$Promoter_type == "Other",], aes(x=Segment_number, y=HvC_DMR_freq  , col=Tissue, fill=Tissue, group = interaction(DE_status_HvC , Tissue), lty=DE_status_HvC ))+
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


ggplot(data=data_HvC_all_DEstat[data_HvC_all_DEstat$Promoter_type == "Other",], aes(x=Segment_number, y=non_CpG_Fst  , col=Tissue, fill=Tissue, group = interaction(DE_status_HvC , Tissue), lty=DE_status_HvC ))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))











#HvP DEstat ####

data_HvP_all_DEstat <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "HvP_DE_status_CGI_all.rda", sep="_"))
  all.res$Tissue <- tissue
  data_HvP_all_DEstat <- rbind(data_HvP_all_DEstat, all.res)
}

#HvP - DMR_freq ####
data_HvP_all_DEstat<-data_HvP_all_DEstat[!is.na(data_HvP_all_DEstat$DE_status_HvP),]

#Promoter
CGI_prom="N"
for ( tissue in c("Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t<-t.test(data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>30 & data_HvP_all_DEstat$Segment_number<51 & data_HvP_all_DEstat$DE_status_HvP == "DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq, 
            data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>30 & data_HvP_all_DEstat$Segment_number<51 & data_HvP_all_DEstat$DE_status_HvP == "non-DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq,paired=T)
  print(t$estimate)
  print(t$p.value)
}



#Gene body
CGI_prom="N"
for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>50 & data_HvP_all_DEstat$Segment_number<150 & data_HvP_all_DEstat$DE_status_HvP == "DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq, 
              data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>50 & data_HvP_all_DEstat$Segment_number<150 & data_HvP_all_DEstat$DE_status_HvP == "non-DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq,paired=T)
  print(t$estimate)
  print(t$p.value)
}


#HvP - non-CpG Fst ####
#Promoter
CGI_prom="N"
for ( tissue in c("Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>30 & data_HvP_all_DEstat$Segment_number<51 & data_HvP_all_DEstat$DE_status_HvP == "DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
              data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>30 & data_HvP_all_DEstat$Segment_number<51 & data_HvP_all_DEstat$DE_status_HvP == "non-DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
  print(t$estimate)
  print(t$p.value)
}


#Gene body
CGI_prom="N"
for ( tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>50 & data_HvP_all_DEstat$Segment_number<150 & data_HvP_all_DEstat$DE_status_HvP == "DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
              data_HvP_all_DEstat[data_HvP_all_DEstat$Segment_number>50 & data_HvP_all_DEstat$Segment_number<150 & data_HvP_all_DEstat$DE_status_HvP == "non-DE" & data_HvP_all_DEstat$Tissue == tissue & data_HvP_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
  print(t$estimate)
  print(t$p.value)
}



#
data_HvP_all_DEstat$Promoter_type <- ifelse(data_HvP_all_DEstat$CGI_prom == "Y", "CGI", "Other")

ggplot(data=data_HvP_all_DEstat[data_HvP_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=HvP_DMR_freq  , col=Tissue, fill=Tissue, group = interaction(DE_status_HvP , Tissue), lty=DE_status_HvP ))+
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


ggplot(data=data_HvP_all_DEstat[data_HvP_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=non_CpG_Fst  , col=Tissue, fill=Tissue, group = interaction(DE_status_HvP , Tissue), lty=DE_status_HvP ))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))









#Transgressive DEstat ####

data_Transgressive_all_DEstat <- data.frame()

for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
  load(paste(tissue, "Transgressive_DE_status_CGI_all.rda", sep="_"))
  all.res$Tissue <- tissue
  data_Transgressive_all_DEstat <- rbind(data_Transgressive_all_DEstat, all.res)
}

#HvP - DMR_freq ####
data_Transgressive_all_DEstat<-data_Transgressive_all_DEstat[!is.na(data_Transgressive_all_DEstat$DE_status_Transgressive ),]

#Promoter
CGI_prom="N"
for ( tissue in c("Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t<-t.test(data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>30 & data_Transgressive_all_DEstat$Segment_number<51 & data_Transgressive_all_DEstat$DE_status_Transgressive == "DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq, 
            data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>30 & data_Transgressive_all_DEstat$Segment_number<51 & data_Transgressive_all_DEstat$DE_status_Transgressive == "non-DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq,paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=30))
}



#Gene body
CGI_prom="N"
for ( tissue in c("Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>50 & data_Transgressive_all_DEstat$Segment_number<150 & data_Transgressive_all_DEstat$DE_status_Transgressive == "DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq, 
              data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>50 & data_Transgressive_all_DEstat$Segment_number<150 & data_Transgressive_all_DEstat$DE_status_Transgressive == "non-DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvP_DMR_freq,paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=30))
}






#HvC - DMR_freq ####
data_Transgressive_all_DEstat<-data_Transgressive_all_DEstat[!is.na(data_Transgressive_all_DEstat$DE_status_Transgressive ),]

#Promoter
CGI_prom="N"
for ( tissue in c("Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t<-t.test(data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>30 & data_Transgressive_all_DEstat$Segment_number<51 & data_Transgressive_all_DEstat$DE_status_Transgressive == "DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq, 
            data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>30 & data_Transgressive_all_DEstat$Segment_number<51 & data_Transgressive_all_DEstat$DE_status_Transgressive == "non-DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq,paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=30))
}



#Gene body
CGI_prom="N"
for ( tissue in c("Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>50 & data_Transgressive_all_DEstat$Segment_number<150 & data_Transgressive_all_DEstat$DE_status_Transgressive == "DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq, 
              data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>50 & data_Transgressive_all_DEstat$Segment_number<150 & data_Transgressive_all_DEstat$DE_status_Transgressive == "non-DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$HvC_DMR_freq,paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=30))
}



#non-CpG Fst ####
#Promoter
CGI_prom="N"
for ( tissue in c( "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>30 & data_Transgressive_all_DEstat$Segment_number<51 & data_Transgressive_all_DEstat$DE_status_Transgressive == "DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
              data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>30 & data_Transgressive_all_DEstat$Segment_number<51 & data_Transgressive_all_DEstat$DE_status_Transgressive == "non-DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=30))
}


#Gene body
CGI_prom="N"
for ( tissue in c( "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  t <- t.test(data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>50 & data_Transgressive_all_DEstat$Segment_number<150 & data_Transgressive_all_DEstat$DE_status_Transgressive == "DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, 
              data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Segment_number>50 & data_Transgressive_all_DEstat$Segment_number<150 & data_Transgressive_all_DEstat$DE_status_Transgressive == "non-DE" & data_Transgressive_all_DEstat$Tissue == tissue & data_Transgressive_all_DEstat$CGI_prom == CGI_prom, ]$non_CpG_Fst, paired=T)
  #print(t$estimate)
  print(t$p.value)
  print(p.adjust(t$p.value, method="bonferroni", n=30))
}





#
data_Transgressive_all_DEstat$Promoter_type <- ifelse(data_Transgressive_all_DEstat$CGI_prom == "Y", "CGI", "Other")

ggplot(data=data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Promoter_type == "Other",], aes(x=Segment_number, y=HvC_DMR_freq, col = Tissue, group = interaction(DE_status_Transgressive , Tissue), lty=DE_status_Transgressive ))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  #geom_smooth(data=data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Promoter_type == "Other",], aes(x=Segment_number, y=HvC_DMR_freq, fill=Tissue), method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_smooth(data=data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Promoter_type == "Other",], aes(x=Segment_number, y=HvP_DMR_freq, fill=Tissue), method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab("DMR frequency")+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  #scale_color_manual(name="DMR comparison", values = c("HYBvCOL"="red", "HYBvPIE"="blue"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))


ggplot(data=data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Promoter_type == "Other" & data_Transgressive_all_DEstat$Tissue != "Testis" & data_Transgressive_all_DEstat$Tissue != "Brain" ,], aes(x=Segment_number, y=non_CpG_Fst  , col=Tissue, fill=Tissue, group = interaction(DE_status_Transgressive , Tissue), lty=DE_status_Transgressive ))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("red", "orange", "olivedrab3")) +
  scale_fill_manual(name="Tissue", values = c( "red", "orange", "olivedrab3")) +
  geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  ylab(expression("non-CpG" ~ italic("F"["ST"])))+
  scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=23, colour="black"), axis.title=element_text(size=26), legend.text=element_text(size=24),  legend.title=element_text(size=26))

ggplot(data=data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=non_CpG_Fst, col=Tissue, fill=Tissue, lty=DE_status_Transgressive))+
  theme_classic()+
  #geom_ribbon(alpha=0.3, colour=NA)+
  #geom_line() +
  scale_color_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  scale_fill_manual(name="Tissue", values = c("honeydew4", "red", "orange", "olivedrab3", "purple")) +
  #geom_smooth(method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  geom_smooth(data=data_Transgressive_all_DEstat[data_Transgressive_all_DEstat$Promoter_type == "CGI",], aes(x=Segment_number, y=non_CpG_Fst), method = "loess", span = 0.15, method.args = list(degree=1), alpha=0.2)+
  
  geom_vline(xintercept=50, alpha=0.3)+
  geom_vline(xintercept=150, alpha=0.3)+
  #ylab("DMR frequency")+
  facet_grid(~Tissue)+
  scale_x_continuous(name ="", labels=c("","","","",""))+
  #scale_x_continuous(name ="", labels=c("-5 kb","TSS","Gene Body","TTS","5 kb"))+
  theme(strip.text.x = element_text(size = 20), aspect.ratio=1, legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14), legend.text=element_text(size=14),  legend.title=element_text(size=14))



