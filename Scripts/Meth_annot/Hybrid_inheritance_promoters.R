#!/usr/bin/env Rscript

library(ggplot2)
library(lmodel2)
library(reshape2)

load("meth_annot_data_wide_CGI_noNA") #From Meth_PCA_and_BCA.R
load("meth_annot_data_wide_Other_noNA") #From Meth_PCA_and_BCA.R



set.seed(1000)
inheritance_data <- data.frame()
firstdatacol=4

for (promoter_type in c("Other", "CGI")){
  for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
    if(promoter_type == "CGI"){
      meth_annot_tissue <- meth_annot_data_wide_CGI_noNA[meth_annot_data_wide_CGI_noNA$Tissue==tissue ,]
      #meth_annot_tissue <- meth_annot_data_wide_CGI[meth_annot_data_wide_CGI$Tissue==tissue ,]
      #meth_annot_tissue <- meth_annot_tissue[ , colSums(is.na(meth_annot_tissue)) == 0]
    }
    if(promoter_type == "Other"){
      meth_annot_tissue <- meth_annot_data_wide_Other_noNA[meth_annot_data_wide_Other_noNA$Tissue==tissue ,]
      #meth_annot_tissue <- meth_annot_data_wide_Other[meth_annot_data_wide_Other$Tissue==tissue ,]
      #meth_annot_tissue <- meth_annot_tissue[ , colSums(is.na(meth_annot_tissue)) == 0]
    }
   
     HYBdata <- meth_annot_tissue[meth_annot_tissue$Species == "HYB",]
     
     if(tissue=="Kidney"){ HYBnum <- 2} else {HYBnum <- 3}
     HYB_sample <- sample(HYBnum, ncol(HYBdata), replace=T)
                  
    
     HYBdata_hyperCOL <- c()
     HYBdata_hyperPIE <- c()
     for (i in 1:length(HYB_sample)){
       indFirst <- HYB_sample[i]
       indSecond <- sample((1:HYBnum)[-indFirst], 1)
        
       HYBdata_hyperCOL <- c(HYBdata_hyperCOL, HYBdata[indFirst, firstdatacol+i-1 ])
       HYBdata_hyperPIE <- c(HYBdata_hyperPIE, HYBdata[indSecond, firstdatacol+i-1 ])
     }
    
    
    H<-colMeans(meth_annot_tissue[meth_annot_tissue$Species == "HYB",firstdatacol:ncol(meth_annot_tissue)]) 
    C<-colMeans(meth_annot_tissue[meth_annot_tissue$Species == "COL",firstdatacol:ncol(meth_annot_tissue)]) 
    P<-colMeans(meth_annot_tissue[meth_annot_tissue$Species == "PIE",firstdatacol:ncol(meth_annot_tissue)])
    
    
    H_v_C <- colMeans(meth_annot_tissue[meth_annot_tissue$Species == "HYB",firstdatacol:ncol(meth_annot_tissue)]) - colMeans(meth_annot_tissue[meth_annot_tissue$Species == "COL",firstdatacol:ncol(meth_annot_tissue)])
    H_v_P <- colMeans(meth_annot_tissue[meth_annot_tissue$Species == "HYB",firstdatacol:ncol(meth_annot_tissue)]) - colMeans(meth_annot_tissue[meth_annot_tissue$Species == "PIE",firstdatacol:ncol(meth_annot_tissue)])
    
    H_v_C_hyper <- HYBdata_hyperCOL - colMeans(meth_annot_tissue[meth_annot_tissue$Species == "COL",firstdatacol:ncol(meth_annot_tissue)])
    H_v_P_hyper <- HYBdata_hyperPIE - colMeans(meth_annot_tissue[meth_annot_tissue$Species == "PIE",firstdatacol:ncol(meth_annot_tissue)])
    
    new_chunk <- cbind(tissue, promoter_type, H_v_C, H_v_P, H, C, P, H_v_C_hyper, H_v_P_hyper, HYBdata_hyperCOL, HYBdata_hyperPIE)
    inheritance_data <- rbind(inheritance_data, new_chunk)
  }
}

inheritance_data$HYBdata_hyperCOL <- as.double(inheritance_data$HYBdata_hyperCOL)
inheritance_data$HYBdata_hyperPIE <- as.double(inheritance_data$HYBdata_hyperPIE)
inheritance_data$H_v_C_hyper <- as.double(inheritance_data$H_v_C_hyper)
inheritance_data$H_v_P_hyper <- as.double(inheritance_data$H_v_P_hyper)
inheritance_data$H_v_C <- as.double(inheritance_data$H_v_C)
inheritance_data$H_v_P <- as.double(inheritance_data$H_v_P)
inheritance_data$H <- as.double(inheritance_data$H)
inheritance_data$P <- as.double(inheritance_data$P)
inheritance_data$C <- as.double(inheritance_data$C)

inheritance_data$Chromosome <- sub('(Chr.*)_.*_.*', '\\1', rownames(inheritance_data))
inheritance_data$Chr_type <- ifelse(inheritance_data$Chromosome == "ChrZ", "Z", "A") 

inheritance_data$C_v_P <- inheritance_data$H_v_C - inheritance_data$H_v_P




inheritance_data$Inpat <- "Conserved"

#Overdominant
inheritance_data$Inpat <- ifelse(sqrt(inheritance_data$H_v_C_hyper^2 + inheritance_data$H_v_P_hyper^2) > 0.1 & inheritance_data$H_v_P_hyper > (pi/8)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper < (8/pi)*inheritance_data$H_v_C_hyper, "Overdominant", inheritance_data$Inpat)
#Underdominant
inheritance_data$Inpat <- ifelse(sqrt(inheritance_data$H_v_C_hyper^2 + inheritance_data$H_v_P_hyper^2) > 0.1 & inheritance_data$H_v_P_hyper < (pi/8)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper > (8/pi)*inheritance_data$H_v_C_hyper, "Underdominant", inheritance_data$Inpat)

#Additive
inheritance_data$Inpat <- ifelse(sqrt(inheritance_data$H_v_C_hyper^2 + inheritance_data$H_v_P_hyper^2) > 0.1 & ((inheritance_data$H_v_P_hyper > (-pi/8)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper < (-8/pi)*inheritance_data$H_v_C_hyper) |  (inheritance_data$H_v_P_hyper < (-pi/8)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper > (-8/pi)*inheritance_data$H_v_C_hyper)), "Additive", inheritance_data$Inpat)

#Pied-dominant
inheritance_data$Inpat <- ifelse(sqrt(inheritance_data$H_v_C_hyper^2 + inheritance_data$H_v_P_hyper^2) > 0.1 & ((inheritance_data$H_v_P_hyper < (-pi/8)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper >(pi/8)*inheritance_data$H_v_C_hyper) |  (inheritance_data$H_v_P_hyper >(-pi/8)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper < (pi/8)*inheritance_data$H_v_C_hyper)), "Pied-dominant", inheritance_data$Inpat)

#Collared-dominant
inheritance_data$Inpat <- ifelse(sqrt(inheritance_data$H_v_C_hyper^2 + inheritance_data$H_v_P_hyper^2) > 0.1 & ((inheritance_data$H_v_P_hyper > (-8/pi)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper >(8/pi)*inheritance_data$H_v_C_hyper) |  (inheritance_data$H_v_P_hyper < (-8/pi)*inheritance_data$H_v_C_hyper & inheritance_data$H_v_P_hyper < (8/pi)*inheritance_data$H_v_C_hyper)), "Collared-dominant", inheritance_data$Inpat)




inheritance_data$Inpat <- factor(inheritance_data$Inpat ,levels = c("Conserved", "Collared-dominant", "Pied-dominant", "Additive", "Overdominant", "Underdominant"))


circles <- data.frame(
  x0 = c(0),
  y0 = c(0),
  r = c(0.1)
)


for (promoter_type in c("Other", "CGI")){
  #
  
  for (tissue in c("Heart", "Kidney", "Liver", "Brain", "Testis")){
    
    
    Areg <- lmodel2(H_v_P_hyper~H_v_C_hyper, data=inheritance_data[inheritance_data$tissue == tissue & inheritance_data$Chr_type == "A" & inheritance_data$promoter_type == promoter_type,])
    Zreg <- lmodel2(H_v_P_hyper~H_v_C_hyper, data=inheritance_data[inheritance_data$tissue == tissue & inheritance_data$Chr_type == "Z" & inheritance_data$promoter_type == promoter_type,])
    test=paste(expression(italic(R)^2 ~ "="), round(Areg$rsquare,2), sep=" ")
    
    ggplot()+
      geom_circle(aes(x0 = x0, y0 = y0, r = r, fill =r ), alpha=0, data = circles)+
      geom_point(data=inheritance_data[inheritance_data$tissue == tissue & inheritance_data$promoter_type == promoter_type & inheritance_data$Inpat != "Conserved",], aes(x=H_v_C_hyper, y=H_v_P_hyper, col=Inpat))+
      geom_point(data=inheritance_data[inheritance_data$tissue == tissue & inheritance_data$promoter_type == promoter_type, ], aes(x=H_v_C_hyper, y=H_v_P_hyper), shape=1, colour = "black", alpha=0.5)+
      
      xlim(-0.5,0.5)+
      ylim(-0.5,0.5)+
      theme_classic()+
      scale_color_viridis_d(begin=1/6, option="inferno") +
      xlab("Hybrid - Collared")+
      ylab("Hybrid - Pied")+
      geom_hline(yintercept=0, alpha = 0.7)+
      geom_vline(xintercept=0, alpha = 0.7)+
      geom_abline(slope=pi/8, intercept=0, lty=2, alpha=0.3)+
      geom_abline(slope=-pi/8, intercept=0, lty=2, alpha=0.3)+
      geom_abline(slope=-8/pi, intercept=0, lty=2, alpha=0.3)+
      geom_abline(slope=8/pi, intercept=0, lty=2, alpha=0.3)+
      
      geom_abline(slope=Areg$regression.results[2,3], intercept=Areg$regression.results[2,2], col="black", size=1)+ #Major axis regression
      geom_abline(slope=Zreg$regression.results[2,3], intercept=Zreg$regression.results[2,2], col="red", size=1)+
      annotate("text",x=-0.3, y= 0.45, size=6, label=bquote(italic(p) ~ "=" ~ .(round(Areg$P.param,3)) ~ ","~ italic(R)^2 ~ "="~.(round(Areg$rsquare,3))), col="black")+
      annotate("text",x=-0.3, y= 0.4, size=6, label=bquote(italic(p) ~ "=" ~ .(round(Zreg$P.param,3)) ~ ","~ italic(R)^2 ~ "="~.(round(Zreg$rsquare,3))), col="red")+
      
      
      #geom_hline(yintercept=0.1, lty=2, alpha=0.3)+
      #geom_hline(yintercept=-0.1, lty=2, alpha=0.3)+
      #geom_vline(xintercept=0.1, lty=2, alpha=0.3)+
      #geom_vline(xintercept=-0.1, lty=2, alpha=0.3)+
    theme(aspect.ratio=1, legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=18, colour="black"), axis.title=element_text(size=21), legend.text=element_text(size=24),  legend.title=element_text(size=26))
    ggsave(filename = paste("Nov_Promoter_HybSplit", promoter_type, tissue,  "Inheritance_pattern.png", sep="_"), plot = last_plot())
    
  }


inheritance_class_data_full_freq<- as.data.frame(table(inheritance_data$tissue, inheritance_data$Inpat,  inheritance_data$promoter_type, inheritance_data$Chr_type))
colnames(inheritance_class_data_full_freq) <-  c("Tissue", "Inpat", "promoter_type", "Chr_type", "Freq")


inheritance_class_data_full_freq$Inpat <- factor(inheritance_class_data_full_freq$Inpat ,levels = c("Conserved", "Collared-dominant", "Pied-dominant", "Additive", "Overdominant", "Underdominant"))


#Proportion conserved per tissue and promoter type
prop.table(table(inheritance_data$Inpat, inheritance_data$promoter_type,  inheritance_data$tissue),c(2,3))
#


prop.table(table(inheritance_data$Inpat, inheritance_data$promoter_type,  inheritance_data$tissue, inheritance_data$Chr_type),c(2,3,4))
inheritance_class_data_full_prop <- as.data.frame(prop.table(table(inheritance_data$Inpat, inheritance_data$promoter_type,  inheritance_data$tissue, inheritance_data$Chr_type),c(2,3,4)))
colnames(inheritance_class_data_full_prop) <-  c("Inpat", "promoter_type", "Tissue", "Chr_type", "Freq")

ggplot(inheritance_class_data_full_prop[inheritance_class_data_full_prop$Inpat != "Conserved",], aes(fill=Inpat, y=Freq, x=Tissue)) + 
  theme_classic()+
  ylab("Proportion")+
  scale_fill_viridis(name="Inheritance pattern", discrete = T, begin=1/6, option="inferno") +
  geom_bar(position = "stack", stat="identity", colour="black")+
  facet_wrap(~promoter_type+Chr_type, nrow=1)+
  theme(strip.text = element_text(size = 20), panel.border = element_rect(colour = "black", fill=NA, size=0.5), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=16))

}


#Statistics ####

df <- data.frame()
for(tissue in unique(inheritance_data$tissue)){
  
  inTiss <- inheritance_data[inheritance_data$tissue == tissue,]
  
  #Does the distribution of inheritance patterns differ between promoter types?
  fp1 <- fisher.test(inTiss$Inpat, inTiss$promoter_type, simulate.p.value = T, B = 1e6)
  xsq <- chisq.test(inTiss$Inpat, inTiss$promoter_type)
  xsq$residuals
  
  #Does the distribution of inheritance patterns differ between promoter types when excluding the conservative class?
  inTiss <- inheritance_data_noCons[inheritance_data_noCons$tissue == tissue,]

  fp2 <- fisher.test(inTiss$Inpat, inTiss$promoter_type, simulate.p.value = T, B = 1e6)
  xsq <- chisq.test(inTiss$Inpat, inTiss$promoter_type)
  print(tissue)
  print(xsq$residuals)
  
  df <- rbind(df, cbind(tissue, pval_w_cons=fp1$p.value, pval_n_cons=fp2$p.value))
  
}


df[df$pval < 0.05,]


fisher.test(inheritance_data$Inpat, inheritance_data$Chr_type, simulate.p.value = T, B = 1e6)
xsq <- chisq.test(inheritance_data$Inpat, inheritance_data$Chr_type)
xsq$residuals

for(tissue in c("Brain", "Heart", "Kidney", "Liver", "Testis")){
  print(tissue)
  print(mean(abs(inheritance_data[inheritance_data$tissue == tissue & inheritance_data$Chr_type == "Z",]$C_v_P))
        -mean(abs(inheritance_data[inheritance_data$tissue == tissue & inheritance_data$Chr_type == "A",]$C_v_P)))
  print(wilcox.test(abs(inheritance_data[inheritance_data$tissue == tissue & inheritance_data$Chr_type == "Z",]$C_v_P),
                    abs(inheritance_data[inheritance_data$tissue == tissue & inheritance_data$Chr_type == "A",]$C_v_P)))
  
}




