library(readr)
library(tidyverse)
library(patchwork)

#read in ----
norm_wide<-read_csv("clean_data/norm_wide_clean.csv")

#read in ----
express_table_full<-read_csv("clean_data/gene_interest_express_GC_all.csv")

express_table_full$gene<-as.factor(express_table_full$gene)
express_table_full$level<-as.factor(express_table_full$level)
express_table_full$Block<-as.factor(express_table_full$Block)
express_table_full$protID<-as.factor(express_table_full$protID)

express_table_full<- express_table_full[which(express_table_full$gene %in% colnames(norm_wide[,55:108])),] #142 unique proteins in 54 differnt GH / cazy combos 

copynum<-as.data.frame(table(express_table_full[express_table_full$levelblock == "lowblock11",]$gene))
copynum<- copynum[which(copynum$Var1 %in% colnames(norm_wide[,55:108])),]

norm_all_prots<- pivot_wider(express_table_full[,-c(2,4,6,7,9:15)], names_from = "prot_gene",values_from = TPM)
norm_all_prots<-norm_all_prots[-which(norm_all_prots$Block == "block19"),]
norm_all_prots<-norm_all_prots[match(norm_wide$levelblock,paste(norm_all_prots$level,norm_all_prots$Block,sep = "")),]
norm_all_prots$MnP<- norm_wide$MnP

## all GH grouped vs MnP ----

GH16_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GH16)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "GH16")+
  xlab(label = "MnP")+
  theme_classic()

GH128_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GH128)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "GH128")+
  xlab(label = "MnP")+
  theme_classic()

GH_combo_MnP_cor<-cor(sqrt(norm_wide$MnP),sqrt(norm_wide[,55:108]))
GH_combo_express_sum<-colSums(norm_wide[,55:108])
GH_combo_count <- copynum$Freq
GH_short_combo<-as.data.frame(t(rbind(GH_combo_MnP_cor,GH_combo_express_sum,GH_combo_count)))
GH_short_combo$names<-row.names(GH_short_combo)

anno_plus_GH<-read_delim("clean_data/GH_annotation_v2.csv",delim = ";")
anno_plus_GH$combo <- if_else(is.na(anno_plus_GH$ecNum),anno_plus_GH$cazy,anno_plus_GH$ecNum)
anno_plus_GH$combo <- paste(anno_plus_GH$ecNum,anno_plus_GH$cazy,sep = "_")
#one double entry bc I annotated some GH manually to EC numbers (V2 of annotation file - whereas are NA_EC in V1) - remove row 11 where 3.2.1.4 is duped
#anno_plus_GH<- anno_plus_GH[-7,]

#problem is now some are called by GH in the GH_short file not the EC nums I assigned 
GH_short_combo[which(GH_short_combo$names == "NA_GH85"),]$names <- "3.2.1.96_GH85"
GH_short_combo[which(GH_short_combo$names == "NA_GH152"),]$names <- "3.2.1.39_GH152"
GH_short_combo[which(GH_short_combo$names == "NA_GH125"),]$names <- "3.2.1.163_GH125"
GH_short_combo[which(GH_short_combo$names == "3.2.1.3_CBM20"),]$names<-"3.2.1.3_GH15_CBM20"

GH_short_combo <- GH_short_combo[match(anno_plus_GH$combo,GH_short_combo$names),] #this will short but also remove some that i didnt think needed to be included
#removed 3.2.1.106, 3.2.1.113, 2.4.1.18, 3.6.4.12, GH18, GH72_CBM43, GH79
GH_short_combo <- cbind(GH_short_combo,anno_plus_GH[,1:7])
GH_short_combo<- GH_short_combo[-7,]

GH_long_combo<- pivot_longer(GH_short_combo,cols = c(1:3))
GH_long_combo_sum<-GH_long_combo[which(GH_long_combo$name == "GH_combo_express_sum"),]
GH_long_combo_sum$names<-factor(GH_long_combo_sum$names,levels = (GH_long_combo_sum$names)[order(GH_long_combo_sum$bonds)])
GH_long_combo_cor<-GH_long_combo[which(GH_long_combo$name == "V1"),]
GH_long_combo_cor$names<-factor(GH_long_combo_cor$names,levels = (GH_long_combo_cor$names)[order(GH_long_combo_cor$bonds)])
GH_long_combo_count<-GH_long_combo[which(GH_long_combo$name == "GH_combo_count"),]
GH_long_combo_count$names<-factor(GH_long_combo_count$names,levels = (GH_long_combo_count$names)[order(GH_long_combo_count$bonds)])


cor<-ggplot(GH_long_combo_cor, aes(name, names, fill= value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "darkblue",mid = "white",high = "darkred",midpoint = 0)+
  theme_classic()


sum<-ggplot(GH_long_combo_sum, aes(name, names, fill= value)) + 
  geom_tile()+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())
  
copy<- ggplot(GH_long_combo_count, aes(name, names)) + 
  geom_text(aes(label=value),size = 5)+
  geom_point(aes(alpha = value),size=6,colour="darkblue")+
  scale_alpha(range = c(0.1,0.8))+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank())

cor+sum+copy
#scale_fill_gradient2(low = "darkblue",mid = "white",high = "darkred",midpoint = 0)
#[-which(GH_long_sum$names %in% c("GH_GH16","GH_GH128","GH_3.2.1.58")),]

## ---- all prots

GH_prots_cor<-as.data.frame(cor(sqrt(norm_all_prots[,-c(1,2,74)]),sqrt(norm_all_prots[,74])))
