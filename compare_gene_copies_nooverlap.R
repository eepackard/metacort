library(readr)
library(tidyverse)

#read in ----
express_table_full_no<-read_csv("clean_data/gene_interest_express_nooverlap_all.csv")
express_table_full<-read_csv("clean_data/gene_interest_express_GC_all.csv")

express_table_full_no$gene<-as.factor(express_table_full_no$gene)
express_table_full_no$level<-as.factor(express_table_full_no$level)
express_table_full_no$Block<-as.factor(express_table_full_no$Block)
express_table_full_no$protID<-as.factor(express_table_full_no$protID)
express_table_full$gene<-as.factor(express_table_full$gene)
express_table_full$level<-as.factor(express_table_full$level)
express_table_full$Block<-as.factor(express_table_full$Block)
express_table_full$protID<-as.factor(express_table_full$protID)

express_table_full$set <- rep("full",nrow(express_table_full))
express_table_full_no$set <- rep("nooverlapp",nrow(express_table_full_no))

express_table_full_no <- rbind.data.frame(express_table_full_no,express_table_full)

#tiff("figures/MnP_expression_per_copy.tiff")
ggplot(express_table_full_no[which(express_table_full_no$gene == "MnP"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "KGD"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#only 3 of 8 are expressed

ggplot(express_table_full_no[which(express_table_full_no$gene == "ALAS"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#single copy

ggplot(express_table_full_no[which(express_table_full_no$gene == "LACC" & express_table_full_no$sigP > 0.8),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #2079975 which has low expression does not have sigP

ggplot(express_table_full_no[which(express_table_full_no$gene == "LACC"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "CAT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #neither have sigP but no surprise as they are generally intracellular - only 1 is expressed

ggplot(express_table_full_no[which(express_table_full_no$gene == "SOD"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#no signalP 

ggplot(express_table_full_no[which(express_table_full_no$gene == "PROT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "PROT" & express_table_full_no$sigP > 0.8),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() # 33 of 49 protiens have signP > 0.5, and 29 of 49 > 0.8 - driven largely by one protein 1600311

ggplot(express_table_full_no[which(express_table_full_no$gene == "APEP"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "APEP" & express_table_full_no$sigP > 0.5),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #only 1825072 has sigP (out of 18) not highly expressed - this is from 3.4.11.10 

ggplot(express_table_full_no[which(express_table_full_no$gene == "GLY"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "GLY" & express_table_full_no$sigP > 0.8),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #2 of 5 have no sigP and they are not highly expressed (only 1594518 (CRO2), 1779395 (CRO1), and 1984128 (CRO3/4/5) have sigP > 0.8)

ggplot(express_table_full_no[which(grepl("GMC",express_table_full_no$gene)),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# only one has sigP 1861140 (GMC 1.1.3.13 = alcohol oxidase)

ggplot(express_table_full_no[which(grepl("GMC",express_table_full_no$gene)),],aes(x=set,y=TPM,fill = gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "GMC_GDH"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AAO/GDH/PDH

ggplot(express_table_full_no[which(express_table_full_no$gene == "GMC_AOx"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AoX

ggplot(express_table_full_no[which(express_table_full_no$gene == "GMC_AAO"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AAO

ggplot(express_table_full_no[which(express_table_full_no$gene == "CHIT_3.2.1.14" & express_table_full_no$sigP > 0.8),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full_no[which(express_table_full_no$gene == "CHIT_3.2.1.14" ),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full_no[which(express_table_full_no$gene == "CHIT_3.2.1.52" & express_table_full_no$sigP > 0.8),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full_no[which(express_table_full_no$gene == "CHIT_3.2.1.52" ),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full_no[which(express_table_full_no$gene == "CHIT_3.5.1.41" & express_table_full_no$sigP > 0.8),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full_no[which(express_table_full_no$gene == "CHIT_3.5.1.41" ),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_OPT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#2 that are expressed - the rest is mostly noise? 

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_POT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#veryyyy low expression - only 1 gene

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_AAAP"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#interestingly it is really only the two that didnt have KOG classification that are expressed - one other lowly expressed 1726313

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_ACT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_LAT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# all expressed a little but not much relative to others

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_YAT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly two expressed

ggplot(express_table_full_no[which(express_table_full_no$gene == "OrgN_AAT"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# not nearly as expressed as YAT

ggplot(express_table_full_no[grepl("N_ammonium",express_table_full_no$gene),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# 

ggplot(express_table_full_no[which(express_table_full_no$gene == "N_Nitrate"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# 

ggplot(express_table_full_no[which(express_table_full_no$gene == "N_DHA"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# 

ggplot(express_table_full_no[which(express_table_full_no$gene == "GLUC"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly two expressed

ggplot(express_table_full_no[which(express_table_full_no$gene == "GLUC"& express_table_full_no$sigP > 0.9),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# and only 2 are secreted - maybe slightly more in high than low for prot 1858453

ggplot(express_table_full_no[which(express_table_full_no$gene == "BETA"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "BETA" & express_table_full_no$sigP > 0.9),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# some almost all are expressed but only two are secreted - and no obvious difference in expression between high and low 

ggplot(express_table_full_no[which(express_table_full_no$gene == "NA_GH16"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly one expressed

ggplot(express_table_full_no[which(express_table_full_no$gene == "NA_GH16"& express_table_full_no$sigP > 0.9),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# and only 8 are secreted - but maybe leaning towards positive correlation for the highly expressed gene? 

ggplot(express_table_full_no[which(express_table_full_no$gene == "TASE"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "TASE" & express_table_full_no$sigP > 0.9),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() 

ggplot(express_table_full_no[which(express_table_full_no$gene == "TPP"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "TPP" & express_table_full_no$sigP > 0.9),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_no[which(express_table_full_no$gene == "TPS"),],aes(x=set,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#not excreted

