library(readr)
library(tidyverse)

#read in ----
express_table_full<-read_csv("clean_data/gene_interest_express_GC_all.csv")

express_table_full$gene<-as.factor(express_table_full$gene)
express_table_full$level<-as.factor(express_table_full$level)
express_table_full$Block<-as.factor(express_table_full$Block)
express_table_full$protID<-as.factor(express_table_full$protID)

ggplot(express_table_full[which(express_table_full$gene == "MnP"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  geom_dotplot(binaxis = 'y',stackdir = 'center',position = position_dodge(0.8),binwidth = 0.1)+
  theme_classic()

tiff("figures/MnP_expression_per_copy.tiff")
ggplot(express_table_full[which(express_table_full$gene == "MnP"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
 theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "MnP" & express_table_full$sigP > 0.5),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# something wrong with the sigP

ggplot(express_table_full[which(express_table_full$gene == "KGD"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#only 3 of 8 are expressed

ggplot(express_table_full[which(express_table_full$gene == "ALAS"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#single copy

ggplot(express_table_full[which(express_table_full$gene == "LACC" & express_table_full$sigP > 0.8),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #2079975 which has low expression does not have sigP

ggplot(express_table_full[which(express_table_full$gene == "LACC"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "CAT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #neither have sigP but no surprise as they are generally intracellular - only 1 is expressed

ggplot(express_table_full[which(express_table_full$gene == "SOD"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#no signalP 

ggplot(express_table_full[which(express_table_full$gene == "PROT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "PROT" & express_table_full$sigP > 0.8),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() # 33 of 49 protiens have signP > 0.5, and 29 of 49 > 0.8 - driven largely by one protein 1600311

ggplot(express_table_full[which(express_table_full$gene == "APEP"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "APEP" & express_table_full$sigP > 0.5),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #only 1825072 has sigP (out of 18) not highly expressed - this is from 3.4.11.10 

ggplot(express_table_full[which(express_table_full$gene == "GLY"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "GLY" & express_table_full$sigP > 0.8),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #2 of 5 have no sigP and they are not highly expressed (only 1594518 (CRO2), 1779395 (CRO1), and 1984128 (CRO3/4/5) have sigP > 0.8)

ggplot(express_table_full[which(grepl("GMC",express_table_full$gene)),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# only one has sigP 1861140 (GMC 1.1.3.13 = alcohol oxidase)

ggplot(express_table_full[which(grepl("GMC",express_table_full$gene)),],aes(x=level,y=TPM,fill = gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "GMC_GDH"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AAO/GDH/PDH

ggplot(express_table_full[which(express_table_full$gene == "GMC_AOx"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AoX

ggplot(express_table_full[which(express_table_full$gene == "GMC_AAO"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AAO

ggplot(express_table_full[which(express_table_full$gene == "CHIT_3.2.1.14" & express_table_full$sigP > 0.8),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full[which(express_table_full$gene == "CHIT_3.2.1.14" ),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full[which(express_table_full$gene == "CHIT_3.2.1.52" & express_table_full$sigP > 0.8),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full[which(express_table_full$gene == "CHIT_3.2.1.52" ),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full[which(express_table_full$gene == "CHIT_3.5.1.41" & express_table_full$sigP > 0.8),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full[which(express_table_full$gene == "CHIT_3.5.1.41" ),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#chitinase

ggplot(express_table_full[which(express_table_full$gene == "OrgN_OPT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#2 that are expressed - the rest is mostly noise? 

ggplot(express_table_full[which(express_table_full$gene == "OrgN_POT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#veryyyy low expression - only 1 gene

ggplot(express_table_full[which(express_table_full$gene == "OrgN_AAAP"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#interestingly it is really only the two that didnt have KOG classification that are expressed - one other lowly expressed 1726313

ggplot(express_table_full[which(express_table_full$gene == "OrgN_ACT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "OrgN_LAT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# all expressed a little but not much relative to others

ggplot(express_table_full[which(express_table_full$gene == "OrgN_YAT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly two expressed

ggplot(express_table_full[which(express_table_full$gene == "OrgN_AAT"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# not nearly as expressed as YAT

ggplot(express_table_full[grepl("N_ammonium",express_table_full$gene),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# 

ggplot(express_table_full[which(express_table_full$gene == "N_Nitrate"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# 

ggplot(express_table_full[which(express_table_full$gene == "N_DHA"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# 

ggplot(express_table_full[which(express_table_full$gene == "GLUC"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly two expressed

ggplot(express_table_full[which(express_table_full$gene == "GLUC"& express_table_full$sigP > 0.9),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# and only 2 are secreted - maybe slightly more in high than low for prot 1858453

ggplot(express_table_full[which(express_table_full$gene == "BETA"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "BETA" & express_table_full$sigP > 0.9),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# some almost all are expressed but only two are secreted - and no obvious difference in expression between high and low 

ggplot(express_table_full[which(express_table_full$gene == "GH_NA_GH16"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly one expressed

ggplot(express_table_full[which(express_table_full$gene == "GH16"& express_table_full$sigP > 0.9),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# and only 8 are secreted - but maybe leaning towards positive correlation for the highly expressed gene? 

ggplot(express_table_full[which(express_table_full$gene == "TASE"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "TASE" & express_table_full$sigP > 0.9),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() 

ggplot(express_table_full[which(express_table_full$gene == "TPP"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "TPP" & express_table_full$sigP > 0.9),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "TPS"),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#not excreted


#tripep from explore
Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

tripept<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.4.14.9"),]
express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")
tripe.express<- express_all[which(express_all$protID %in% tripept$proteinId),]
Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab",delim = "\t", escape_double = FALSE,
                                            trim_ws = TRUE)
tripe.express$sigP <- Coromn1_FilteredModels1_sigp6[match(tripe.express$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob

ggplot(tripe.express,aes(x=level,y=TPM,fill = as.factor(protID)))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #only three are strongly expressed - 1726572 , 1760118, 1985389

ggplot(tripe.express[which(tripe.express$sigP > 0.6),],aes(x=level,y=TPM,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()



