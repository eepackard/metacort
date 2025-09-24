library(readr)
library(tidyverse)

#read in ----
express_table_full<-read_csv("clean_data/gene_interest_express_GC_all.csv")

express_table_full$gene<-as.factor(express_table_full$gene)
express_table_full$level<-as.factor(express_table_full$level)
express_table_full$Block<-as.factor(express_table_full$Block)
express_table_full$protID<-as.factor(express_table_full$protID)

ggplot(express_table_full[which(express_table_full$gene == "MnP"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  geom_dotplot(binaxis = 'y',stackdir = 'center',position = position_dodge(0.8),binwidth = 0.1)+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "MnP"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
 theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "MnP" & express_table_full$sigP > 0.5),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# something wrong with the sigP

ggplot(express_table_full[which(express_table_full$gene == "LACC" & express_table_full$sigP > 0.8),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #2079975 which has low expression does not have sigP

ggplot(express_table_full[which(express_table_full$gene == "LACC"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "CAT"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #neither have sigP but no surprise as they are generally intracellular

ggplot(express_table_full[which(express_table_full$gene == "SOD"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#no signalP 

ggplot(express_table_full[which(express_table_full$gene == "PROT"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "PROT" & express_table_full$sigP > 0.8),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() # 33 of 49 protiens have signP > 0.5, and 29 of 49 > 0.8 

ggplot(express_table_full[which(express_table_full$gene == "APEP"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "APEP" & express_table_full$sigP > 0.5),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #only 1825072 has sigP (out of 18) not highly expressed

ggplot(express_table_full[which(express_table_full$gene == "GLY"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "GLY" & express_table_full$sigP > 0.8),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #2 of 5 have no sigP and they are not highly expressed (only 1594518 (CRO2), 1779395 (CRO1), and 1984128 (CRO3/4/5) have sigP > 0.8)

ggplot(express_table_full[which(grepl("GMC",express_table_full$gene)),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# only one has sigP 1861140 (GMC 1.1.3.13 = alcohol oxidase)

ggplot(express_table_full[which(grepl("GMC",express_table_full$gene)),],aes(x=level,y=norm_reads,fill = gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "GMC_1.1.3.7"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AAO

ggplot(express_table_full[which(express_table_full$gene == "GMC_1.1.3.13"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()#AoX

ggplot(express_table_full[which(express_table_full$gene == "MnP"),],aes(x=level,y=BTub))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

norm_wide<-express_table[,-c(4,5)] |>  pivot_wider(names_from = "gene",values_from = norm_reads)


