library(readr)
library(tidyverse)
library(gridExtra)
library(corrgram)


#read in data
express_table<-read_csv("clean_data/gene_interest_express_GC.csv")
express_ref_wide<-read_csv("clean_data/ref_clean_wide.csv")
express_ref_wide_TPM<-read_csv("clean_data/ref_clean_wide_TPM.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

ggplot(express_table)+
  geom_boxplot(aes(x=levelblock,y=sum_TPM))+
  theme_classic()

ggplot(express_table)+
  geom_boxplot(aes(x=level,y=sum_TPM))+
  theme_classic()

express_ref_wide<- express_ref_wide[-which(express_ref_wide$Block == "block19"),]
express_ref_wide_TPM<- express_ref_wide_TPM[-which(express_ref_wide_TPM$Block == "block19"),]

ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=sum_reads,colour = level))+
  theme_classic()

ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=sum_TPM,colour = level))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=UBC,fill = level))+
  geom_line(aes(group = Block,x=level,y=UBC))+
  geom_point(aes(fill = level,group = Block,x=level,y=UBC))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=BTub_1,fill = level))+
  geom_line(aes(group = Block,x=level,y=BTub_1))+
  geom_point(aes(fill = level,group = Block,x=level,y=BTub_1))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=TEF_1,fill = level))+
  geom_line(aes(group = Block,x=level,y=TEF_1))+
  geom_point(aes(fill = level,group = Block,x=level,y=TEF_1))+
  theme_classic()


ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=TEF_1,fill = level))+
  geom_line(aes(group = Block,x=level,y=TEF_1))+
  geom_point(aes(fill = level,group = Block,x=level,y=TEF_1))+
  theme_classic()


corrgram(express_table[,-c(1,2,3,6)],lower.panel = panel.conf)
corrgram(express_ref_wide,lower.panel = panel.conf)
corrgram(express_ref_wide_TPM,lower.panel = panel.conf)


