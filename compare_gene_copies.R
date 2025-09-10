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

ggplot(express_table_full[which(express_table_full$gene == "LACC"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "PROT"),],aes(x=level,y=norm_reads,fill = protID))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()


norm_wide<-express_table[,-c(4,5)] |>  pivot_wider(names_from = "gene",values_from = norm_reads)


