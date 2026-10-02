library(readr)
library(tidyverse)

#read in ----
express_table_full<-read_csv("clean_data/gene_interest_express_spruce_GC_all.csv")

express_table_full$gene<-as.factor(express_table_full$gene)
express_table_full$level<-as.factor(express_table_full$level)
express_table_full$Block<-as.factor(express_table_full$Block)
express_table_full$prot_gene<-as.factor(express_table_full$prot_gene)


ggplot(express_table_full[which(express_table_full$gene == "Sweet"),],aes(x=level,y=TPM,fill = prot_gene))+
  geom_boxplot(position = position_dodge(0.8))+
 theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "actin"),],aes(x=level,y=TPM,fill = prot_gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic() #only really three (expressed but one wayyy more)

ggplot(express_table_full[which(express_table_full$gene == "tub"),],aes(x=level,y=TPM,fill = prot_gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full[which(express_table_full$gene == "invert"),],aes(x=level,y=TPM,fill = prot_gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

##select sweets more specifically

norm_wide_full<-express_table_full[,-c(2,4,6)] |>  pivot_wider(names_from = "prot_gene",values_from = TPM)
norm_wide_full_sweet <- norm_wide_full[,which(grepl("Sweet",colnames(norm_wide_full)))]

hist(colSums(norm_wide_full_sweet))
#pick out with > 200

sweet_select <- norm_wide_full_sweet[,which(colSums(norm_wide_full_sweet) > 200)]
sweet_rest <- norm_wide_full_sweet[,which(colSums(norm_wide_full_sweet) < 200)]
write_csv(as.data.frame(colnames(sweet_rest)),"clean_data/low_sweets_spruce.csv")

express_table_full_v2 <- express_table_full[-which(express_table_full$prot_gene %in% colnames(sweet_rest)),]
express_table_full_v3 <- express_table_full[-which(express_table_full$prot_gene %in% colnames(sweet_select)),]

ggplot(express_table_full_v2[which(express_table_full_v2$gene == "Sweet"),],aes(x=level,y=TPM,fill = prot_gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()

ggplot(express_table_full_v3[which(express_table_full_v3$gene == "Sweet"),],aes(x=level,y=TPM,fill = prot_gene))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()



