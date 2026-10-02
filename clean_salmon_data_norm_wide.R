library(readr)
library(tidyverse)


#read in and clean ----
express_table<-read_csv("clean_data/gene_interest_express_GC.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

express_table<-express_table[-which(express_table$Block == "block19"),]

##wide dfs ----
norm_wide<-express_table[,-c(5,7:11)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

##combine ---- 
norm_wide$AA<-rowSums(norm_wide[,which(colnames(norm_wide) %in% c("OrgN_AAAP","OrgN_ACT","OrgN_LAT","OrgN_YAT","OrgN_AAT"))])
norm_wide$APET<-rowSums(norm_wide[,which(colnames(norm_wide) %in% c("OrgN_POT","OrgN_OPT"))])
norm_wide$GMC_sum<- rowSums(norm_wide[,grepl("GMC",colnames(norm_wide))])
norm_wide$CHIT_sum<- rowSums(norm_wide[,grepl("CHIT",colnames(norm_wide))])
norm_wide$APEP_sum <- rowSums(norm_wide[,grepl("APEP",colnames(norm_wide))])
norm_wide$CPEP_sum <- rowSums(norm_wide[,grepl("CPEP",colnames(norm_wide))])
norm_wide$MCPEP_sum <- rowSums(norm_wide[,grepl("MCPEP",colnames(norm_wide))])
norm_wide$DTPEP_sum <- rowSums(norm_wide[,grepl("DTPEP",colnames(norm_wide))])


##CUE ----

#sum_reads
#reads_wide<-express_table[,-c(4,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)
#reads_wide_sub<- express_table.2[,-c(4,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)
#make ratio
norm_wide$CUE<-norm_wide$GT48/norm_wide$TCA_1.2.4.2
norm_wide$TASE_TPP<-norm_wide$TASE/norm_wide$TPP

write_csv(norm_wide,"clean_data/norm_wide_clean.csv")

