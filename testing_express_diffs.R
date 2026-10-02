library(readr)
library(tidyverse)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)

##interested in 
##MnP = EC 1.11.1.13
##1,3-beta-glucan synthase = 2.4.1.34
##2-oxoglutarate dehydrogenase, E1 subunit = 1.2.4.2
## alcohol axidase (H2O2 production) = 1.1.3.13
## catalase = 1.11.1.6

##potential reference
##beta tubulin = prot numbers 1845132 and 1851235


MnP_prots<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.13"),]
GT48<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.34"),]
KGD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.4.2"),]
GMC<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.1.3.13"),]
CAT<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.6"),]


quant_coromn_low <- read_delim("raw_data/quant_coromn2.sf", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
quant_coromn_high <- read_delim("raw_data/quant_coromn3.sf", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

quant_coromn_low<- separate(quant_coromn_low,"Name",into=c("jgi","Genome","protID","CEnum"),sep="\\|")
quant_coromn_high<- separate(quant_coromn_high,"Name",into=c("jgi","Genome","protID","CEnum"),sep="\\|")

##now pull out tables for each cooresponding to each

quant_coromn_low_MnP_TPM<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% MnP_prots$proteinId),7])
quant_coromn_low_MnP_reads<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% MnP_prots$proteinId),8])
quant_coromn_low_GT48_TPM<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% GT48$proteinId),7])
quant_coromn_low_GT48_reads<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% GT48$proteinId),8])
quant_coromn_low_KGD_TPM<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% KGD$proteinId),7])
quant_coromn_low_KGD_reads<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% KGD$proteinId),8])
quant_coromn_low_GMC_TPM<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% GMC$proteinId),7])
quant_coromn_low_GMC_reads<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% GMC$proteinId),8])
quant_coromn_low_CAT_TPM<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% CAT$proteinId),7])
quant_coromn_low_CAT_reads<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% CAT$proteinId),8])
quant_coromn_low_Btub_TPM<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% c("1845132","1851235")),7])
quant_coromn_low_Btub_reads<-sum(quant_coromn_low[which(quant_coromn_low$protID %in% c("1845132","1851235")),8])


quant_coromn_high_MnP_TPM<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% MnP_prots$proteinId),7])
quant_coromn_high_MnP_reads<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% MnP_prots$proteinId),8])
quant_coromn_high_GT48_TPM<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% GT48$proteinId),7])
quant_coromn_high_GT48_reads<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% GT48$proteinId),8])
quant_coromn_high_KGD_TPM<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% KGD$proteinId),7])
quant_coromn_high_KGD_reads<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% KGD$proteinId),8])
quant_coromn_high_GMC_TPM<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% GMC$proteinId),7])
quant_coromn_high_GMC_reads<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% GMC$proteinId),8])
quant_coromn_high_CAT_TPM<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% CAT$proteinId),7])
quant_coromn_high_CAT_reads<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% CAT$proteinId),8])
quant_coromn_high_Btub_TPM<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% c("1845132","1851235")),7])
quant_coromn_high_Btub_reads<-sum(quant_coromn_high[which(quant_coromn_high$protID %in% c("1845132","1851235")),8])

quant_coromn_high_MnP_norm<-quant_coromn_high_MnP_reads/quant_coromn_high_Btub_reads
quant_coromn_high_GT48_norm<-quant_coromn_high_GT48_reads/quant_coromn_high_Btub_reads
quant_coromn_high_KGD_norm<-quant_coromn_high_KGD_reads/quant_coromn_high_Btub_reads
quant_coromn_high_GMC_norm<-quant_coromn_high_GMC_reads/quant_coromn_high_Btub_reads
quant_coromn_high_CAT_norm<-quant_coromn_high_CAT_reads/quant_coromn_high_Btub_reads

quant_coromn_low_MnP_norm<-quant_coromn_low_MnP_reads/quant_coromn_low_Btub_reads
quant_coromn_low_GT48_norm<-quant_coromn_low_GT48_reads/quant_coromn_low_Btub_reads
quant_coromn_low_KGD_norm<-quant_coromn_low_KGD_reads/quant_coromn_low_Btub_reads
quant_coromn_low_GMC_norm<-quant_coromn_low_GMC_reads/quant_coromn_low_Btub_reads
quant_coromn_low_CAT_norm<-quant_coromn_low_CAT_reads/quant_coromn_low_Btub_reads

#make df
gene<-c(rep("MnP",2),rep("GT48",2),rep("KGD",2),rep("GMC",2),rep("CAT",2),rep("Btub",2))
level<-c(rep(c("low","high"),6))
TPM<-c(quant_coromn_low_MnP_TPM,quant_coromn_high_MnP_TPM,quant_coromn_low_GT48_TPM,quant_coromn_high_GT48_TPM,quant_coromn_low_KGD_TPM,quant_coromn_high_KGD_TPM,quant_coromn_low_GMC_TPM,quant_coromn_high_GMC_TPM,quant_coromn_low_CAT_TPM,quant_coromn_high_CAT_TPM,quant_coromn_low_Btub_TPM,quant_coromn_high_Btub_TPM)
reads<-c(quant_coromn_low_MnP_reads,quant_coromn_high_MnP_reads,quant_coromn_low_GT48_reads,quant_coromn_high_GT48_reads,quant_coromn_low_KGD_reads,quant_coromn_high_KGD_reads,quant_coromn_low_GMC_reads,quant_coromn_high_GMC_reads,quant_coromn_low_CAT_reads,quant_coromn_high_CAT_reads,quant_coromn_low_Btub_reads,quant_coromn_high_Btub_reads)
norm<-c(quant_coromn_low_MnP_norm,quant_coromn_high_MnP_norm,quant_coromn_low_GT48_norm,quant_coromn_high_GT48_norm,quant_coromn_low_KGD_norm,quant_coromn_high_KGD_norm,quant_coromn_low_GMC_norm,quant_coromn_high_GMC_norm,quant_coromn_low_CAT_norm,quant_coromn_high_CAT_norm,NA,NA)

express_table<-as.data.frame(cbind(gene,level,TPM,reads,norm))
express_table$TPM<-as.numeric(express_table$TPM)
express_table$reads<-as.numeric(express_table$reads)


ggplot(express_table[-which(express_table$gene == "Btub"),])+
  geom_boxplot(aes(x=gene,y=TPM,colour = level))+
  theme_classic()
