library(readr)
library(tidyverse)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

##interested in 
##MnP = EC 1.11.1.13
##1,3-beta-glucan synthase = 2.4.1.34
##2-oxoglutarate dehydrogenase, E1 subunit = 1.2.4.2
## alcohol axidase (H2O2 production?) = 1.1.3.13
## catalase = 1.11.1.6
##glyoxal oxidase = 1.2.3.15
##super oxide dismutase = 1.15.1.1
#ALAS = 2.3.1.37 
##chitinase = 3.2.1.14

##potential reference
##beta tubulin = prot numbers 1845132 and 1851235


MnP_prots<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.13"),]
GT48<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.34"),]
KGD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.4.2"),]
GMC<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.1.3.13"),]
CAT<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.6"),]
GLY<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.3.15"),]
SOD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.15.1.1"),]
ALAS<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.3.1.37"),]
CHIT<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.2.1.14"),]

GENE_list<-list(MnP_prots,GT48,KGD,GMC,CAT,GLY,SOD,ALAS,CHIT)
names(GENE_list)<-c("MNP","GT48","KGD","GMC","CAT","GLY","SOD","ALAS","CHIT")

#samples 1 A1, 2 D2, 3 B1, 4 D3, 7 A4, 8 D1, 9 D4, 11 B1, 12 A4, 13 C2, 15 A2, 17 B3, 18 B1, 19 D2, 20 D1
high_list<-read_delim("raw_data/high.list",delim = "\t",trim_ws = TRUE,col_names = FALSE)
high_list<-as.list(t(high_list))
high_list_files<-as.list(gsub("quant_salmon","salmon_quant/quant_salmon",high_list))#add path

#samples 1 C1, 2 A4, 3 D4, 4 B2, 7 C1, 8 D3, 9 B2, 11 B4, 12 C2, 13 B3, 15 D1, 17 C4, 18 A4, 19 A2, 20 C3
low_list<-read_delim("raw_data/low.list",delim = "\t",trim_ws = TRUE,col_names = FALSE)
low_list<-as.list(t(low_list))
low_list_files<-as.list(gsub("quant_salmon","salmon_quant/quant_salmon",low_list))#add path

quant_low<-list()
quant_high<-list()

for (i in 1:15) {
quant_low[[i]] <- read_delim(low_list_files[[i]], 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
quant_high[[i]] <- read_delim(high_list_files[[i]], 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
}

#name the first column to level so that I can sub in jgi for high or low of MnP
for (i in 1:15) {
quant_low[[i]]<- separate(quant_low[[i]],"Name",into=c("level","Genome","protID","CEnum"),sep="\\|")
quant_low[[i]]$level<-gsub("jgi","low",quant_low[[i]]$level)
quant_high[[i]]<- separate(quant_high[[i]],"Name",into=c("level","Genome","protID","CEnum"),sep="\\|")
quant_high[[i]]$level<-gsub("jgi","high",quant_high[[i]]$level)
}

quant_full<-list()
quant<-list()

for (i in 1:15) {
quant_full[[i]]<-bind_rows(quant_low[[i]],quant_high[[i]])  
quant[[i]]<-quant_full[[i]][,c(1,3,7,8)] 
names(quant[[i]])<-c(paste("block",[[i]]))
}



group_by(level) %>% 

quant_low_MnP_TPM<-list()
quant_low_MnP_reads<-list()
quant_low_GT48_TPM<-list()
quant_low_GT48_reads<-list()
quant_low_KGD_TPM<-list()
quant_low_KGD_reads<-list()
quant_low_GMC_TPM<-list()
quant_low_GMC_reads<-list()
quant_low_CAT_TPM<-list()
quant_low_CAT_reads<-list()
quant_low_Btub_TPM<-list()
quant_low_Btub_reads<-list()
quant_low_GLY_TPM<-list()
quant_low_GLY_reads<-list()
quant_low_SOD_TPM<-list()
quant_low_SOD_reads<-list()
quant_low_ALAS_TPM<-list()
quant_low_ALAS_reads<-list()
quant_low_CHIT_TPM<-list()
quant_low_CHIT_reads<-list()

quant_high_MnP_TPM<-list()
quant_high_MnP_reads<-list()
quant_high_GT48_TPM<-list()
quant_high_GT48_reads<-list()
quant_high_KGD_TPM<-list()
quant_high_KGD_reads<-list()
quant_high_GMC_TPM<-list()
quant_high_GMC_reads<-list()
quant_high_CAT_TPM<-list()
quant_high_CAT_reads<-list()
quant_high_Btub_TPM<-list()
quant_high_Btub_reads<-list()
quant_high_GLY_TPM<-list()
quant_high_GLY_reads<-list()
quant_high_SOD_TPM<-list()
quant_high_SOD_reads<-list()
quant_high_ALAS_TPM<-list()
quant_high_ALAS_reads<-list()
quant_high_CHIT_TPM<-list()
quant_high_CHIT_reads<-list()

quant_high_MnP_norm<-list()
quant_high_GT48_norm<-list()
quant_high_KGD_norm<-list()
quant_high_GMC_norm<-list()
quant_high_CAT_norm<-list()
quant_high_GLY_norm<-list()
quant_high_SOD_norm<-list()
quant_high_ALAS_norm<-list()
quant_high_CHIT_norm<-list()

quant_low_MnP_norm<-list()
quant_low_GT48_norm<-list()
quant_low_KGD_norm<-list()
quant_low_GMC_norm<-list()
quant_low_CAT_norm<-list()
quant_low_GLY_norm<-list()
quant_low_SOD_norm<-list()
quant_low_ALAS_norm<-list()
quant_low_CHIT_norm<-list()

for (i in 1:15) {
quant_low_MnP_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% MnP_prots$proteinId),7])
quant_low_MnP_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% MnP_prots$proteinId),8])
quant_low_GT48_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% GT48$proteinId),7])
quant_low_GT48_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% GT48$proteinId),8])
quant_low_KGD_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% KGD$proteinId),7])
quant_low_KGD_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% KGD$proteinId),8])
quant_low_GMC_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% GMC$proteinId),7])
quant_low_GMC_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% GMC$proteinId),8])
quant_low_CAT_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% CAT$proteinId),7])
quant_low_CAT_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% CAT$proteinId),8])
quant_low_SOD_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% SOD$proteinId),7])
quant_low_SOD_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% SOD$proteinId),8])
quant_low_ALAS_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% ALAS$proteinId),7])
quant_low_ALAS_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% ALAS$proteinId),8])
quant_low_CHIT_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% CHIT$proteinId),7])
quant_low_CHIT_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% CHIT$proteinId),8])
quant_low_GLY_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% GLY$proteinId),7])
quant_low_GLY_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% GLY$proteinId),8])
quant_low_Btub_TPM[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% c("1845132","1851235")),7])
quant_low_Btub_reads[[i]]<-sum(quant_low[[i]][which(quant_low[[i]]$protID %in% c("1845132","1851235")),8])


quant_high_MnP_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% MnP_prots$proteinId),7])
quant_high_MnP_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% MnP_prots$proteinId),8])
quant_high_GT48_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% GT48$proteinId),7])
quant_high_GT48_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% GT48$proteinId),8])
quant_high_KGD_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% KGD$proteinId),7])
quant_high_KGD_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% KGD$proteinId),8])
quant_high_GMC_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% GMC$proteinId),7])
quant_high_GMC_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% GMC$proteinId),8])
quant_high_CAT_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% CAT$proteinId),7])
quant_high_CAT_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% CAT$proteinId),8])
quant_high_SOD_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% SOD$proteinId),7])
quant_high_SOD_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% SOD$proteinId),8])
quant_high_ALAS_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% ALAS$proteinId),7])
quant_high_ALAS_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% ALAS$proteinId),8])
quant_high_CHIT_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% CHIT$proteinId),7])
quant_high_CHIT_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% CHIT$proteinId),8])
quant_high_GLY_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% GLY$proteinId),7])
quant_high_GLY_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% GLY$proteinId),8])
quant_high_Btub_TPM[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% c("1845132","1851235")),7])
quant_high_Btub_reads[[i]]<-sum(quant_high[[i]][which(quant_high[[i]]$protID %in% c("1845132","1851235")),8])

quant_high_MnP_norm[[i]]<-quant_high_MnP_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_GT48_norm[[i]]<-quant_high_GT48_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_KGD_norm[[i]]<-quant_high_KGD_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_GMC_norm[[i]]<-quant_high_GMC_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_SOD_norm[[i]]<-quant_high_SOD_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_ALAS_norm[[i]]<-quant_high_ALAS_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_CHIT_norm[[i]]<-quant_high_CHIT_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_CAT_norm[[i]]<-quant_high_CAT_reads[[i]]/quant_high_Btub_reads[[i]]
quant_high_GLY_norm[[i]]<-quant_high_GLY_reads[[i]]/quant_high_Btub_reads[[i]]

quant_low_MnP_norm[[i]]<-quant_low_MnP_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_GT48_norm[[i]]<-quant_low_GT48_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_KGD_norm[[i]]<-quant_low_KGD_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_GMC_norm[[i]]<-quant_low_GMC_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_SOD_norm[[i]]<-quant_low_SOD_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_ALAS_norm[[i]]<-quant_low_ALAS_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_CHIT_norm[[i]]<-quant_low_CHIT_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_CAT_norm[[i]]<-quant_low_CAT_reads[[i]]/quant_low_Btub_reads[[i]]
quant_low_GLY_norm[[i]]<-quant_low_GLY_reads[[i]]/quant_low_Btub_reads[[i]]
}


#make df
gene<-c(rep(c(rep("MnP",2),rep("GT48",2),rep("KGD",2),rep("GMC",2),rep("SOD",2),rep("ALAS",2),rep("CHIT",2),rep("CAT",2),rep("GLY",2),rep("Btub",2)),15))
level<-c(rep(c("low","high"),10*15))
TPM<-list()
reads<-list()
norm<-list()

for (i in 1:15){
  TPM[[i]]<-c(quant_low_MnP_TPM[[i]],quant_high_MnP_TPM[[i]],quant_low_GT48_TPM[[i]],quant_high_GT48_TPM[[i]],quant_low_KGD_TPM[[i]],quant_high_KGD_TPM[[i]],quant_low_GMC_TPM[[i]],quant_high_GMC_TPM[[i]],quant_low_SOD_TPM[[i]],quant_high_SOD_TPM[[i]],quant_low_ALAS_TPM[[i]],quant_high_ALAS_TPM[[i]],quant_low_CHIT_TPM[[i]],quant_high_CHIT_TPM[[i]],quant_low_CAT_TPM[[i]],quant_high_CAT_TPM[[i]],quant_low_GLY_TPM[[i]],quant_high_GLY_TPM[[i]],quant_low_Btub_TPM[[i]],quant_high_Btub_TPM[[i]])
  reads[[i]]<-c(quant_low_MnP_reads[[i]],quant_high_MnP_reads[[i]],quant_low_GT48_reads[[i]],quant_high_GT48_reads[[i]],quant_low_KGD_reads[[i]],quant_high_KGD_reads[[i]],quant_low_GMC_reads[[i]],quant_high_GMC_reads[[i]],quant_low_SOD_reads[[i]],quant_high_SOD_reads[[i]],quant_low_ALAS_reads[[i]],quant_high_ALAS_reads[[i]],quant_low_CHIT_reads[[i]],quant_high_CHIT_reads[[i]],quant_low_CAT_reads[[i]],quant_high_CAT_reads[[i]],quant_low_GLY_reads[[i]],quant_high_GLY_reads[[i]],quant_low_Btub_reads[[i]],quant_high_Btub_reads[[i]])
  norm[[i]]<-c(quant_low_MnP_norm[[i]],quant_high_MnP_norm[[i]],quant_low_GT48_norm[[i]],quant_high_GT48_norm[[i]],quant_low_KGD_norm[[i]],quant_high_KGD_norm[[i]],quant_low_GMC_norm[[i]],quant_high_GMC_norm[[i]],quant_low_SOD_norm[[i]],quant_high_SOD_norm[[i]],quant_low_ALAS_norm[[i]],quant_high_ALAS_norm[[i]],quant_low_CHIT_norm[[i]],quant_high_CHIT_norm[[i]],quant_low_CAT_norm[[i]],quant_high_CAT_norm[[i]],quant_low_GLY_norm[[i]],quant_high_GLY_norm[[i]],NA,NA)
}

express_table<-as.data.frame(cbind(gene,level, c(TPM[[1]],TPM[[2]],TPM[[3]],TPM[[4]],TPM[[5]],TPM[[6]],TPM[[7]],TPM[[8]],TPM[[9]],TPM[[10]],TPM[[11]],TPM[[12]],TPM[[13]],TPM[[14]],TPM[[15]]),
                                              c(reads[[1]],reads[[2]],reads[[3]],reads[[4]],reads[[5]],reads[[6]],reads[[7]],reads[[8]],reads[[9]],reads[[10]],reads[[11]],reads[[12]],reads[[13]],reads[[14]],reads[[15]]),
                                              c(norm[[1]],norm[[2]],norm[[3]],norm[[4]],norm[[5]],norm[[6]],norm[[7]],norm[[8]],norm[[9]],norm[[10]],norm[[11]],norm[[12]],norm[[13]],norm[[14]],norm[[15]])))

express_table$rep<-c(rep(1:15,each=20))
colnames(express_table)<-c("gene","level","TPM","reads","norm","rep")
express_table$TPM<-as.numeric(express_table$TPM)
express_table$reads<-as.numeric(express_table$reads)
express_table$norm<-as.numeric(express_table$norm)
express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$rep<-as.factor(express_table$rep)

express_table_no7low<-express_table[-which(express_table$rep == "7" & express_table$level == "low"),]

ggplot(express_table[-which(express_table$gene == "Btub"),])+
  geom_boxplot(aes(x=gene,y=norm,colour = level))+
  theme_classic()

ggplot(express_table[-which(express_table$gene == "Btub"),])+
  geom_boxplot(aes(x=gene,y=TPM,colour = level))+
  theme_classic()

ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=reads,colour = level))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "KGD"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()


ggplot(express_table[which(express_table$gene == "GT48"),])+#& express_table$rep < 15
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "SOD"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "GMC"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "GLY"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "SOD"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "ALAS"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "CHIT"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

##there is some rep/blocks where the difference is much stronger - lets limit to those where difference in normalized difference is greater than 2 
## reps - 1,2,4,8,9,12,13,14

express_table.2<-express_table[which(express_table$rep %in% c(1,2,4,8,9,12,13,14)),]

ggplot(express_table.2[-which(express_table.2$gene == "Btub"),])+
  geom_boxplot(aes(x=gene,y=norm,colour = level))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "GMC"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "GLY"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "GT48"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "ALAS"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "SOD"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "CAT"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "CHIT"),])+
  geom_boxplot(aes(x=level,y=norm,fill = level))+
  geom_line(aes(group = rep,x=level,y=norm))+
  geom_point(aes(fill = level,group = rep,x=level,y=norm))+
  theme_classic()
