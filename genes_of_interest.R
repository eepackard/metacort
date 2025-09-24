library(readr)
library(tidyverse)

##interested in 
##MnP = EC 1.11.1.13
##1,3-beta-glucan synthase = 2.4.1.34
##2-oxoglutarate dehydrogenase, E1 subunit = 1.2.4.2
## alcohol oxidase (H2O2 production?) = 1.1.3.13
## catalase = 1.11.1.6
##glyoxal oxidase = 1.2.3.15
##super oxide dismutase = 1.15.1.1
#ALAS = 2.3.1.37 
##chitinase = 3.2.1.14

##potential reference
##beta tubulin = prot numbers 1845132 and 1851235
express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)


A1_prot_Coromn <- read_delim("raw_data/A1_prot_Coromn.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                                `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                                Location = col_skip(), Scaffold = col_skip(), 
                                                                                                                Start = col_skip(), End = col_skip(), 
                                                                                                                Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)

GMC<-read_csv("clean_data/GMC_clean.csv")
OrgN<-read_csv("clean_data/N_trans_clean.csv")
REF<-read_csv("clean_data/ref_clean.csv")


MnP<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.13"),]
GT48<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.34"),]
KGD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.4.2"),]
CAT<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.6"),]
GLY<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.3.15"),]
SOD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.15.1.1"),]
ALAS<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.3.1.37"),]
CHIT<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.2.1.14"),]#endochitinase
NAG<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.2.1.52"),]#exochitinase
CHSN<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.16"),]
APEP<-Coromn1_FilteredModels1_ec[which(grepl("3.4.11.",Coromn1_FilteredModels1_ec$ecNum)),]

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}
GENES<-AppendMe(c("MnP","GT48","KGD","CAT","GLY","SOD","ALAS","CHIT","NAG","CHSN","APEP"))

###---- based on EC 
express_interest<-express_all[which(express_all$protID %in% GENES$proteinId),] #find only the genes of interest

###---- proteases
express_prot_all<-express_all[which(express_all$protID %in% A1_prot_Coromn$`Protein Id`),]
plot(as.factor(express_prot_all$protID),express_prot_all$TPM) ## mostly just one main protein - 1600311
express_prot_by_ID<- express_prot_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_prot<- express_prot_all |>  group_by(level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_prot$gene <- rep("PROT",nrow(express_prot))
express_prot<-express_prot[,c(1,2,5,3,4)]
express_prot_all$gene <- rep("PROT",nrow(express_prot_all))

###----laccases
express_lac_all<- express_all[which(express_all$protID %in% c("2079975","1581758","1821865","1205939","1668803")),]
plot(as.factor(express_lac_all$protID),express_lac_all$TPM)
express_lac_by_ID<- express_lac_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_lac<- express_lac_all |>  group_by(level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_lac$gene <- rep("LACC",nrow(express_lac))
express_lac<-express_lac[,c(1,2,5,3,4)]
express_lac_all$gene <- rep("LACC",nrow(express_lac_all))

### ---- GMC
express_GMC_all<- express_all[which(express_all$protID %in% GMC$proteinId),]
express_GMC_all$ecNum<- GMC[match(express_GMC_all$protID,GMC$proteinId),]$ecNum
plot(as.factor(express_GMC_all$protID),express_GMC_all$TPM)
express_GMC_by_ID<- express_GMC_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_GMC<- express_GMC_all |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_GMC$gene <- paste("GMC_",express_GMC$ecNum,sep = "")
express_GMC<-express_GMC[,c(1,2,6,4,5)]
express_GMC_all$gene <- paste("GMC_",express_GMC_all$ecNum,sep = "")
express_GMC_all <- express_GMC_all[,-6]

### ---- GMC
express_OrgN_all<- express_all[which(express_all$protID %in% OrgN$proteinId),]
express_OrgN_all$ecNum<- OrgN[match(express_OrgN_all$protID,OrgN$proteinId),]$gene
plot(as.factor(express_OrgN_all$protID),express_OrgN_all$TPM)
express_OrgN_by_ID<- express_OrgN_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_OrgN<- express_OrgN_all |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_OrgN$gene <- paste("OrgN_",express_OrgN$ecNum,sep = "")
express_OrgN<-express_OrgN[,c(1,2,6,4,5)]
express_OrgN_all$gene <- paste("OrgN_",express_OrgN_all$ecNum,sep = "")
express_OrgN_all <- express_OrgN_all[,-6]


###----reference genes for normalization
express_REF<- express_all[which(express_all$protID %in% REF$protID),] #
express_REF$gene <- REF[match(express_REF$protID,REF$protID),]$gene

### ---- match and combine

express_interest$gene<-GENES[match(express_interest$protID,GENES$proteinId),]$gene#add a column that matches the protienIDs to the gene names

express_interest_full <- rbind(express_interest,express_GMC_all,express_prot_all,express_lac_all,express_OrgN_all)
express_interest_full$sigP <- Coromn1_FilteredModels1_sigp6[match(express_interest_full$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_interest_full$prot_gene <-paste(express_interest_full$protID,express_interest_full$gene,sep = "_")
express_interest <- express_interest |>  group_by(level,Block,gene) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))  #add together reads from several gene copies
express_interest <- rbind(express_interest,express_lac,express_prot,express_GMC,express_OrgN)

express_interest$levelblock<-paste(express_interest$level,express_interest$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_interest_full$levelblock<-paste(express_interest_full$level,express_interest_full$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_REF$levelblock<-paste(express_REF$level,express_REF$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping

express_ref_wide<-express_REF[,-c(2,3)] |>  pivot_wider(names_from = "gene",values_from = NumReads)

express_interest<-cbind(express_interest,express_ref_wide[match(express_interest$levelblock,express_ref_wide$levelblock),4:8]) # add column with the number of BTub reads per level and block
express_interest_full<-cbind(express_interest_full,express_ref_wide[match(express_interest_full$levelblock,express_ref_wide$levelblock),4:8])

express_interest$UBC <- express_interest$sum_reads/express_interest$UBC
express_interest$TEF_2 <- express_interest$sum_reads/express_interest$TEF_2
express_interest$TEF_1 <- express_interest$sum_reads/express_interest$TEF_1
express_interest$BTub_1 <- express_interest$sum_reads/express_interest$BTub_1
express_interest$BTub_2 <- express_interest$sum_reads/express_interest$BTub_2

express_interest_full$UBC <- express_interest_full$NumReads/express_interest_full$UBC
express_interest_full$TEF_2 <- express_interest_full$NumReads/express_interest_full$TEF_2
express_interest_full$TEF_1 <- express_interest_full$NumReads/express_interest_full$TEF_1
express_interest_full$BTub_1 <- express_interest_full$NumReads/express_interest_full$BTub_1
express_interest_full$BTub_2 <- express_interest_full$NumReads/express_interest_full$BTub_2


express_interest_wide<-express_interest[,-c(4)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)
express_interest_wide_TPM<-express_interest[,-c(5)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

express_interest_wide_full<-express_interest_full[,-c(2,3,6,7)] |>  pivot_wider(names_from = "prot_gene",values_from = NumReads)
express_interest_wide_TPM_full<-express_interest_full[,-c(2,4,6,7)] |>  pivot_wider(names_from = "prot_gene",values_from = TPM)

express_interest_wide$levelblock<-paste(express_interest_wide$level,express_interest_wide$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_interest_wide_TPM$levelblock<-paste(express_interest_wide_TPM$level,express_interest_wide_TPM$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_interest_wide_full$levelblock<-paste(express_interest_wide_full$level,express_interest_wide_full$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_interest_wide_TPM_full$levelblock<-paste(express_interest_wide_TPM_full$level,express_interest_wide_TPM_full$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_ref_wide$levelblock<-paste(express_ref_wide$level,express_ref_wide$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping

express_interest_wide<-express_interest_wide[match(express_ref_wide$levelblock,express_interest_wide$levelblock),] # add column with the number of BTub reads per level and block
express_interest_wide_TPM<-express_interest_wide_TPM[match(express_ref_wide$levelblock,express_interest_wide_TPM$levelblock),]
express_interest_wide_full<-express_interest_wide_full[match(express_ref_wide$levelblock,express_interest_wide_full$levelblock),]
express_interest_wide_TPM_full<-express_interest_wide_TPM_full[match(express_ref_wide$levelblock,express_interest_wide_TPM_full$levelblock),]


write_csv(express_interest,"clean_data/gene_interest_express_GC.csv") #write data
write_csv(express_interest_full,"clean_data/gene_interest_express_GC_all.csv")
