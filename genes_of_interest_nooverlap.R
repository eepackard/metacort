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
## Glucan 1,3-beta-glucosidase. = 3.2.1.58 (exo)
## Beta-glucosidase = 3.2.1.21
#ALAS = 2.3.1.37 
##chitinase = 3.2.1.14

##potential reference
##beta tubulin = prot numbers 1845132 and 1851235
express_all<-read_csv("clean_data/expression_salmon_quant_nooverlap_all.csv")

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)


A1_prot_Coromn <- read_delim("raw_data/A1_prot_Coromn.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                                `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                                Location = col_skip(), Scaffold = col_skip(), 
                                                                                                                Start = col_skip(), End = col_skip(), 
                                                                                                                Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

OPT_Coromn <- read_delim("raw_data/OPT_2.A.67.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                        `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                        Location = col_skip(), Scaffold = col_skip(), 
                                                                                                        Start = col_skip(), End = col_skip(), 
                                                                                                        Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

GH_Coromn <- read_delim("raw_data/GH_all_Coromn.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                          `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                          Location = col_skip(), Scaffold = col_skip(), 
                                                                                                          Start = col_skip(), End = col_skip(), 
                                                                                                          Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

GH_Coromn$ecnum <- Coromn1_FilteredModels1_ec[match(GH_Coromn$`Protein Id`,Coromn1_FilteredModels1_ec$proteinId),]$ecNum

AAAP_Coromn <- read_delim("raw_data/AAAP_2.A.18.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                          `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                          Location = col_skip(), Scaffold = col_skip(), 
                                                                                                          Start = col_skip(), End = col_skip(), 
                                                                                                          Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)

transporters <- read_delim("raw_data/transporters.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)

GMC<-read_csv("clean_data/GMC_clean.csv")
OrgN<-read_csv("clean_data/N_trans_clean.csv")
CHIT<-read_csv("clean_data/CHIT_clean.csv")
REF<-read_csv("clean_data/ref_clean.csv")


MnP<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.13"),]## now with updated gene models I will need to add/remove two manually 
MnP[which(MnP$proteinId == "2105191"),]$proteinId <- 1753449
MnP[which(MnP$proteinId == "2105192"),]$proteinId <- 2121138
GT48<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.34"),]
KGD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.4.2"),]
CAT<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.11.1.6"),]
GLY<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.2.3.15"),]## same here
GLY[which(GLY$proteinId == "1725000"),]$proteinId <- 1579362
SOD<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "1.15.1.1"),]
ALAS<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.3.1.37"),]
CHSN<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.16"),]
GLUC<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.2.1.58"),]
BETA<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.2.1.21"),]
APEP<-Coromn1_FilteredModels1_ec[which(grepl("3.4.11.",Coromn1_FilteredModels1_ec$ecNum)),]
TASE<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.2.1.28"),]
TPP<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "3.1.3.12"),]
TPS<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum == "2.4.1.15"),]

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}
GENES<-AppendMe(c("MnP","GT48","KGD","CAT","GLY","SOD","ALAS","CHSN","GLUC","BETA","APEP","TASE","TPP","TPS"))
#96 unique proteins - some repeats because of annotation info
GENES<-GENES[match(unique(GENES$proteinId),GENES$proteinId),]

###---- based on EC 
express_interest<-express_all[which(express_all$protID %in% GENES$proteinId),] #find only the genes of interest
express_interest<-rbind(express_interest,express_all[which(express_all$protID == "1688323"),])#this adds NAgT

# add signal peptide
express_interest$sigP<-Coromn1_FilteredModels1_sigp6[match(express_interest$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob

###---- proteases
express_prot_all<-express_all[which(express_all$protID %in% A1_prot_Coromn$`Protein Id`),]
plot(as.factor(express_prot_all$protID),express_prot_all$TPM) ## mostly just one main protein - 1600311
express_prot_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_prot_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
#express_prot_sec <- express_prot_all[-which(is.na(express_prot_all$sigP)),]
#express_prot_sec <- express_prot_sec[-which(express_prot_sec$sigP < 0.8),]
express_prot_by_ID<- express_prot_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_prot<- express_prot_all |>  group_by(level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_prot$gene <- rep("PROT",nrow(express_prot))
express_prot<-express_prot[,c(1,2,5,3,4)]
express_prot_all$gene <- rep("PROT",nrow(express_prot_all))

###---- OPT
express_OPT_all<-express_all[which(express_all$protID %in% OPT_Coromn$`Protein Id`),]
plot(as.factor(express_OPT_all$protID),express_OPT_all$TPM) 
express_OPT_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_OPT_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_OPT_by_ID<- express_OPT_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_OPT<- express_OPT_all |>  group_by(level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_OPT$gene <- rep("OrgN_OPT",nrow(express_OPT))
express_OPT<-express_OPT[,c(1,2,5,3,4)]
express_OPT_all$gene <- rep("OrgN_OPT",nrow(express_OPT_all))

###---- AAAP
express_AAAP_all<-express_all[which(express_all$protID %in% AAAP_Coromn$`Protein Id`),]
plot(as.factor(express_AAAP_all$protID),express_AAAP_all$TPM) 
express_AAAP_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_AAAP_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_AAAP_by_ID<- express_AAAP_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_AAAP<- express_AAAP_all |>  group_by(level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_AAAP$gene <- rep("OrgN_AAAP",nrow(express_AAAP))
express_AAAP<-express_AAAP[,c(1,2,5,3,4)]
express_AAAP_all$gene <- rep("OrgN_AAAP",nrow(express_AAAP_all))

###----laccases
express_lac_all<- express_all[which(express_all$protID %in% c("2079975","1581758","1821865","1205939","1668803")),]
plot(as.factor(express_lac_all$protID),express_lac_all$TPM)
express_lac_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_lac_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
#express_lac_sec <- express_lac_all[-which(express_lac_all$sigP < 0.8),]
express_lac_by_ID<- express_lac_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_lac<- express_lac_all |>  group_by(level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_lac$gene <- rep("LACC",nrow(express_lac))
express_lac<-express_lac[,c(1,2,5,3,4)]
express_lac_all$gene <- rep("LACC",nrow(express_lac_all))

### ---- GMC
express_GMC_all<- express_all[which(express_all$protID %in% GMC$proteinId),]
express_GMC_all$ecNum<- GMC[match(express_GMC_all$protID,GMC$proteinId),]$based_on_tree
plot(as.factor(express_GMC_all$protID),express_GMC_all$TPM)
express_GMC_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_GMC_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_GMC_by_ID<- express_GMC_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_GMC<- express_GMC_all |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_GMC$gene <- paste("GMC_",express_GMC$ecNum,sep = "")
express_GMC<-express_GMC[,c(1,2,6,4,5)]
express_GMC_all$gene <- paste("GMC_",express_GMC_all$ecNum,sep = "")
express_GMC_all <- express_GMC_all[,-6]

### ---- GH
express_GH_all<- express_all[which(express_all$protID %in% GH_Coromn$`Protein Id`),]
express_GH_all$ecNum<- GH_Coromn[match(express_GH_all$protID,GH_Coromn$`Protein Id`),]$ecnum
express_GH_all$cazy<- GH_Coromn[match(express_GH_all$protID,GH_Coromn$`Protein Id`),]$Cazy
express_GH_all$combo<- if_else(is.na(express_GH_all$ecNum),express_GH_all$cazy,express_GH_all$ecNum)
plot(as.factor(express_GH_all$protID),express_GH_all$TPM)
express_GH_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_GH_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_GH_by_ID<- express_GH_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_GH<- express_GH_all |>  group_by(level,Block,combo) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_GH$gene <- paste("GH_",express_GH$combo,sep = "")
express_GH<-express_GH[,c(1,2,6,4,5)]
express_GH_all$gene <- paste("GH",express_GH_all$ecNum,express_GH_all$cazy,sep = "_")
express_GH_all <- express_GH_all[,-c(6,7,8)]

### ---- OrgN transporters
express_OrgN_all<- express_all[which(express_all$protID %in% OrgN$proteinId),]
express_OrgN_all$ecNum<- OrgN[match(express_OrgN_all$protID,OrgN$proteinId),]$gene
plot(as.factor(express_OrgN_all$protID),express_OrgN_all$TPM)
express_OrgN_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_OrgN_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_OrgN_by_ID<- express_OrgN_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_OrgN<- express_OrgN_all |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_OrgN$gene <- paste("OrgN_",express_OrgN$ecNum,sep = "")
express_OrgN<-express_OrgN[,c(1,2,6,4,5)]
express_OrgN_all$gene <- paste("OrgN_",express_OrgN_all$ecNum,sep = "")
express_OrgN_all <- express_OrgN_all[,-6]

##first there will be overlap in my full transporters DB so I should remove dups
transporters$type_TCDB<-paste(transporters$type,transporters$TCDB,sep = "_")
transporters<-transporters[-which(transporters$`Protein Id` %in% c(express_OrgN_all$protID,express_AAAP_all$protID,express_OPT_all$protID)),]
express_tother_all<- express_all[which(express_all$protID %in% transporters$`Protein Id`),]
express_tother_all$ecNum<- transporters[match(express_tother_all$protID,transporters$`Protein Id`),]$type_TCDB
express_tother_all$sigP<-Coromn1_FilteredModels1_sigp6[match(express_tother_all$protID,Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob
express_tother_by_ID<- express_tother_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_tother<- express_tother_all |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_tother$gene <- paste("N_",express_tother$ecNum,sep = "")
express_tother<-express_tother[,c(1,2,6,4,5)]
express_tother_all$gene <- paste("N_",express_tother_all$ecNum,sep = "")
express_tother_all <- express_tother_all[,-6]


### ---- Chitin breakdown
express_CHIT_all<- express_all[which(express_all$protID %in% CHIT$proteinId),] 
express_CHIT_all$ecNum<- CHIT[match(express_CHIT_all$protID,CHIT$proteinId),]$ecNum
express_CHIT_all$sigP<- CHIT[match(express_CHIT_all$protID,CHIT$proteinId),]$sigP
express_CHIT_all$gene<- CHIT[match(express_CHIT_all$protID,CHIT$proteinId),]$gene
plot(as.factor(express_CHIT_all$protID),express_CHIT_all$TPM)
#express_CHIT_sec <- express_CHIT_all[-which(is.na(express_CHIT_all$sigP)),]
#express_CHIT_sec <- express_CHIT_sec[-which(express_CHIT_sec$sigP < 0.8),]
express_CHIT_by_ID<- express_CHIT_all |>  group_by(protID) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_CHIT<- express_CHIT_all |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))
express_CHIT$gene <- paste("CHIT_",express_CHIT$ecNum,sep = "")
express_CHIT<-express_CHIT[,c(1,2,6,4,5)]
express_CHIT_all$gene <- paste("CHIT_",express_CHIT_all$ecNum,sep = "")
express_CHIT_all <- express_CHIT_all[,-6]

### ---- match and combine

express_interest$gene<-GENES[match(express_interest$protID,GENES$proteinId),]$gene#add a column that matches the protienIDs to the gene names
express_interest[which(express_interest$protID == "1688323"),7]<-rep("NAGt",nrow(express_interest[which(express_interest$protID == "1688323"),7]))

express_interest_full <- rbind(express_interest,express_GMC_all,express_prot_all,express_AAAP_all,express_OPT_all,express_lac_all,express_OrgN_all,express_tother_all,express_CHIT_all,express_GH_all)
express_interest_full$prot_gene <-paste(express_interest_full$protID,express_interest_full$gene,sep = "_")

##remove BETA, GLUC, that are not secreted 
#express_interest<- express_interest[-which(express_interest$gene == "BETA" & is.na(express_interest$sigP)),]## all that have sigP are above 0.9
#express_interest<- express_interest[-which(express_interest$gene == "GLUC" & is.na(express_interest$sigP)),]## 
#express_interest<- express_interest[-which(express_interest$gene == "GLUC" & express_interest$sigP < 0.8),]


express_interest <- express_interest |>  group_by(level,Block,gene) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))  #add together reads from several gene copies
express_interest <- rbind(express_interest,express_lac,express_prot,express_AAAP,express_OPT,express_GMC,express_OrgN,express_tother,express_CHIT,express_GH)

###----reference genes for normalization
express_REF<- express_all[which(express_all$protID %in% REF$protID),] #
express_REF$gene <- REF[match(express_REF$protID,REF$protID),]$gene
colnames(express_REF)<-c("level","protID","sum_TPM","sum_reads","Block","gene")
express_REF <- rbind(express_REF[,-2],express_interest[which(express_interest$gene == "KGD"),c(1,4,5,2,3)])

express_interest$levelblock<-paste(express_interest$level,express_interest$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_interest_full$levelblock<-paste(express_interest_full$level,express_interest_full$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_REF$levelblock<-paste(express_REF$level,express_REF$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping

express_ref_wide<-express_REF[,-c(2)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)

express_interest<-cbind(express_interest,express_ref_wide[match(express_interest$levelblock,express_ref_wide$levelblock),4:9]) # add column with the number of BTub reads per level and block
express_interest_full<-cbind(express_interest_full,express_ref_wide[match(express_interest_full$levelblock,express_ref_wide$levelblock),4:9])

express_interest$UBC <- express_interest$sum_reads/express_interest$UBC
express_interest$TEF_2 <- express_interest$sum_reads/express_interest$TEF_2
express_interest$TEF_1 <- express_interest$sum_reads/express_interest$TEF_1
express_interest$BTub_1 <- express_interest$sum_reads/express_interest$BTub_1
express_interest$BTub_2 <- express_interest$sum_reads/express_interest$BTub_2
express_interest$KGD <- express_interest$sum_reads/express_interest$KGD

express_interest_full$UBC <- express_interest_full$NumReads/express_interest_full$UBC
express_interest_full$TEF_2 <- express_interest_full$NumReads/express_interest_full$TEF_2
express_interest_full$TEF_1 <- express_interest_full$NumReads/express_interest_full$TEF_1
express_interest_full$BTub_1 <- express_interest_full$NumReads/express_interest_full$BTub_1
express_interest_full$BTub_2 <- express_interest_full$NumReads/express_interest_full$BTub_2
express_interest_full$KGD <- express_interest_full$NumReads/express_interest_full$KGD


write_csv(express_interest,"clean_data/gene_interest_express_nooverlap.csv") #write data
write_csv(express_interest_full,"clean_data/gene_interest_express_nooverlap_all.csv")
