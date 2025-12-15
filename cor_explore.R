library(readr)
library(tidyverse)

#read in the expression data and the data from JGI

express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")


Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)


A1_prot_Coromn <- read_delim("raw_data/A1_prot_Coromn.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                                `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                                Location = col_skip(), Scaffold = col_skip(), 
                                                                                                                Start = col_skip(), End = col_skip(), 
                                                                                                                Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

Coromn1_FilteredModels1_go_2025_04_22 <- read_delim("raw_data/Coromn1_FilteredModels1_go_2025-04-22.tab", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)

Coromn1_FilteredModels1_domain_2025_04_22 <- read_delim("raw_data/Coromn1_FilteredModels1_domain_2025-04-22.tab", 
                                                        delim = "\t", escape_double = FALSE, 
                                                        trim_ws = TRUE)

OrgN<-read_csv("clean_data/N_trans_clean.csv")
GMC<-read_csv("clean_data/GMC_clean.csv")

#need to remove block19
express_all<-express_all[-which(express_all$Block == "block19"),]

#want to pull out genes that either have EC, I have pulled out otherwise (GMC, OrgN), or are proteases (maybe have ec)

express_all_EC<- express_all[which(express_all$protID %in% Coromn1_FilteredModels1_ec$proteinId),]
Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% A1_prot_Coromn$`Protein Id`),] #EC 3.4.23. are aspartic proteases
A1_prot_Coromn[-which(A1_prot_Coromn$`Protein Id` %in% Coromn1_FilteredModels1_ec$proteinId),] ## there is five prots that do not have EC but are A1
express_all_EC$ecNum <- Coromn1_FilteredModels1_ec[match(express_all_EC$protID,Coromn1_FilteredModels1_ec$proteinId),]$ecNum

express_all_exprot<-express_all[which(express_all$protID %in% A1_prot_Coromn[-which(A1_prot_Coromn$`Protein Id` %in% Coromn1_FilteredModels1_ec$proteinId),]$`Protein Id`),]
express_all_exprot$ecNum <- rep("3.4.23.-",nrow(express_all_exprot))# not sure exactly which type they fall into

Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% GMC$proteinId),] #all 1.1.3.13 are already there
express_all_GMC<- express_all[which(express_all$protID %in% GMC[-which(GMC$proteinId %in% Coromn1_FilteredModels1_ec$proteinId),]$proteinId),]
express_all_GMC$ecNum <- GMC[match(express_all_GMC$protID,GMC$proteinId),]$ecNum# not sure exactly which type they fall into

Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% OrgN$proteinId),] #non have EC
express_all_orgN<- express_all[which(express_all$protID %in% OrgN$proteinId),]
express_all_orgN$ecNum <- OrgN[match(express_all_orgN$protID,OrgN$proteinId),]$gene# not sure exactly which type they fall into


#bind all
express_select<-rbind(express_all_EC,express_all_exprot,express_all_orgN,express_all_GMC)

#will try to methods one that is prot by prot and one with them combined by EC

#make summed df 
express_select_sum<- express_select |>  group_by(ecNum,level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))

#add levelblock factor
express_all$levelblock <- paste(express_all$level,express_all$Block,sep = "")
express_select_sum$levelblock <- paste(express_select_sum$level,express_select_sum$Block,sep = "")

#make wide dfs
express_wide_sum<- express_select_sum[,-c(2,3,5)] |>  pivot_wider(names_from = "ecNum",values_from = sum_TPM)
express_wide<- express_all[,-c(1,4,5)] |>  pivot_wider(names_from = "protID",values_from = TPM)


#pull out MnP
MnP_sum<-express_wide_sum[,which(colnames(express_wide_sum) == "1.11.1.13")]
colnames(MnP_sum)<-"MnP"
express_wide_sum<-express_wide_sum[,-which(colnames(express_wide_sum) == "1.11.1.13")]

MnP<-as.data.frame(rowSums(express_wide[,which(colnames(express_wide) %in% unique(express_all_EC[which(express_all_EC$ecNum == "1.11.1.13"),]$protID))]))
colnames(MnP)<-"MnP"
express_wide<-express_wide[,-which(colnames(express_wide) %in% unique(express_all_EC[which(express_all_EC$ecNum == "1.11.1.13"),]$protID))]

#make first column row names
samples_sum<-express_wide_sum$levelblock
express_wide_sum<-express_wide_sum[,-1]
rownames(express_wide_sum)<-samples_sum

samples<- express_wide$levelblock
express_wide<-express_wide[,-1]
rownames(express_wide)<-samples

#remove any columns where the abudnace is less than 100
express_wide_sum <- express_wide_sum[,-which(colSums(express_wide_sum)< 100)]
express_wide <- express_wide[,-which(colSums(express_wide)< 100)]

#should log transform all to be safe

MnP$MnP<-log10(MnP$MnP)
MnP_sum$MnP<-log10(MnP_sum$MnP)

express_wide <- express_wide |> mutate_all(log10)
express_wide_sum <- express_wide_sum |> mutate_all(log10)

#now can test
cordfsum<-as.data.frame(cor(express_wide_sum,MnP_sum,method = "pearson"))
cordfsum$ec<-rownames(cordfsum)
cordfsum<-cordfsum[order(cordfsum$MnP,decreasing = TRUE),]
cordfsum$rank<- c(1:length(cordfsum$MnP))

cordf<-as.data.frame(cor(express_wide,MnP,method = "pearson"))
cordf$prot<-rownames(cordf)
cordf<-cordf[order(cordf$MnP,decreasing = TRUE),]
cordf$rank<- c(1:length(cordf$MnP))

#lets match some into to make it easier to inspect
#but some have multiple so i will combine
iprID<-aggregate(iprId ~ proteinId, Coromn1_FilteredModels1_domain_2025_04_22,FUN = "paste")
iprDesc<-aggregate(iprDesc ~ proteinId, Coromn1_FilteredModels1_domain_2025_04_22,FUN = "paste")
ipr<-cbind(iprID,iprDesc)
ipr<-ipr[,-3]

goname<-aggregate(go_name ~ proteinId, Coromn1_FilteredModels1_go_2025_04_22,FUN = "paste")
goacc<-aggregate(go_acc ~ proteinId, Coromn1_FilteredModels1_go_2025_04_22,FUN = "paste")
go<-cbind(goname,goacc)
go<-go[,-3]


cordf$goname<-go[match(cordf$prot,go$proteinId),]$go_name
cordf$goacc<-go[match(cordf$prot,go$proteinId),]$go_acc
cordf$iprID<-ipr[match(cordf$prot,ipr$proteinId),]$iprId
cordf$iprDesc<-ipr[match(cordf$prot,ipr$proteinId),]$iprDesc

cordfsum$def<- Coromn1_FilteredModels1_ec[match(cordfsum$ec,Coromn1_FilteredModels1_ec$ecNum),]$definition

##lets try to sort out the ones i looked at a priori
express_interest<-read_csv("clean_data/gene_interest_express_GC_all.csv")
express_interest$ecNum<-Coromn1_FilteredModels1_ec[match(express_interest$protID,Coromn1_FilteredModels1_ec$proteinId),]$ecNum

cordf$apriori<-if_else(cordf$prot %in% express_interest$protID, "red","black")
cordfsum$apriori<-if_else(cordfsum$ec %in% express_interest$ecNum, "red",if_else(cordfsum$ec %in% c("OPT","ACT","AAAP","YAT","LAT"),"red",if_else(cordfsum$ec %in% c("1.1.3.7","1.1.99.18"),"red","black")))


ggplot(cordf)+
  geom_point(aes(x=rank,y=MnP,colour = apriori),alpha=0.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()

ggplot(cordfsum)+
  geom_point(aes(x=rank,y=MnP,colour = apriori),alpha=0.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()
