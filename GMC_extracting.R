library(readr)
library(tidyverse)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)
Coromn1_FilteredModels1_kog <- read_delim("raw_data/Coromn1_FilteredModels1_kog_2025-04-22.tab", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)

jgi_alignment_hits_GMC <- read_csv("raw_data/jgi_alignment_hits_GMC.csv")

## alcohol oxidase AOx (H2O2 production?) = 1.1.3.13 , AA3_3
## cellobiose dehydrogenase CDH = 1.1.99.18, AA3_1
## arly-alcohol oxidase AAO = 1.1.3.7, AA3_2
## glucose oxidase GOx 1.1.3.4 AA3_2
##glucose dehydrogenase GDH = 1.1.5.9 AA3_2
## pyranose dehydrogenase PDH 1.1.99.29 AA3_2 
## pyranose oxidase POx = 1.1.3.13 AA3_4

##KOG1238, 

EC_GMC<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum %in% c("1.1.3.13","1.1.99.18","1.1.3.7","1.1.3.4","1.1.5.9","1.1.99.29","1.1.3.13")),]
KOG_GMC<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1238"),]

#clean and combine
overlap<-EC_GMC[which(EC_GMC$proteinId %in% KOG_GMC$proteinId),]
extra<-EC_GMC[-which(EC_GMC$proteinId %in% KOG_GMC$proteinId),]
overlap$kogid <- rep("KOG1238",nrow(overlap))
overlap<- overlap[,which(colnames(overlap) %in% c("proteinId","ecNum" ,"definition","kogid"))]

extra$kogid <- rep(NA,nrow(extra))
extra <- extra[,which(colnames(extra) %in% c("proteinId","ecNum" ,"definition","kogid"))]

EC_GMC<- rbind(extra,overlap)

kog_extra<- KOG_GMC[-which(KOG_GMC$proteinId %in% EC_GMC$proteinId),]
kog_extra$ecNum <- rep(NA,nrow(kog_extra))
kog_extra <- kog_extra[,c(2,7,4,3)]
colnames(kog_extra) <- c("proteinId","ecNum" ,"definition","kogid")

GMC<- rbind(EC_GMC,kog_extra)
length(unique(GMC$proteinId)) == nrow(GMC)

##check blast results
#using the reference sequences identified in Matilla

#keep only with pident greater than 35 %
jgi_alignment_hits_GMC$Hit <- gsub("Coromn1\\|","",jgi_alignment_hits_GMC$Hit)
jgi_alignment_hits_GMC$`% Hit Identity`<-as.numeric(gsub("%","",jgi_alignment_hits_GMC$`% Hit Identity`))

jgi_alignment_hits_GMC<-jgi_alignment_hits_GMC[which(jgi_alignment_hits_GMC$`% Hit Identity` >35),]
length(unique(jgi_alignment_hits_GMC$Hit)) #so only 23 actually 

length(which(GMC$proteinId %in% jgi_alignment_hits_GMC$Hit))#and all of those I already identifed

##now I can select which GMC it is based on highest blast - changing to score
GMC_protID_list<-as.list(GMC$proteinId)
Blast_result<-list()
for (i in 1:nrow(GMC)){
Blast_result[[i]]<-jgi_alignment_hits_GMC[which(jgi_alignment_hits_GMC$Hit == GMC_protID_list[[i]]),]
Blast_result[[i]]<-Blast_result[[i]][which(Blast_result[[i]]$Score == max(Blast_result[[i]]$Score)),c(1,3,8,9)]
}

blast<-bind_rows(Blast_result)
GMC$blast<- blast[match(GMC$proteinId,blast$Hit),]$`Query Name`

GMC$ecNum <-if_else(grepl("Alcohol oxidase",GMC$blast),"1.1.3.13",if_else(grepl("Cellobiose dehydrogenase",GMC$blast),"1.1.99.18",if_else(grepl("Aryl-alcohol",GMC$blast),"1.1.3.7",GMC$ecNum)))


#not in blast but assigned by EC and/or KOG
nohit_GMC<-GMC[-which(GMC$proteinId %in% jgi_alignment_hits_GMC$Hit),]

#now see if they are excreted

GMC_sigP<-Coromn1_FilteredModels1_sigp6[which(Coromn1_FilteredModels1_sigp6$protein_id %in% GMC$proteinId),]
GMC$sigP<-rep(NA,nrow(GMC))
GMC[match(GMC_sigP$protein_id,GMC$proteinId),]$sigP <- GMC_sigP$sp_prob
GMC[which(GMC$proteinId %in% GMC_sigP$protein_id),]$sigP <- c(rep("sigP",11))

## I made updates to gene models for protein 1984530 and 1749753, which had a chimera and some exon/intron gaps missed, respectively
##manually change the proteinID numbers so that they will match with the new salmon results - this doesn't change their tree placement 

GMC[which(GMC$proteinId == "1984530"),]$proteinId <- 1850569
GMC[which(GMC$proteinId == "1749753"),]$proteinId <- 1682818

write_csv(GMC,"clean_data/GMC_clean.csv")

##I then also made a tree with these sequences and based on how this was I made some adjustments to the classificatoin

GMC<-read_delim("clean_data/GMC_clean_plus_tree_annotation.csv",delim = ";")
#the ones with NA were short and kinda have one motif but not confident of there annotation - fragements or allelic variants 
GMC<-GMC[-which(is.na(GMC$`AA_#`)),]

write_csv(GMC,"clean_data/GMC_clean.csv")
