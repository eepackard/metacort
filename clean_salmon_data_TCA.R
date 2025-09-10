library(readr)
library(tidyverse)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

##interested in pulling out all genes in the TCA cycle 

TCA<-Coromn1_FilteredModels1_ec[grep("TCA", Coromn1_FilteredModels1_ec$pathway),]

#samples 1 A1, 2 D2, 3 B1, 4 D3, 7 A4, 8 D1, 9 D4, 11 B1, 12 A4, 13 C2, 15 A2, 17 B3, 18 B1, 19 D2, 20 D1
high_list<-read_delim("raw_data/high.list",delim = "\t",trim_ws = TRUE,col_names = FALSE)
high_list<-as.list(t(high_list))
high_list_files<-as.list(gsub("quant_salmon_quantYC-4261-","salmon_quant_GC/quant_",high_list))#add path

#samples 1 C1, 2 A4, 3 D4, 4 B2, 7 C1, 8 D3, 9 B2, 11 B4, 12 C2, 13 B3, 15 D1, 17 C4, 18 A4, 19 A2, 20 C3
low_list<-read_delim("raw_data/low.list",delim = "\t",trim_ws = TRUE,col_names = FALSE)
low_list<-as.list(t(low_list))
low_list_files<-as.list(gsub("quant_salmon_quantYC-4261-","salmon_quant_GC/quant_",low_list))#add path

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
#combine high and low from same blocks and then clear out a few unneeded variables
for (i in 1:15) {
quant_full[[i]]<-bind_rows(quant_low[[i]],quant_high[[i]])  
quant[[i]]<-quant_full[[i]][,c(1,3,7,8)] 
}

##name the items in list as reps - note that it does not actually match the block numbers assigned during collection!
names(quant)<-gsub("-.*","",gsub("salmon_quant_GC/quant_","block",low_list_files))

quant <- quant |> map2(names(quant), ~.x |> mutate(Block = .y)) #add column that has the block name
express_all<-bind_rows(quant)#merge all reps into one dataframe

express_interest<-express_all[which(express_all$protID %in% TCA$proteinId),] #find only the genes of interest
express_BTub<- express_all[which(express_all$protID %in% c("1845132","1851235")),] #find the data for Beta tubulin
express_BTub <- express_BTub |>  group_by(level,Block) |> summarise(sum_reads = sum(NumReads)) #there is two copies so I can add together data

express_interest$ecNum<-TCA[match(express_interest$protID,TCA$proteinId),]$ecNum#add a column that matches the protienIDs to the gene names

express_interest_TPM <- express_interest |>  group_by(level,Block,ecNum) |> summarise(sum_TPM = sum(TPM))  #add together reads from several gene copies
express_interest_reads <- express_interest |>  group_by(level,Block,ecNum) |> summarise(sum_reads = sum(NumReads)) #add together reads from several gene copies

express_interest_reads$levelblock<-paste(express_interest_reads$level,express_interest_reads$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping
express_BTub$levelblock<-paste(express_BTub$level,express_BTub$Block,sep = "") #to match the BTub data and reads data need to have a column for grouping

express_interest_reads$BTub<-express_BTub[match(express_interest_reads$levelblock,express_BTub$levelblock),]$sum_reads # add column with the number of BTub reads per level and block

express_interest_norm <-express_interest_reads |> mutate(norm_reads = sum_reads/BTub) #normalize by dividing the num of reads by num of BTub reads

express_interest_summary<-cbind(express_interest_TPM,express_interest_reads$sum_reads,express_interest_norm$norm_reads) #bind back the TPM, reads, and norm
colnames(express_interest_summary)<-c(colnames(express_interest_TPM),"sum_reads","norm_reads") #rename columns

write_csv(express_interest_summary,"clean_data/gene_interest_express_TCA.csv") #write data

