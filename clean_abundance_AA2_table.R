library(readr)
library(tidyverse)
library(plyr)

#this is to clean and inspect the abundance table (salmon quantifited TPM of HMM assembled contigs of AA2 transcripts)

#load data
abun<- read_delim("raw_data/merged_abundance_table_blasted_AA2.csv",delim = ";")
abun<- abun[-which(abun$...1 == "_"),]
row.names(abun)<- abun$...1

abun[abun < 2] <- 0
abun <- abun[-which(rowSums(abun[,-1]) < 11),]

abun <- separate(abun,col = "...1",sep = "_",into = c("node","node_num","prot_blast","extra","extra2","extra3"))
abun$prot_blast <- paste(abun$prot_blast,abun$extra,abun$extra2,abun$extra3,sep = "")
abun$node <- paste(abun$node,abun$node_num,sep = "_")
abun$prot_blast <- gsub("NA","",abun$prot_blast)
abun<-abun[,-grep("extra|node_num",colnames(abun))]
abun_agg<- ddply(abun, "prot_blast", numcolwise(sum))

abun_agg_ge <-abun_agg
abun_agg_ge <- separate(abun_agg_ge,col = "prot_blast",into = c("genus","sp","other"),sep = c(3,6))
abun_agg_ge_sum<- aggregate(.~genus,data = abun_agg_ge[,-c(2,3)], sum,na.rm=TRUE)
abun_agg_ge_sum$genus <- abun_agg_ge[match(abun_agg_ge_sum$sp,abun_agg_ge$sp),]$genus

abun_coromn <- abun_agg[grep("Coromn",abun_agg$prot_blast),]

total_AA2 <- c("tot_AA2",colSums(abun[,3:ncol(abun)]))

total_coromn <- c("tot_coromn",colSums(abun_coromn[,-1]))

HMM_results<-as.data.frame(rbind(abun_agg_ge_sum,total_AA2,total_coromn))
coln<- HMM_results$genus
HMM_results <- as.data.frame(t(HMM_results[,-1]))
colnames(HMM_results)<-coln
HMM_results <- HMM_results |> mutate_if(is.character,as.numeric)
HMM_results$prop_corom <- HMM_results$tot_coromn/HMM_results$tot_AA2
HMM_results$prop_cort <- HMM_results$Cor/HMM_results$tot_AA2
HMM_results$sample <- rownames(HMM_results)
HMM_results <- separate(HMM_results,col = "sample",sep = "-",into = c("level","block","extra","extra2","extra3")) 
HMM_results$block <- paste("block",HMM_results$block,sep = "")
HMM_results$levelblock <- paste(HMM_results$level,HMM_results$block,sep = "")
HMM_results <- HMM_results[,-grep("extra",colnames(HMM_results))]

write_csv(abun_agg,"clean_data/filtered_sum_AA2_abun.csv")
write_csv(HMM_results,"clean_data/HMM_results.csv")
