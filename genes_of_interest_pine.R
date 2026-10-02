library(readr)
library(tidyverse)

##interested in 


express_all<-read_csv("clean_data/expression_salmon_quant_pine_all.csv")

Sweet<-read_delim("raw_data/pine_sweet_transcripts.txt",delim = "/t",col_names = c("trans"))
tub<-read_delim("raw_data/gamma-tubulin_transcripts.txt",delim= "/t",col_names = c("trans"))
actin<-read_delim("raw_data/actin_transcripts.txt", delim= "/t",col_names = c("trans"))
invert<-read_delim("raw_data/invertase_transcripts.txt", delim= "/t",col_names = c("trans"))

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}

GENES<-AppendMe(c("Sweet","tub","actin","invert"))
GENES<-GENES[match(unique(GENES$trans),GENES$trans),]

GENES$trans<-gsub("Pinsy_PS_","",GENES$trans)
GENES$trans<-gsub("PS_","",GENES$trans)

#pull out only highly expressed sweets
sweet_names<- read_csv("clean_data/low_sweets_pine.csv",skip = 1,col_names = c("names"))
sweet_names$names<- gsub("_Sweet","",sweet_names$names)

GENES<- GENES[-which(GENES$trans %in% sweet_names$names),]

###---- based on name match and combine
express_interest<-express_all[which(express_all$Transcript %in% GENES$trans),] #find only the genes of interest

express_interest$gene<-GENES[match(express_interest$Transcript,GENES$trans),]$gene#add a column that matches the protienIDs to the gene names

express_interest_full <- express_interest
express_interest_full$prot_gene <-paste(express_interest_full$Transcript,express_interest_full$gene,sep = "_")


express_interest <- express_interest |>  group_by(level,Block,gene) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))  #add together reads from several gene copies


write_csv(express_interest,"clean_data/gene_interest_express_pine_GC.csv") #write data
write_csv(express_interest_full,"clean_data/gene_interest_express_pine_GC_all.csv")
