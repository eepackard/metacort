library(readr)
library(tidyverse)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)

#3.4.11

APEP<-Coromn1_FilteredModels1_ec[which(grepl("3.4.11.",Coromn1_FilteredModels1_ec$ecNum)),]
CPEP<-Coromn1_FilteredModels1_ec[which(grepl("3.4.16.",Coromn1_FilteredModels1_ec$ecNum)),]#serine type carbo
MCPEP<-Coromn1_FilteredModels1_ec[which(grepl("3.4.17.",Coromn1_FilteredModels1_ec$ecNum)),]
DTPEP<-Coromn1_FilteredModels1_ec[which(grepl("3.4.14.",Coromn1_FilteredModels1_ec$ecNum)),]

#combine
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}

PEP<-AppendMe(c("APEP","CPEP","MCPEP","DTPEP"))


#now see if they are excreted

PEP$sigP<-Coromn1_FilteredModels1_sigp6[match(PEP$proteinId, Coromn1_FilteredModels1_sigp6$protein_id),]$sp_prob

write_csv(PEP,"clean_data/PEP_clean.csv")
