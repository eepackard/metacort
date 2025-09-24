library(readr)
library(tidyverse)


Coromn1_FilteredModels1_kog <- read_delim("raw_data/Coromn1_FilteredModels1_kog_2025-04-22.tab", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)


LAT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1287"),]
ACT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1289"),]
YAT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1286"),]
AAAP<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid %in% c("KOG1303","KOG1305")),]

OPT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG2262"),]
POT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG2504"),]

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}

OrgN_trans<-AppendMe(c("LAT","ACT","YAT","AAAP","OPT","POT"))
OrgN_trans$size<- if_else(OrgN_trans$gene %in% c("OPT","POT"),"Peptide","AminoAcid")



write_csv(OrgN_trans,"clean_data/N_trans_clean.csv")
