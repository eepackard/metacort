library(readr)
library(tidyverse)


Coromn1_FilteredModels1_kog <- read_delim("raw_data/Coromn1_FilteredModels1_kog_2025-04-22.tab", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)


OPT_Coromn <- read_delim("raw_data/OPT_2.A.67.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                                `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                                Location = col_skip(), Scaffold = col_skip(), 
                                                                                                                Start = col_skip(), End = col_skip(), 
                                                                                                                Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

AAAP_Coromn <- read_delim("raw_data/AAAP_2.A.18.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                        `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                        Location = col_skip(), Scaffold = col_skip(), 
                                                                                                        Start = col_skip(), End = col_skip(), 
                                                                                                        Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)




LAT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1287"),] ## this aligns perfectly with transporters DB assignment
ACT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1289"),] ## this aligns perfectly with transporters DB assignment
YATAAT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$kogid == "KOG1286"),] ## three are YAT and rest are AAT - should just call YAT/AAT
AAT<-YATAAT[-which(YATAAT$proteinId %in% c("1750305","1780949","1817633")),]
YAT<-YATAAT[which(YATAAT$proteinId %in% c("1750305","1780949","1817633")),]
POT<-Coromn1_FilteredModels1_kog[which(Coromn1_FilteredModels1_kog$proteinId == "1982532"),] #only one protein classified

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}

OrgN_trans<-AppendMe(c("LAT","ACT","AAT","YAT","POT"))
OrgN_trans$size<- if_else(OrgN_trans$gene %in% c("OPT","POT"),"Peptide","AminoAcid")

### I will pull in OPT and AAAP seperatley

write_csv(OrgN_trans,"clean_data/N_trans_clean.csv")
