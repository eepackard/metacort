library(readr)
library(tidyverse)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

Coromn1_FilteredModels1_sigp6 <- read_delim("raw_data/Coromn1_FilteredModels1_sigp6_2025-04-22.tab", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)

#endo-acting
## Chitinase = 3.2.1.14 
## Chitin deacetylase = 3.5.1.41
## chitin active lytic polymonosaccharide monooxygenases (LPMOs) encoding Lytic chitin monooxygenases = 1.14.99.53 - not present based on EC
#exo-acting
## B-N-acetylhexosaminidase = 3.2.1.52 (NAG)

ENDO<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum %in% c("3.2.1.14","3.5.1.41")),]
EXO<-Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$ecNum %in% c("3.2.1.52")),]
EXO<-EXO[c(1,6),]#only two unique

#combine
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}

CHIT<-AppendMe(c("ENDO","EXO"))


#now see if they are excreted

CHIT_sigP<-Coromn1_FilteredModels1_sigp6[which(Coromn1_FilteredModels1_sigp6$protein_id %in% CHIT$proteinId),]
CHIT$sigP<-rep(NA,nrow(CHIT))
CHIT[match(CHIT_sigP$protein_id,CHIT$proteinId),]$sigP <- CHIT_sigP$sp_prob



write_csv(CHIT,"clean_data/CHIT_clean.csv")
