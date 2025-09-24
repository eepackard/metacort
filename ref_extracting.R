library(readr)
library(tidyverse)


##potential reference
##beta tubulin = prot numbers 1845132 and 1851235
##TEF-1a = prot number 1980816 and 1753879
##ubc = prot 	949941

express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")


Coromn1_FilteredModels1_kog <- read_delim("raw_data/Coromn1_FilteredModels1_kog_2025-04-22.tab", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)


BTub_1<- express_all[which(express_all$protID %in% c("1845132")),] #,"1851235"find the data for Beta tubulin
BTub_2<- express_all[which(express_all$protID %in% c("1851235")),]

TEF_1<- express_all[which(express_all$protID %in% c("1980816")),] #
TEF_2<- express_all[which(express_all$protID %in% c("1753879")),]

UBC<- express_all[which(express_all$protID %in% c("949941")),] #

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), gene = x)
  }))
}
GENES<-AppendMe(c("BTub_1","BTub_2","TEF_1","TEF_2","UBC"))

plot(as.factor(BTub_1$level),BTub_1$NumReads)
plot(as.factor(BTub_2$level),BTub_2$NumReads)
plot(as.factor(BTub_1$Block),BTub_1$NumReads)
plot(as.factor(BTub_2$Block),BTub_2$NumReads)
plot(as.factor(BTub_1$level),BTub_1$TPM)
plot(as.factor(BTub_2$level),BTub_2$TPM)

plot(as.factor(TEF_1$level),TEF_1$NumReads)
plot(as.factor(TEF_2$level),TEF_2$NumReads)
plot(as.factor(TEF_1$level),TEF_1$TPM)
plot(as.factor(TEF_2$level),TEF_2$TPM)

plot(as.factor(UBC$level),UBC$NumReads)
plot(as.factor(UBC$level),UBC$TPM)

plot(BTub_1$NumReads,TEF_1$NumReads)
plot(BTub_2$NumReads,TEF_1$NumReads)
plot(BTub_2$NumReads,TEF_2$NumReads)
plot(TEF_1$NumReads,TEF_2$NumReads)
plot(BTub_1$NumReads,BTub_2$NumReads)
plot(UBC$NumReads,TEF_1$NumReads)
plot(UBC$NumReads,TEF_2$NumReads)
plot(UBC$NumReads,BTub_1$NumReads)
plot(UBC$NumReads,BTub_2$NumReads)

plot(UBC$NumReads,UBC$TPM)
plot(TEF_1$NumReads,TEF_1$TPM)
plot(TEF_2$NumReads,TEF_2$TPM)
plot(BTub_1$NumReads,BTub_1$TPM)
plot(BTub_2$NumReads,BTub_2$TPM)

GENES$gene<-as.factor(GENES$gene)
GENES$level<-as.factor(GENES$level)
GENES$Block<-as.factor(GENES$Block)

wide_TPM<-GENES[,-c(2,4)] |>  pivot_wider(names_from = "gene",values_from = TPM)
wide_TPM$levelblock <- paste(wide_TPM$level,wide_TPM$Block,sep = "")
wide_reads<-GENES[,-c(2,3)] |>  pivot_wider(names_from = "gene",values_from = NumReads)
wide_reads$levelblock <- paste(wide_reads$level,wide_reads$Block,sep = "")

write_csv(GENES,"clean_data/ref_clean.csv")
