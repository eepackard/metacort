library(readr)
library(tidyverse)



#run statistics on the following

# MnP vs PROT, GLy/CRO, GT48, CHIT/NAG, SOD, CAT, ALAS, APEP

norm_wide<-read_csv("clean_data/norm_wide_clean.csv")
norm_wide_sub<-read_csv("clean_data/norm_wide_sub_clean.csv")

#paired t-test

#assumptions
list_gene<-list(c(colnames(norm_wide[,-c(1:3,10)])))
results<-list()
p_value<-list()
for (i in 1:26){
  results[[i]]<-shapiro.test(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene[[1]][[i]])]))
  p_value[[i]]<-results[[i]]$p.value
}

results<-data.frame(round(unlist(p_value),3),unlist(list_gene))
non_norm<-results[which(results$round.unlist.p_value...3. < 0.05),]

wilcox.test(norm_wide[which(norm_wide$level == "high"),]$MnP,norm_wide[which(norm_wide$level == "low"),]$MnP,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$PROT,norm_wide[which(norm_wide$level == "low"),]$PROT,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$GLY,norm_wide[which(norm_wide$level == "low"),]$GLY,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$GT48,norm_wide[which(norm_wide$level == "low"),]$GT48,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$CHIT_ENDO,norm_wide[which(norm_wide$level == "low"),]$CHIT_ENDO,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$CHIT_EXO,norm_wide[which(norm_wide$level == "low"),]$CHIT_EXO,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$APEP,norm_wide[which(norm_wide$level == "low"),]$APEP,paired = TRUE)
wilcox.test(norm_wide[which(norm_wide$level == "high"),]$ALAS,norm_wide[which(norm_wide$level == "low"),]$ALAS,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$SOD,norm_wide[which(norm_wide$level == "low"),]$SOD,paired = TRUE)
wilcox.test(norm_wide[which(norm_wide$level == "high"),]$LACC,norm_wide[which(norm_wide$level == "low"),]$LACC,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$CAT,norm_wide[which(norm_wide$level == "low"),]$CAT,paired = TRUE)
wilcox.test(norm_wide[which(norm_wide$level == "high"),]$GMC_sum,norm_wide[which(norm_wide$level == "low"),]$GMC_sum,paired = TRUE)
wilcox.test(norm_wide[which(norm_wide$level == "high"),]$APET,norm_wide[which(norm_wide$level == "low"),]$APET,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$AA,norm_wide[which(norm_wide$level == "low"),]$AA,paired = TRUE)
t.test(norm_wide[which(norm_wide$level == "high"),]$CUE,norm_wide[which(norm_wide$level == "low"),]$CUE,paired = TRUE)


#reduced to top 8 ----

#assumptions
list_gene<-list(c(colnames(norm_wide_sub[,-c(1:3,10)])))
results<-list()
p_value<-list()
for (i in 1:26){
results[[i]]<-shapiro.test(unlist(norm_wide_sub[,which(colnames(norm_wide_sub) %in% list_gene[[1]][[i]])]))
p_value[[i]]<-results[[i]]$p.value
}

results<-data.frame(round(unlist(p_value),3),unlist(list_gene))
non_norm<-results[which(results$round.unlist.p_value...3. < 0.05),]

wilcox.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$MnP,norm_wide_sub[which(norm_wide_sub$level == "low"),]$MnP,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$PROT,norm_wide_sub[which(norm_wide_sub$level == "low"),]$PROT,paired = TRUE)
wilcox.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$GLY,norm_wide_sub[which(norm_wide_sub$level == "low"),]$GLY,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$GT48,norm_wide_sub[which(norm_wide_sub$level == "low"),]$GT48,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$CHSN,norm_wide_sub[which(norm_wide_sub$level == "low"),]$CHSN,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$CHIT_ENDO,norm_wide_sub[which(norm_wide_sub$level == "low"),]$CHIT_ENDO,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$CHIT_EXO,norm_wide_sub[which(norm_wide_sub$level == "low"),]$CHIT_EXO,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$APEP,norm_wide_sub[which(norm_wide_sub$level == "low"),]$APEP,paired = TRUE)
wilcox.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$ALAS,norm_wide_sub[which(norm_wide_sub$level == "low"),]$ALAS,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$SOD,norm_wide_sub[which(norm_wide_sub$level == "low"),]$SOD,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$CAT,norm_wide_sub[which(norm_wide_sub$level == "low"),]$CAT,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$GMC_sum,norm_wide_sub[which(norm_wide_sub$level == "low"),]$GMC_sum,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$AA,norm_wide_sub[which(norm_wide_sub$level == "low"),]$AA,paired = TRUE)
t.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$APET,norm_wide_sub[which(norm_wide_sub$level == "low"),]$APET,paired = TRUE)
wilcox.test(norm_wide_sub[which(norm_wide_sub$level == "high"),]$CUE,norm_wide_sub[which(norm_wide_sub$level == "low"),]$CUE,paired = TRUE)
