library(readr)
library(tidyverse)
library(lme4)
library(car)
library(ggpubr)

#run statistics
norm_wide<-read_csv("clean_data/norm_wide_clean.csv")#this is TPM

list_gene<-list(c("ALAS" , "APEP", "BETA" ,"GLUC" , "GLY","GT48" ,"KGD" ,  "LACC" , "PROT" , "GMC_AOx",  "GMC_AAO", "GMC_GDH"  
   , "OrgN_AAAP"  ,   "OrgN_ACT","OrgN_LAT","OrgN_OPT","OrgN_YAT","OrgN_AAT","CHIT_3.2.1.14" ,"CHIT_3.2.1.52"
     ,"CHIT_3.5.1.41", "TPS","TPP","TASE","TASE_TPP" ,"CUE","MnP" ))# exlcuding POT because it is barely expressed (only in 7 of 28 samples)...

##check normality and transform ----
results<-list()
p_value_ass<-list()

for (i in 1:length(list_gene[[1]])){
  results[[i]]<-shapiro.test(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene[[1]][[i]])]))
  p_value_ass[[i]]<-results[[i]]$p.value
  hist(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene[[1]][[i]])]),main = list_gene[[1]][[i]])
  hist(log10(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene[[1]][[i]])])),main = paste("log", list_gene[[1]][[i]]))
  hist(sqrt(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene[[1]][[i]])])),main = paste("sqrt", list_gene[[1]][[i]]))
}  
##based on histograms - seems chit3.2.1.14 and chit 3.2.1.52 shouldnt be logged it only makes them worse - they are already quite normal
#CUE is also already okay
## MnP is very right skew and log kinda just makes left skew - think a sqrt would be better

results<-data.frame(round(unlist(p_value_ass),3),unlist(list_gene))
non_norm<-results[which(results$round.unlist.p_value_ass...3. < 0.05),]
##13 are not normal - 7 are okay - APEP, BETA, KGD, PROT,CHIT, CUE - but prot, BETA, and APEP still look much better if logged, KGD doesnt make too much difference

unlist(list_gene)

norm_wide_trans <- norm_wide |> mutate_at(vars(matches(unlist(list_gene)[-c(25,26)])),sqrt)# - only one i will not transform is CUE and TASE/TPP - ratios

results<-list()
p_value_ass<-list()

for (i in 1:length(list_gene[[1]])){
  results[[i]]<-shapiro.test(unlist(norm_wide_trans[,which(colnames(norm_wide_trans) %in% list_gene[[1]][[i]])]))
  p_value_ass[[i]]<-results[[i]]$p.value
}  

results<-data.frame(round(unlist(p_value_ass),3),unlist(list_gene))
logged_non_norm<-results[which(results$round.unlist.p_value_ass...3. < 0.05),] 
## after looging GLY, GT48 are still not normal... 

##run linear models ----
mods<-list()
summs<-list()
coeff<-list()
est<-list()
Stder<-list()
t_val<-list()
atable<-list()
p_value_ano<-list()
plots<-list()
for (i in 1:26){ ##MnP is 27 in the list so i should be able to just run loop only 26 times
  mods[[i]]<-lmer(unlist(norm_wide_trans[,which(colnames(norm_wide_trans) %in% list_gene[[1]][[i]])])~MnP+(1|Block),data= norm_wide_trans)
  summs[[i]]<-summary(mods[[i]])
  coeff[[i]]<-summs[[i]]$coefficients
  est[[i]]<-coeff[[i]][2]
  Stder[[i]]<-coeff[[i]][4]
  t_val[[i]]<-coeff[[i]][6]
  atable[[i]]<-Anova(mods[[i]])
  p_value_ano[[i]]<-atable[[i]]$`Pr(>Chisq)`
  plots[[i]]<-plot(mods[[i]])
  }

ggarrange(plotlist =plots[20:25])
##not unreasonable - including GT48 which is only one not "normal"

unlist(list_gene)[grep("singular",summs)]#LACC, OrgN_ACT, OrgN_AAT, OrgN_AAAP, GT48 

mod_resuls<-data.frame(unlist(list_gene)[-26],round(unlist(est),digits = 4),round(unlist(Stder),digits = 4),round(unlist(t_val),digits = 4),round(unlist(p_value_ano),digits = 4))
colnames(mod_resuls)<-c("gene","estimate","Std.error","t_value","p_value")

mod_resuls$adjusted_p<- p.adjust(mod_resuls$p_value,method = "fdr")

write_csv(mod_resuls,file = "clean_data/lm_results.csv")

