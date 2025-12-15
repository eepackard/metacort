library(readr)
library(tidyverse)
library(lme4)
library(car)
library(ggpubr)

#run statistics
norm_wide<-read_csv("clean_data/norm_wide_clean.csv")#this is TPM
norm_wide_sub<-read_csv("clean_data/norm_wide_sub_clean.csv")


list_gene<-list(c("ALAS" , "APEP", "BETA" ,"GLUC" , "GLY","GT48" ,"KGD" ,  "LACC" , "PROT" , "GMC_1.1.3.13",  "GMC_1.1.3.7"  
   ,"GMC_1.1.99.18", "OrgN_AAAP"  ,   "OrgN_ACT","OrgN_LAT","OrgN_OPT","OrgN_YAT","CHIT_3.2.1.14" ,"CHIT_3.2.1.52"
     ,"CHIT_3.5.1.41" ,"CUE","MnP" ))# exlcuding POT because it is barely expressed...

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
##15 are not normal

results<-data.frame(round(unlist(p_value_ass),3),unlist(list_gene))
non_norm<-results[which(results$round.unlist.p_value_ass...3. < 0.05),]
##15 are not normal - 7 are okay - APEP, BETA, KGD, PROT,CHIT, CUE - but prot, BETA, and APEP still look much better if logged, KGD doesnt make too much difference

norm_wide_trans <- norm_wide |> mutate_at(vars(matches(unlist(list_gene)[-c(18,19,21,22)])),log10)
norm_wide_trans <- norm_wide_trans |> mutate_at(vars(matches(unlist(list_gene)[c(22)])),sqrt)

results<-list()
p_value_ass<-list()

for (i in 1:length(list_gene[[1]])){
  results[[i]]<-shapiro.test(unlist(norm_wide_trans[,which(colnames(norm_wide_trans) %in% list_gene[[1]][[i]])]))
  p_value_ass[[i]]<-results[[i]]$p.value
}  

results<-data.frame(round(unlist(p_value_ass),3),unlist(list_gene))
logged_non_norm<-results[which(results$round.unlist.p_value_ass...3. < 0.05),] 
## after looging GLY, GT48, AAAP are still not normal... 

##un linear models ----
mods<-list()
summs<-list()
coeff<-list()
est<-list()
Stder<-list()
t_val<-list()
atable<-list()
p_value_ano<-list()
plots<-list()
for (i in 1:21){ ##MnP is 22 in the list so i should be able to just run loop only 21 times
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

ggarrange(plotlist =plots)
##not unreasonable

unlist(list_gene)[grep("singular",summs)]#LACC, OrgN_ACT, AAAP, GT48 

mod_resuls<-data.frame(unlist(list_gene)[-22],round(unlist(est),digits = 3),round(unlist(Stder),digits = 3),round(unlist(t_val),digits = 3),round(unlist(p_value_ano),digits = 4))
colnames(mod_resuls)<-c("gene","estimate","Std.error","t_value","p_value")

mod_resuls$adjusted_p<- p.adjust(mod_resuls$p_value,method = "fdr")

PROT_v_GT48_lm<-lmer(GT48~PROT+(1|Block),data= norm_wide)
summary(PROT_v_GT48_lm)$coefficients
Anova(PROT_v_GT48_lm)

GMC_1_v_GT48_lm<-lmer(GT48~GMC_1.1.3.13+(1|Block),data= norm_wide)
summary(GMC_1_v_GT48_lm)$coefficients
Anova(GMC_1_v_GT48_lm)
GMC_2_v_GT48_lm<-lmer(GT48~GMC_1.1.3.7+(1|Block),data= norm_wide)
summary(GMC_2_v_GT48_lm)$coefficients
Anova(GMC_2_v_GT48_lm)
GMC_3_v_GT48_lm<-lmer(GT48~GMC_1.1.99.18+(1|Block),data= norm_wide)
summary(GMC_3_v_GT48_lm)$coefficients
Anova(GMC_3_v_GT48_lm)

