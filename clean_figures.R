library(readr)
library(tidyverse)
library(patchwork)

#read in and clean ----
express_table<-read_csv("clean_data/gene_interest_express_GC.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

express_table<-express_table[-which(express_table$Block == "block19"),]

##reduce ----
##there is some Block/blocks where the difference is much stronger - lets limit to those where difference in KGD is greater than 3  
select<-which(express_table[which(express_table$gene == "MnP" & express_table$level == "high"),]$KGD-express_table[which(express_table$gene == "MnP" & express_table$level == "low"),]$KGD > 3)
block<-as.data.frame(express_table[which(express_table$gene == "MnP" & express_table$level == "high"),]$Block)
block[c(select),]

express_table.2<-express_table[which(express_table$Block %in% block[c(select),]),]

##wide dfs ----
norm_wide<-express_table[,-c(5,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)
norm_wide_sub<- express_table.2[,-c(5,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

##combine ---- 
norm_wide$AA<-rowSums(norm_wide[,which(colnames(norm_wide) %in% c("OrgN_AAAP","OrgN_ACT","OrgN_LAT","OrgN_YAT"))])
norm_wide_sub$AA<-rowSums(norm_wide_sub[,which(colnames(norm_wide_sub) %in% c("OrgN_AAAP","OrgN_ACT","OrgN_LAT","OrgN_YAT"))])
norm_wide$APET<-rowSums(norm_wide[,which(colnames(norm_wide) %in% c("OrgN_POT","OrgN_OPT"))])
norm_wide_sub$APET<-rowSums(norm_wide_sub[,which(colnames(norm_wide_sub) %in% c("OrgN_POT","OrgN_OPT"))])
norm_wide$GMC_sum<- rowSums(norm_wide[,grepl("GMC",colnames(norm_wide))])
norm_wide_sub$GMC_sum<- rowSums(norm_wide_sub[,grepl("GMC",colnames(norm_wide_sub))])
norm_wide$CHIT_sum<- rowSums(norm_wide[,grepl("CHIT",colnames(norm_wide))])
norm_wide_sub$CHIT_sum<- rowSums(norm_wide_sub[,grepl("CHIT",colnames(norm_wide_sub))])

##CUE ----

#sum_reads
#reads_wide<-express_table[,-c(4,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)
#reads_wide_sub<- express_table.2[,-c(4,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)
#make ratio
norm_wide$CUE<-norm_wide$GT48/norm_wide$KGD
norm_wide_sub$CUE<-norm_wide_sub$GT48/norm_wide_sub$KGD



#Paired plots ----

## MnP ----

MnP_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=MnP,fill = level))+
  geom_line(aes(group = Block,x=level,y=MnP))+
  geom_point(aes(fill = level,group = Block,x=level,y=MnP))+
  annotate("segment",x=1,xend = 2,y=30,yend =30,linewidth=1)+
  annotate("text",x=1.5,y=31,label="p = 0.004",size=5)+
  theme_classic()

MnP_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=MnP,fill = level))+
  geom_line(aes(group = Block,x=level,y=MnP))+
  geom_point(aes(fill = level,group = Block,x=level,y=MnP))+
  annotate("segment",x=1,xend = 2,y=30,yend =30,linewidth=1)+
  annotate("text",x=1.5,y=31,label="p = 0.002",size=5)+
  theme_classic()


MnP_all+MnP_sub+plot_annotation(tag_levels = "a")

## growth----

GT48_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=GT48,fill = level))+
  geom_line(aes(group = Block,x=level,y=GT48))+
  geom_point(aes(fill = level,group = Block,x=level,y=GT48))+
  annotate("segment",x=1,xend = 2,y=9,yend =9,linewidth=1)+
  annotate("text",x=1.5,y=9.2,label="*",size=5)+
  theme_classic()

CHSN_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=CHSN,fill = level))+
  geom_line(aes(group = Block,x=level,y=CHSN))+
  geom_point(aes(fill = level,group = Block,x=level,y=CHSN))+
  annotate("segment",x=1,xend = 2,y=17.5,yend =17.5,linewidth=1)+
  annotate("text",x=1.5,y=18,label="**",size=5)+
  theme_classic()

ALAS_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=ALAS,fill = level))+
  geom_line(aes(group = Block,x=level,y=ALAS))+
  geom_point(aes(fill = level,group = Block,x=level,y=ALAS))+
  annotate("segment",x=1,xend = 2,y=1.1,yend =1.1,linewidth=1)+
  annotate("text",x=1.5,y=1.15,label="**",size=5)+
  theme_classic()

CUE_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=CUE,fill = level))+
  geom_line(aes(group = Block,x=level,y=CUE))+
  geom_point(aes(fill = level,group = Block,x=level,y=CUE))+
  theme_classic()

GT48_sub+CUE_sub+CHSN_sub+ALAS_sub+plot_annotation(tag_levels = "a")


## H2O2----

GLY_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=GLY,fill = level))+
  geom_line(aes(group = Block,x=level,y=GLY))+
  geom_point(aes(fill = level,group = Block,x=level,y=GLY))+
  theme_classic()

GMC_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=GMC_sum,fill = level))+
  geom_line(aes(group = Block,x=level,y=GMC_sum))+
  geom_point(aes(fill = level,group = Block,x=level,y=GMC_sum))+
  theme_classic()

GLY_sub+GMC_sub+plot_annotation(tag_levels = "a")

## N mobilise and transport ----

PROT_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=PROT,fill = level))+
  geom_line(aes(group = Block,x=level,y=PROT))+
  geom_point(aes(fill = level,group = Block,x=level,y=PROT))+
  theme_classic()

APEP_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=APEP,fill = level))+
  geom_line(aes(group = Block,x=level,y=APEP))+
  geom_point(aes(fill = level,group = Block,x=level,y=APEP))+
  theme_classic()
AA_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=AA,fill = level))+
  geom_line(aes(group = Block,x=level,y=AA))+
  geom_point(aes(fill = level,group = Block,x=level,y=AA))+
  theme_classic()

APET_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=APET,fill = level))+
  geom_line(aes(group = Block,x=level,y=APET))+
  geom_point(aes(fill = level,group = Block,x=level,y=APET))+
  theme_classic()


PROT_sub+AA_sub+APEP_sub+APET_sub+plot_layout(nrow=2)+plot_annotation(tag_levels = "a")


## chitin----

CHIT_inase_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=CHIT_3.2.1.14,fill = level))+
  geom_line(aes(group = Block,x=level,y=CHIT_3.2.1.14))+
  geom_point(aes(fill = level,group = Block,x=level,y=CHIT_3.2.1.14))+
  theme_classic()

CHIT_deace_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=CHIT_3.5.1.41,fill = level))+
  geom_line(aes(group = Block,x=level,y=CHIT_3.5.1.41))+
  geom_point(aes(fill = level,group = Block,x=level,y=CHIT_3.5.1.41))+
  theme_classic()

CHIT_NAG_sub<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=CHIT_3.2.1.52,fill = level))+
  geom_line(aes(group = Block,x=level,y=CHIT_3.2.1.52))+
  geom_point(aes(fill = level,group = Block,x=level,y=CHIT_3.2.1.52))+
  theme_classic()

CHIT_inase_sub+CHIT_deace_sub+CHIT_NAG_sub+plot_annotation(tag_levels = "a")

####sup 1 high v low all ----

GT48_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=GT48,fill = level))+
  geom_line(aes(group = Block,x=level,y=GT48))+
  geom_point(aes(fill = level,group = Block,x=level,y=GT48))+
  theme_classic()

CHSN_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=CHSN,fill = level))+
  geom_line(aes(group = Block,x=level,y=CHSN))+
  geom_point(aes(fill = level,group = Block,x=level,y=CHSN))+
  theme_classic()

GLY_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=GLY,fill = level))+
  geom_line(aes(group = Block,x=level,y=GLY))+
  geom_point(aes(fill = level,group = Block,x=level,y=GLY))+
  theme_classic()

ALAS_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=ALAS,fill = level))+
  geom_line(aes(group = Block,x=level,y=ALAS))+
  geom_point(aes(fill = level,group = Block,x=level,y=ALAS))+
  theme_classic()

LACC_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=LACC,fill = level))+
  geom_line(aes(group = Block,x=level,y=LACC))+
  geom_point(aes(fill = level,group = Block,x=level,y=LACC))+
  theme_classic()

PROT_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=PROT,fill = level))+
  geom_line(aes(group = Block,x=level,y=PROT))+
  geom_point(aes(fill = level,group = Block,x=level,y=PROT))+
  theme_classic()

APEP_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=APEP,fill = level))+
  geom_line(aes(group = Block,x=level,y=APEP))+
  geom_point(aes(fill = level,group = Block,x=level,y=APEP))+
  theme_classic()

GMC_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=GMC_sum,fill = level))+
  geom_line(aes(group = Block,x=level,y=GMC_sum))+
  geom_point(aes(fill = level,group = Block,x=level,y=GMC_sum))+
  theme_classic()

CHIT_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=CHIT_sum,fill = level))+
  geom_line(aes(group = Block,x=level,y=CHIT_sum))+
  geom_point(aes(fill = level,group = Block,x=level,y=CHIT_sum))+
  theme_classic()

AA_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=AA,fill = level))+
  geom_line(aes(group = Block,x=level,y=AA))+
  geom_point(aes(fill = level,group = Block,x=level,y=AA))+
  theme_classic()

APET_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=APET,fill = level))+
  geom_line(aes(group = Block,x=level,y=APET))+
  geom_point(aes(fill = level,group = Block,x=level,y=APET))+
  theme_classic()

GT48_all+CHSN_all+ALAS_all+LACC_all+GLY_all+GMC_all+PROT_all+APEP_all+AA_all+APET_all+CHIT_all+plot_annotation(tag_levels = "a")



#Scatter plots ----

#assumptions
list_gene<-list(c(colnames(norm_wide[,-c(1:3,7,8,12,13,14,29:32)])))
results<-list()
p_value_ass<-list()

for (i in 1:20){
  results[[i]]<-shapiro.test(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene[[1]][[i]])]))
  p_value_ass[[i]]<-results[[i]]$p.value
}  

results<-data.frame(round(unlist(p_value_ass),3),unlist(list_gene))
non_norm<-results[which(results$round.unlist.p_value_ass...3. < 0.05),] #OrgN POT is barely expressed so I will exclude

##check non_normal distributions
hist(log10(norm_wide$ALAS))
hist(norm_wide$GMC_1.1.3.13)#seems fine
hist(log10(norm_wide$GLUC))
hist(log10(norm_wide$LACC))
hist(log10(norm_wide$OrgN_LAT))
hist(log10(norm_wide$OrgN_OPT))
hist(log10(norm_wide$CHIT_3.5.1.41))

norm_wide$ALAS_log<-log10(norm_wide$ALAS)
norm_wide$LACC_log<-log10(norm_wide$LACC)
norm_wide$OrgN_LAT_log<-log10(norm_wide$OrgN_LAT)
norm_wide$OrgN_OPT_log<-log10(norm_wide$OrgN_OPT)
norm_wide$CHIT_deace_log<-log10(norm_wide$CHIT_3.5.1.41)
norm_wide$GLUC_log<-log10(norm_wide$GLUC)

list_gene_norm<-list(c(colnames(norm_wide[,c(33:38)])))
results<-list()
p_value_ass<-list()

for (i in 1:6){
  results[[i]]<-shapiro.test(unlist(norm_wide[,which(colnames(norm_wide) %in% list_gene_norm[[1]][[i]])]))
  p_value_ass[[i]]<-results[[i]]$p.value
}  

results<-data.frame(round(unlist(p_value_ass),3),unlist(list_gene_norm))
logged_non_norm<-results[which(results$round.unlist.p_value_ass...3. < 0.05),] #OrgN POT is barely expressed so I will exclude
#fixed

##growth/CUE vs MnP ----

GT48_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(GT48)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "GT48 / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

KGD_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(KGD)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "KGD (log-transformed)")+
  xlab(label = "sqrt(MnP)")+
  theme_classic()

CUE_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=CUE))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "CUE = GT48/KGD")+
  theme_classic()


##CUE vs others ----

CUE_v_PROT_sub<-ggplot(norm_wide,aes(y=CUE,x=PROT))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "CUE / KGD ")+
  xlab(label = "A1 Proteases / KGD")+
  theme_classic()

CUE_v_GMC_1_sub<-ggplot(norm_wide,aes(y=KGD,x=GMC_1.1.3.13))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "CUE / KGD ")+
  xlab(label = "Alcohol oxidase / KGD")+
  theme_classic()

CUE_v_GMC_2_sub<-ggplot(norm_wide,aes(y=KGD,x=GMC_1.1.3.7))+
  geom_point()+
  geom_smooth(method = "lm")+
  #ylab(label = "CUE / KGD ")+
  #xlab(label = "Aryl-alcohol oxidase / KGD")+
  theme_classic()

CUE_v_GMC_3_sub<-ggplot(norm_wide,aes(y=CUE,x=GMC_1.1.99.18))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "CUE / KGD ")+
  xlab(label = "Cellobiose dehydrogenase / KGD")+
  theme_classic()


##chitin break down vs MnP ----

CHIT_deace_reg<-ggplot(norm_wide,aes(x=MnP,y=CHIT_3.5.1.41))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Chitin deacetylase / KGD (log-transformed)")+
  xlab(label = "MnP / KGD")+
  theme_classic()

CHIT_inase_reg<-ggplot(norm_wide,aes(x=MnP,y=CHIT_3.2.1.14))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Chitinase / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

CHIT_NAG_reg<-ggplot(norm_wide,aes(x=MnP,y=CHIT_3.2.1.52))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "NAG / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()



##N breakdown and transport vs mnp ----
PROT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(PROT)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "A1 Proteases (log10-transformed)")+
  xlab(label = "sqrt(MnP)")+
  theme_classic()

APEP_reg<-ggplot(norm_wide,aes(x=MnP,y=APEP))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Aminopeptidases / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

OrgN_AAAP_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(OrgN_AAAP)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "AAAP amino acid transport (2.A.18) (log10-transformed) ")+
  xlab(label = "sqrt(MnP)")+
  theme_classic()

OrgN_ACT_reg<-ggplot(norm_wide,aes(x=MnP,y=OrgN_ACT))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

OrgN_LAT_reg<-ggplot(norm_wide,aes(x=MnP,y=OrgN_LAT))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  theme_classic()

OrgN_YAT_reg<-ggplot(norm_wide,aes(x=MnP,y=OrgN_YAT))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  theme_classic()

OrgN_OPT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(OrgN_OPT)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "OPT oligopeptide transport (2.A.67) (log10-transformed) ")+
  xlab(label = "sqrt(MnP)")+
  theme_classic()


##laccases vs mnp ----
LACC_reg<-ggplot(norm_wide,aes(x=MnP,y=LACC))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  theme_classic()

##Mnp vs heme ----
ALAS_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log(ALAS)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ALAS / KGD (log-transformed)")+
  xlab(label = "MnP / KGD")+
  theme_classic()

##GMC vs MnP ----
GMC_1_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(GMC_1.1.3.13)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Alcohol oxidase (log-transformed) ")+
  xlab(label = "sqrt(MnP) ")+
  theme_classic()

GMC_2_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(GMC_1.1.3.7)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Aryl-alcohol oxidase (log-transformed) ")+
  xlab(label = "sqrt(MnP)")+
  theme_classic()

GMC_3_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(GMC_1.1.99.18)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Cellobiose dehydrogenase / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

##CRO vs MnP ----
GLY_reg<-ggplot(norm_wide,aes(x=MnP,y=GLY))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Copper radical oxidases / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

##C breakdown vs MnP ----
BETA_reg<-ggplot(norm_wide,aes(x=MnP,y=BETA))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Beta-glucosidase (EC 3.2.1.21) / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

GLUC_reg<-ggplot(norm_wide,aes(x=MnP,y=GLUC))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Glucanases (EC 3.2.1.58) / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()


#figure 1 ----

tiff(filename = "figures/figure_1.tiff")
GT48_reg

#figure 2 ----
tiff(filename = "figures/figure_2.tiff")
PROT_reg+APEP_reg+OrgN_AAAP_reg+OrgN_OPT_reg+plot_annotation(tag_levels = "a")


#figure 3 ----

tiff(filename = "figures/figure_3.tiff")
GMC_1_reg+GMC_2_reg+GMC_3_reg+GLY_reg+plot_annotation(tag_levels = "a")


#figures sub ----

tiff(filename = "figures/figures_s1.tiff")
BETA_reg+GLUC_reg+plot_annotation(tag_levels = "a")

tiff(filename = "figures/figures_s2.tiff")
GT48_v_PROT_sub|(GT48_v_GMC_1_sub/GT48_v_GMC_2_sub/GT48_v_GMC_3_sub)+plot_annotation(tag_levels = "a")

tiff(filename = "figures/figures_s3.tiff")
ALAS_reg

tiff(filename = "figures/figures_s4.tiff")
CHIT_deace_reg+CHIT_inase_reg+CHIT_NAG_reg+plot_annotation(tag_levels = "a")



write_csv(norm_wide,"clean_data/norm_wide_clean.csv")
write_csv(norm_wide_sub,"clean_data/norm_wide_sub_clean.csv")
