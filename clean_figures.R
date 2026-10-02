library(readr)
library(tidyverse)
library(patchwork)

#read in and clean ----
norm_wide<-read_csv("clean_data/norm_wide_clean.csv")
norm_wide_sub<-read_csv("clean_data/norm_wide_sub_clean.csv")

#Paired plots ----

## MnP ----

MnP_all<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=MnP,fill = level))+
  geom_line(aes(group = Block,x=level,y=MnP))+
  geom_point(aes(fill = level,group = Block,x=level,y=MnP))+
  annotate("segment",x=1,xend = 2,y=30,yend =30,linewidth=1)+
  annotate("text",x=1.5,y=31,label="p = 0.002",size=5)+
  theme_classic()

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

CHIT_NAGt<-ggplot(norm_wide_sub)+
  geom_boxplot(aes(x=level,y=NAGt,fill = level))+
  geom_line(aes(group = Block,x=level,y=NAGt))+
  geom_point(aes(fill = level,group = Block,x=level,y=NAGt))+
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
#test assmuptions in LM script

##growth/CUE vs MnP ----

GT48_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GT48)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "GT48 / KGD ")+
  xlab(label = "MnP / KGD")+
  theme_classic()

KGD_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(KGD)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "KGD")+
  xlab(label = "MnP")+
  theme_classic()

CUE_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=CUE))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "CUE = GT48/KGD")+
  xlab(label = "MnP")+
  theme_classic()


##KGD vs others ----

CUE_v_PROT_sub<-ggplot(norm_wide,aes(y=sqrt(KGD),x=sqrt(PROT)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "KGD ")+
  xlab(label = "A1 Proteases")+
  theme_classic()

CUE_v_GMC_1_sub<-ggplot(norm_wide,aes(y=sqrt(KGD),x=sqrt(GMC_AOx)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "KGD ")+
  xlab(label = "Alcohol oxidases")+
  theme_classic()

CUE_v_GMC_2_sub<-ggplot(norm_wide,aes(y=sqrt(KGD),x=sqrt(GMC_AAO)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "KGD ")+
  xlab(label = "Aryl-alcohol oxidases")+
  theme_classic()

CUE_v_GMC_2_sub<-ggplot(norm_wide,aes(y=sqrt(KGD),x=sqrt(GMC_GDH)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "KGD ")+
  xlab(label = "GDH")+
  theme_classic()

CUE_v_ALAS<-ggplot(norm_wide,aes(y=sqrt(KGD),x=sqrt(ALAS)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "KGD")+
  xlab(label = "ALAS")+
  theme_classic()

##chitin break down vs MnP ----

CHIT_deace_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(CHIT_3.5.1.41)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Chitin deacetylases")+
  xlab(label = "")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

CHIT_inase_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(CHIT_3.2.1.14)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Chitinases")+
  xlab(label = "MnP")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))

CHIT_NAG_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(CHIT_3.2.1.52)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "NAG")+
  xlab(label = "MnP")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))

CHIT_NAGt_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(NAGt)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "NAGt")+
  xlab(label = "MnP")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))

CHIT_NAGt_reg<-ggplot(norm_wide,aes(x=sqrt(CHIT_3.2.1.52),y=sqrt(NAGt)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "NAGt")+
  xlab(label = "NAG")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))


##N breakdown and transport vs mnp ----
#tiff(filename = "figures/asp_prot_pres.tiff",height = 1000,width = 1500,units = "px",res = 300)
PROT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(PROT)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Aspartic Proteases")+
  xlab(label = "")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

APEP_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(APEP)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Aminopeptidases")+
  xlab(label = "")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_AAAP_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_AAAP)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Amino acid/auxin \n porters")+
  xlab(label = "")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_ACT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_ACT)))+
  geom_point()+
  #geom_smooth(method = "lm")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_LAT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_LAT)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "L-type amino acid \n transporters")+
  xlab(label = "")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_YAT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_YAT)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_AAT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_AAT)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_OPT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_OPT)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Oligopeptide \n transporters")+
  xlab(label = "")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 16))

OrgN_POT_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(OrgN_POT)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "POT")+
  xlab(label = "MnP")+
  theme_classic()


##laccases vs mnp ----
LACC_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(LACC)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  theme_classic()

##Mnp vs heme ----
ALAS_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(ALAS)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ALAS")+
  xlab(label = "MnP")+
  theme_classic()

##GMC vs MnP ----
GMC_1_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GMC_AOx)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "GMC AA3_3 (AOx)")+
  xlab(label = "MnP")+
  theme_classic()

GMC_2_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GMC_AAO)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "GMC AA3_2 (AAO/PDH-like)")+
  xlab(label = "MnP")+
  theme_classic()

GMC_3_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GMC_GDH)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "GMC AA3_2 (GDH-like)")+
  xlab(label = "MnP")+
  theme_classic()

## trehalose vs MnP ----
TASE_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(TASE)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Trehalase")+
  xlab(label = "MnP")+
  theme_classic()

TPP_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(TPP)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "TPP")+
  xlab(label = "MnP")+
  theme_classic()

TPS_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(TPS)))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "TPS")+
  xlab(label = "MnP")+
  theme_classic()

TASE_TPP_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=TASE_TPP))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Trehalase/TPP")+
  xlab(label = "MnP")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=sqrt(KGD),y=TASE_TPP))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Trehalase/TPP")+
  xlab(label = "KGD")+
  theme_classic()

TASE_vs_AAO_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=sqrt(GMC_AAO)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "GMC AA3_2 (AAO/PDH-like)")+
  xlab(label = "Trehalose")+
  theme_classic()

TASE_vs_GDH_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=sqrt(GMC_GDH)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "GMC AA3_2 (GDH-like)")+
  xlab(label = "Trehalose")+
  theme_classic()

TASE_vs_AoX_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=sqrt(GMC_AOx)))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "GMC AA3_3")+
  xlab(label = "Trehalose")+
  theme_classic()

TASE/TPS_vs_AAO_reg<-ggplot(norm_wide,aes(x=sqrt(GMC_AAO),y=TASE_TPP))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Trehalase/TPP")+
  xlab(label = "GMC AA3_2 (AAO/PDH-like)")+
  theme_classic()

##CRO vs MnP ----
GLY_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GLY)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Copper radical oxidases")+
  xlab(label = "MnP")+
  theme_classic()

##C breakdown vs MnP ----
BETA_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(BETA)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = expression(""*beta*"-glucosidases (EC 3.2.1.21)"))+
  xlab(label = "MnP")+
  theme_classic()

GLUC_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=sqrt(GLUC)))+
  geom_point()+
  #geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = expression("Glucan1,3-"*beta*"-glucosidases (EC 3.2.1.58)"))+
  xlab(label = "MnP")+
  theme_classic()


## TCA vs TASE/TPS/TPS



#figure 2 ----

tiff(filename = "figures/figure_2.tiff",height = 1000,width = 2500,units = "px",res = 300)
CUE_reg+KGD_reg+plot_annotation(tag_levels = "a")
#CUE_reg+ inset_element(KGD_reg, 0.6,0.01,1,0.4)

#figure x ----
tiff(filename = "figures/figure_NPS_pres_TPS_TASE.tiff",height = 1000,width = 1250,units = "px",res = 300)
TASE_TPP_reg
PROT_reg+OrgN_OPT_reg
PROT_reg+OrgN_OPT_reg+APEP_reg+OrgN_AAAP_reg+plot_annotation(tag_levels = "a")

#figure 3 ----

tiff(filename = "figures/figure_3.tiff",height = 2000,width = 2500,units = "px",res = 300)
GMC_2_reg+GMC_1_reg+GMC_3_reg+GLY_reg+plot_annotation(tag_levels = "a")+plot_layout(nrow = 2)

#figure x ----
tiff(filename = "figures/figure_x2.tiff")
CHIT_deace_reg+CHIT_inase_reg+CHIT_NAG_reg+plot_annotation(tag_levels = "a")+plot_layout(nrow = 2)

#figure 4 ----
tiff(filename = "figures/figure_4.tiff", width = 3000,height = 4000, units = "px",res = 300)
PROT_reg+OrgN_OPT_reg+OrgN_AAAP_reg+OrgN_LAT_reg+APEP_reg+
  CHIT_deace_reg+CHIT_inase_reg+CHIT_NAG_reg+plot_annotation(tag_levels = "a")+plot_layout(ncol = 2)

#figures sub ----

tiff(filename = "figures/figures_s1.tiff",height = 1000,width = 3000,units = "px",res = 300)
BETA_reg+GLUC_reg+plot_annotation(tag_levels = "a")

tiff(filename = "figures/figures_s2.tiff")
ALAS_reg
