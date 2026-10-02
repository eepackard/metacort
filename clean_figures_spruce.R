library(readr)
library(tidyverse)
library(patchwork)

#read in and clean ----
express_table<-read_csv("clean_data/gene_interest_express_spruce_GC.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

express_table<-express_table[-which(express_table$Block == "block19"),]

##wide dfs ----
norm_wide<-express_table[,-c(5)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

##norm ---- 
#make ratio
norm_wide$SWEET_per_tub<-norm_wide$Sweet/norm_wide$tub
norm_wide$SWEET_per_actin<-norm_wide$Sweet/norm_wide$actin
norm_wide$invert_per_tub<-norm_wide$invert/norm_wide$tub
norm_wide$invert_per_actin<-norm_wide$invert/norm_wide$actin


#write_csv(norm_wide,"clean_data/norm_wide_clean_pine.csv")

#bring in Cort data
norm_wide_Coromn<-read_csv("clean_data/norm_wide_clean.csv")

#add to pine just the need variables
#first check same order
paste(norm_wide_Coromn$level,norm_wide_Coromn$Block) == paste(norm_wide$level,norm_wide$Block)

norm_wide$MnP <- norm_wide_Coromn$MnP
norm_wide$TASE <- norm_wide_Coromn$TASE
norm_wide$TPP <- norm_wide_Coromn$TPP
norm_wide$TPS <- norm_wide_Coromn$TPS

##combine actin and tubulin

plot(norm_wide$tub,norm_wide$actin) #10 fold difference in tpm
plot(norm_wide$SWEET_per_tub,norm_wide$SWEET_per_actin) ## but good consistency between normalization methods - except one outlier...
plot(norm_wide$invert_per_tub,norm_wide$invert_per_actin) ## but good consistency between normalization methods 


#Paired plots ----

## refs ----

tub<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=tub,fill = level))+
  geom_line(aes(group = Block,x=level,y=tub))+
  geom_point(aes(fill = level,group = Block,x=level,y=tub))+
  theme_classic()

actin<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=actin,fill = level))+
  geom_line(aes(group = Block,x=level,y=actin))+
  geom_point(aes(fill = level,group = Block,x=level,y=actin))+
  theme_classic()


## sweet----

sweet<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=Sweet,fill = level))+
  geom_line(aes(group = Block,x=level,y=Sweet))+
  geom_point(aes(fill = level,group = Block,x=level,y=Sweet))+
  theme_classic()

## invert ----

invert<-ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=invert,fill = level))+
  geom_line(aes(group = Block,x=level,y=invert))+
  geom_point(aes(fill = level,group = Block,x=level,y=invert))+
  theme_classic()


#Scatter plots ----
#test assmuptions in LM script

##ref vs MnP ----

tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=tub))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "TUB")+
  xlab(label = "MnP")+
  theme_classic()

act_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ACTIN")+
  xlab(label = "MnP")+
  theme_classic()

ggplot(norm_wide,aes(x=TASE/TPS,y=actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ACTIN")+
  xlab(label = "TASE/TPS")+
  theme_classic()

##sweet vs MnP ----

hist(norm_wide$Sweet)
hist(norm_wide$SWEET_per_tub)
hist(norm_wide$SWEET_per_actin)

sweet_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=Sweet))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "SWEETs")+
  xlab(label = "MnP")+
  theme_classic()

SWE_tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=SWEET_per_tub))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "SWEET / TUB")+
  xlab(label = "MnP")+
  theme_classic()

SWE_act_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=SWEET_per_actin))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Pine SWEETs / Actin")+
  xlab(label = "MnP")+
  theme_classic()

##invert vs MnP ----

hist(norm_wide$invert)
hist(norm_wide$invert_per_tub)
hist(norm_wide$invert_per_actin)

invert_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=invert))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Pine invertase")+
  xlab(label = "MnP")+
  theme_classic()

inv_tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=invert_per_tub))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Pine invertase / TUB")+
  xlab(label = "MnP")+
  theme_classic()

inv_act_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=invert_per_actin))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Pine invertase / Actin")+
  xlab(label = "MnP")+
  theme_classic()

##save
tiff(filename = "figures/figure_NPS_pres_spruce.tiff",height = 1000,width = 2500,units = "px",res = 300)
inv_act_reg+SWE_act_reg

##invert vs sweet ----

inv_swe_reg<-ggplot(norm_wide,aes(x=Sweet,y=invert))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet")+
  xlab(label = "invert")+
  theme_classic()

inv_swe_reg<-ggplot(norm_wide,aes(x=SWEET_per_actin,y=invert_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet /actin")+
  xlab(label = "invert /actin")+
  theme_classic()

inv_swe_reg<-ggplot(norm_wide,aes(x=SWEET_per_tub,y=invert_per_tub))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet")+
  xlab(label = "invert")+
  theme_classic()

## trehalose vs sweet ----
TASE_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=SWEET_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet /actin")+
  xlab(label = "Trehalase (sqrt)")+
  theme_classic()

TPP_reg<-ggplot(norm_wide,aes(x=sqrt(TPP),y=SWEET_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Sweet/ actin")+
  xlab(label = "TPP (sqrt)")+
  theme_classic()

TPP/TPS_ave_reg<-ggplot(norm_wide,aes(x=(sqrt(TPP)+sqrt(TPS)/2),y=SWEET_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Sweet/ actin")+
  xlab(label = "TPP (sqrt)")+
  theme_classic()

TPS_reg<-ggplot(norm_wide,aes(x=sqrt(TPS),y=SWEET_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "SWEET /actin")+
  xlab(label = "TPS (sqrt)")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPP,y=SWEET_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet / actin")+
  xlab(label = "Trehalase/TPP")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPS,y=Sweet))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet/actin")+
  xlab(label = "Trehalase/TPS ")+
  theme_classic()

TASE_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=invert))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Trehalase (sqrt)")+
  xlab(label = "invert")+
  theme_classic()

TPP_reg<-ggplot(norm_wide,aes(x=sqrt(TPP),y=invert_per_actin))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert / actin")+
  xlab(label = "TPP (sqrt) ")+
  theme_classic()

TPS_reg<-ggplot(norm_wide,aes(x=sqrt(TPS),y=invert))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert  / actin ")+
  xlab(label = "TPS (sqrt)")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPP,y=invert))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert  / actin")+
  xlab(label = "Trehalase/TPP")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPS,y=invert))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert  / actin")+
  xlab(label = "Trehalase/TPS")+
  theme_classic()


