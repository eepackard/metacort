library(readr)
library(tidyverse)
library(patchwork)

#read in and clean ----
express_table_s<-read_csv("clean_data/gene_interest_express_spruce_GC.csv")
express_table_p<-read_csv("clean_data/gene_interest_express_pine_GC.csv")

express_table_p$species<-rep("pine",nrow(express_table_p))
express_table_s$species<-rep("spruce",nrow(express_table_s))

express_table <- rbind(express_table_p,express_table_s)

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)
express_table$species <-as.factor(express_table$species)

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

#bring in blast data
blast_results <- read_csv("clean_data/blast_1000_results_clean.csv")

#add to pine just the need variables
#first check same order
norm_wide_Coromn$levelblock == paste(norm_wide$level,norm_wide$Block,sep = "")
blast_results$levelblock == paste(norm_wide$level,norm_wide$Block,sep = "")
#need to repeat mnP twice (because data dupped for pine and spruce)

norm_wide$MnP <- c(norm_wide_Coromn$MnP,norm_wide_Coromn$MnP)
norm_wide$TASE <- c(norm_wide_Coromn$TASE,norm_wide_Coromn$TASE)
norm_wide$TPP <- c(norm_wide_Coromn$TPP,norm_wide_Coromn$TPP)
norm_wide$TPS <- c(norm_wide_Coromn$TPS,norm_wide_Coromn$TPS)
norm_wide$prop_cort <- c(blast_results$prop_cort_fun,blast_results$prop_cort_fun)

#normalisation
norm_wide$sweet_norm <- norm_wide$Sweet*norm_wide$prop_cort
norm_wide$invert_norm <- norm_wide$invert*norm_wide$prop_cort
norm_wide$tub_norm <- norm_wide$tub*norm_wide$prop_cort
norm_wide$actin_norm <- norm_wide$actin*norm_wide$prop_cort

##combine actin and tubulin

plot(norm_wide$tub,norm_wide$actin) #10 fold difference in tpm
plot(norm_wide$SWEET_per_tub,norm_wide$SWEET_per_actin) ## but good consistency between normalization methods - except one outlier...
plot(norm_wide$invert_per_tub,norm_wide$invert_per_actin) ## but good consistency between normalization methods

plot(norm_wide$Sweet,norm_wide$sweet_norm)

##compare pine and spruce
plot(norm_wide[which(norm_wide$species == "pine"),]$tub,norm_wide[which(norm_wide$species == "spruce"),]$tub) 
plot(norm_wide[which(norm_wide$species == "pine"),]$actin,norm_wide[which(norm_wide$species == "spruce"),]$actin) 
plot(norm_wide[which(norm_wide$species == "pine"),]$Sweet,norm_wide[which(norm_wide$species == "spruce"),]$Sweet) 
plot(norm_wide[which(norm_wide$species == "pine"),]$invert,norm_wide[which(norm_wide$species == "spruce"),]$invert) 
plot(norm_wide[which(norm_wide$species == "pine"),]$SWEET_per_actin,norm_wide[which(norm_wide$species == "spruce"),]$SWEET_per_actin) 
plot(norm_wide[which(norm_wide$species == "pine"),]$invert_per_actin,norm_wide[which(norm_wide$species == "spruce"),]$invert_per_actin) 



#Paired plots ----

## refs ----

tub<-ggplot(norm_wide,aes(x=level,y=tub))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)

tub<-ggplot(norm_wide,aes(x=level,y=tub_norm))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)

actin<-ggplot(norm_wide,aes(x=level,y=actin))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)

ggplot(norm_wide,aes(x=level,y=actin_norm))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)


## sweet----

sweet<-ggplot(norm_wide,aes(x=level,y=Sweet))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)

sweet<-ggplot(norm_wide,aes(x=level,y=sweet_norm))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)

## invert ----

invert<-ggplot(norm_wide,aes(x=level,y=invert))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)

invert<-ggplot(norm_wide,aes(x=level,y=invert_norm))+
  geom_boxplot(aes(fill=level))+
  geom_line(aes(group = Block))+
  geom_point()+
  theme_classic()+
  facet_wrap(~species)


#Scatter plots ----
#test assmuptions in LM script

##ref vs MnP ----

tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=tub,colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "TUB")+
  xlab(label = "MnP")+
  theme_classic()

tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=tub_norm,colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "TUB")+
  xlab(label = "MnP")+
  theme_classic()

act_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ACTIN")+
  xlab(label = "MnP")+
  theme_classic()

ggplot(norm_wide,aes(x=sqrt(MnP),y=actin_norm,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ACTIN")+
  xlab(label = "MnP")+
  theme_classic()

ggplot(norm_wide,aes(x=TASE/TPS,y=actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ACTIN")+
  xlab(label = "TASE/TPS")+
  theme_classic()

ggplot(norm_wide,aes(x=TASE/TPS,y=actin_norm,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "ACTIN")+
  xlab(label = "TASE/TPS")+
  theme_classic()


##sweet vs MnP ----

hist(norm_wide$Sweet)
hist(norm_wide$sweet_norm)
hist(log10(norm_wide$sweet_norm))
hist(norm_wide$SWEET_per_tub)
hist(norm_wide$SWEET_per_actin)

sweet_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(Sweet),colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "SWEETs")+
  xlab(label = "MnP")+
  theme_classic()

sweet_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=log10(sweet_norm),colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "SWEETs")+
  xlab(label = "MnP")+
  theme_classic()

SWE_tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=SWEET_per_tub,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "SWEET / TUB")+
  xlab(label = "MnP")+
  theme_classic()

ggplot(norm_wide,aes(x=sqrt(MnP),y=SWEET_per_tub*prop_cort,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "SWEET / TUB")+
  xlab(label = "MnP")+
  theme_classic()

SWE_act_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=SWEET_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Pine SWEETs / Actin")+
  xlab(label = "MnP")+
  theme_classic()

##invert vs MnP ----

hist(norm_wide$invert)
hist(norm_wide$invert_per_tub)
hist(norm_wide$invert_per_actin)

invert_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=invert,colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Pine invertase")+
  xlab(label = "MnP")+
  theme_classic()

inv_tub_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=invert_per_tub,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Pine invertase / TUB")+
  xlab(label = "MnP")+
  theme_classic()

inv_act_reg<-ggplot(norm_wide,aes(x=sqrt(MnP),y=invert_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "Pine invertase / Actin")+
  xlab(label = "MnP")+
  theme_classic()

##save
tiff(filename = "figures/figure_NPS_pres_spruce.tiff",height = 1000,width = 2500,units = "px",res = 300)
inv_act_reg+SWE_act_reg

##invert vs sweet ----

inv_swe_reg<-ggplot(norm_wide,aes(x=Sweet,y=invert,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet")+
  xlab(label = "invert")+
  theme_classic()

ggplot(norm_wide,aes(x=sweet_norm,y=invert_norm,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet")+
  xlab(label = "invert")+
  theme_classic()

inv_swe_reg<-ggplot(norm_wide,aes(x=SWEET_per_actin,y=invert_per_actin,colour=species,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet /actin")+
  xlab(label = "invert /actin")+
  theme_classic()

inv_swe_reg<-ggplot(norm_wide,aes(x=SWEET_per_tub,y=invert_per_tub,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet")+
  xlab(label = "invert")+
  theme_classic()

## trehalose vs sweet ----
TASE_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=SWEET_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet /actin")+
  xlab(label = "Trehalase (sqrt)")+
  theme_classic()

TPP_reg<-ggplot(norm_wide,aes(x=sqrt(TPP),y=SWEET_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Sweet/ actin")+
  xlab(label = "TPP (sqrt)")+
  theme_classic()

TPS_reg<-ggplot(norm_wide,aes(x=sqrt(TPS),y=SWEET_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "SWEET /actin")+
  xlab(label = "TPS (sqrt)")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPP,y=SWEET_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet / actin")+
  xlab(label = "Trehalase/TPP")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPS,y=Sweet,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "sweet/actin")+
  xlab(label = "Trehalase/TPS ")+
  theme_classic()

TASE_reg<-ggplot(norm_wide,aes(x=sqrt(TASE),y=invert,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "Trehalase (sqrt)")+
  xlab(label = "invert")+
  theme_classic()

TPP_reg<-ggplot(norm_wide,aes(x=sqrt(TPP),y=invert_per_actin,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert / actin")+
  xlab(label = "TPP (sqrt) ")+
  theme_classic()

TPS_reg<-ggplot(norm_wide,aes(x=sqrt(TPS),y=invert,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert  / actin ")+
  xlab(label = "TPS (sqrt)")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPP,y=invert,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert  / actin")+
  xlab(label = "Trehalase/TPP")+
  theme_classic()

TASE/TPP_reg<-ggplot(norm_wide,aes(x=TASE/TPS,y=invert,colour=species))+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab(label = "invert  / actin")+
  xlab(label = "Trehalase/TPS")+
  theme_classic()


