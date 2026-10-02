library(readr)
library(tidyverse)


#read in and clean
matches <- read_delim("clean_data/count_matched.csv",delim = ";")
matches <- matches[,-1]
HMM_results<- read_csv("clean_data/HMM_results.csv")


matches$block <- paste("block",matches$block, sep = "")
matches$levelblock <- paste(matches$level,matches$block,sep = "")

matches <- matches[match(HMM_results$levelblock,matches$levelblock),] 

#tabulate

matches$prop_omn_pil <- matches$matches_omn_pil/matches$coromn
matches$prop_omn_ful <- matches$matches_omn_ful/matches$coromn
matches$prop_omn_pur <- matches$matches_omn_pur/matches$coromn

matches_long <- pivot_longer(matches[,10:13],cols = 2:4,names_to = "species")

ggplot(matches_long)+
  geom_boxplot(aes(x=species,y=value),)+
  theme_classic()


#read in percent mapping

mapping <- read_delim("clean_data/Mapping_rate.csv",delim = ";")
mapping$mapping_rate <- mapping$mapping_rate*100 #change prop to percent
mapping_coromn <- mapping[which(mapping$sp == "coromn"),]
mapping_coromn$levelblock <- paste(mapping_coromn$level,mapping_coromn$block,sep = "block")
mapping_others <- mapping[-which(mapping$sp == "coromn"),]
mapping_others$levelblock <- paste(mapping_others$level,mapping_others$block,sep = "block")

mapping_others$corom <- mapping_coromn[match(mapping_others$levelblock,mapping_coromn$levelblock),]$mapping_rate

ggplot(mapping_others)+
  geom_point(aes(x=corom,y=mapping_rate,colour = sp,shape = level,size=2))+
  geom_abline(intercept=0,slope=1,col="black")+
  geom_smooth(aes(x=corom,y=mapping_rate,colour = sp),method = "lm")+
  theme_classic()

lmpur<-lm(mapping_rate~corom,data=mapping_others[which(mapping_others$sp == "corpur"),])
lmful<-lm(mapping_rate~corom,data=mapping_others[which(mapping_others$sp == "corful"),])
lmpil<-lm(mapping_rate~corom,data=mapping_others[which(mapping_others$sp == "pilsph"),])

(1-lmpur$coefficients[2])/1 #21 % decrease from coromn (i.e. perfect 1:1 match)
(lmpur$coefficients[2]-lmful$coefficients[2])/lmpur$coefficients[2] #70% decrease from pur to ful
(1-lmful$coefficients[2])/1 #76 % decrease from coromn (i.e. perfect 1:1 match)
(lmpur$coefficients[2]-lmpil$coefficients[2])/lmpur$coefficients[2] #100% decrase - this is possible because slope of pil is essentiall zero (neg)


matches_long$sp <- if_else(grepl("pil",matches_long$species),"pilsph",if_else(grepl("ful",matches_long$species),"corful",if_else(grepl("pur",matches_long$species),"corpur",NA)))
matches_long$levelblocksp <- paste(matches_long$levelblock,matches_long$sp,sep = "")
mapping_others$levelblocksp <- paste(mapping_others$levelblock,mapping_others$sp,sep = "")

mapping_others_v2 <- mapping_others[match(matches_long$levelblocksp,mapping_others$levelblocksp),]
mapping_others_v2$match_prop <- matches_long$value

ggplot(mapping_others_v2)+
  geom_point(aes(x=mapping_rate,y=match_prop))+
  theme_classic()+
  facet_wrap(~sp,scales="free") # no relationship between the proportion of mapping to other genomes and the mapping rate of those genomes - at least by genome

