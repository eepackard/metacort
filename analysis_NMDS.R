library(readr)
library(tidyverse)
library(vegan)

express_table<-read_csv("clean_data/gene_interest_express_GC.csv")

norm_wide<-express_table[,-c(4,5,7:11)] |>  pivot_wider(names_from = "gene",values_from = KGD) #using KGD normalized data
norm_wide$GMC_sum<- rowSums(norm_wide[,grepl("GMC",colnames(norm_wide))])
norm_wide$OrgN_sum<- rowSums(norm_wide[,grepl("OrgN",colnames(norm_wide))])
norm_wide$AA<- rowSums(norm_wide[,which(colnames(norm_wide) %in% c("OrgN_AAAP","OrgN_ACT","OrgN_LAT","OrgN_YAT"))])
norm_wide$APET <- rowSums(norm_wide[,which(colnames(norm_wide) %in% c("OrgN_POT","OrgN_OPT"))])

c_norm_wide <- norm_wide[,-which(grepl("OrgN",colnames(norm_wide)))]
c_norm_wide <- c_norm_wide[,-which(grepl("GMC_1",colnames(c_norm_wide)))]


multi<- metaMDS(c_norm_wide[,-c(1,2,3,12)],autotransform = FALSE)
multi
data.scores <- as.data.frame(scores(multi,"sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- c_norm_wide$level  #  add the grp variable created earlier
head(data.scores)
gene.scores <- as.data.frame(scores(multi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
gene.scores$species <- rownames(gene.scores)  # create a column of species, from the rownames of species.scores
head(gene.scores) 

data.scores<-cbind(c_norm_wide$Block,data.scores)
colnames(data.scores)<-c("Block",colnames(data.scores[,-1]))


ggplot() + 
  geom_text(data=gene.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  geom_line(data=data.scores,aes(group = Block,x=NMDS1,y=NMDS2))+
  scale_colour_manual(values=c("high" = "red", "low" = "blue")) +
  coord_equal() +
  theme_bw()

#subset
select<-which(express_table[which(express_table$gene == "MnP" & express_table$level == "high"),]$KGD-express_table[which(express_table$gene == "MnP" & express_table$level == "low"),]$KGD > 3)
block<-as.data.frame(express_table[which(express_table$gene == "MnP" & express_table$level == "high"),]$Block)
block[c(select),]

express_table.2<-express_table[which(express_table$Block %in% block[c(select),]),]


norm_wide_sub<-express_table.2[,-c(4,5,7:11)] |>  pivot_wider(names_from = "gene",values_from = KGD) #using KGD normalized data
norm_wide_sub$GMC_sum<- rowSums(norm_wide_sub[,grepl("GMC",colnames(norm_wide_sub))])
norm_wide_sub$OrgN_sum<- rowSums(norm_wide_sub[,grepl("OrgN",colnames(norm_wide_sub))])
norm_wide_sub$AA<- rowSums(norm_wide_sub[,which(colnames(norm_wide_sub) %in% c("OrgN_AAAP","OrgN_ACT","OrgN_LAT","OrgN_YAT"))])
norm_wide_sub$APET <- rowSums(norm_wide_sub[,which(colnames(norm_wide_sub) %in% c("OrgN_POT","OrgN_OPT"))])

c_norm_wide_sub <- norm_wide_sub[,-which(grepl("OrgN",colnames(norm_wide_sub)))]
c_norm_wide_sub <- c_norm_wide_sub[,-which(grepl("GMC_1",colnames(c_norm_wide_sub)))]

multi_sub<- metaMDS(c_norm_wide_sub[,-c(1,2,3,12)],autotransform = FALSE)#or without MnP which I suppose shouldnt be included...
multi_sub
(multi_sub)
data.scores_sub <- as.data.frame(scores(multi_sub,"sites"))  #Using the scores_sub function from vegan to extract the site scores_sub and convert to a data.frame
data.scores_sub$site <- rownames(data.scores_sub)  # create a column of site names, from the rownames of data.scores_sub
data.scores_sub$grp <- c_norm_wide_sub$level  #  add the grp variable created earlier
head(data.scores_sub)
gene.scores_sub <- as.data.frame(scores(multi_sub, "species"))  #Using the scores_sub function from vegan to extract the species scores_sub and convert to a data.frame
gene.scores_sub$species <- rownames(gene.scores_sub)  # create a column of species, from the rownames of species.scores_sub
head(gene.scores_sub) 

data.scores_sub<-cbind(c_norm_wide_sub$Block,data.scores_sub)
colnames(data.scores_sub)<-c("Block",colnames(data.scores_sub[,-1]))

ggplot() + 
  geom_text(data=gene.scores_sub,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores_sub,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  #geom_text(data=data.scores_sub,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  geom_line(data=data.scores_sub,aes(group = Block,x=NMDS1,y=NMDS2))+
  scale_colour_manual(values=c("high" = "red", "low" = "blue")) +
  coord_equal() +
  theme_bw()
