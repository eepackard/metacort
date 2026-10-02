library(readr)
library(tidyverse)
library(vegan)


express_table<-read_csv("clean_data/expression_salmon_quant_all_sweet.csv")

express_table$Name<-as.factor(express_table$Name)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)
express_table$level_Block<-paste(express_table$level,express_table$Block)

#because I removed in other analysis I need to be consistent
express_table<-express_table[-which(express_table$Block == "block19"),]

norm_wide<-express_table[,-c(2,3,5,8)] |>  pivot_wider(names_from = "Name",values_from = TPM)


#some are not expressed at all

unexpressed<-norm_wide[,which(colSums(norm_wide[,-c(1:2)]) == 0)+2]
sum(colSums(unexpressed)) #double check
expressed<-norm_wide[,-which(colnames(norm_wide) %in% colnames(unexpressed))]
ncol(unexpressed)+ncol(expressed)==ncol(norm_wide)#double check nothing missed

#is the clade III sweet expressed? 
"Picab_PA_chr02_G000629.mRNA.1" %in% colnames(expressed)
"Picab_PA_chr02_G000629.mRNA.1" %in% colnames(unexpressed) #not expressed... 

# maybe I make a little ordination - just to help see if there is some patterens between high and low 


multi<- metaMDS(expressed[,-c(1,2)],autotransform = TRUE)
multi
data.scores <- as.data.frame(scores(multi,"sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- expressed$level  #  add the grp variable created earlier
head(data.scores)
gene.scores <- as.data.frame(scores(multi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
gene.scores$species <- rownames(gene.scores)  # create a column of species, from the rownames of species.scores
head(gene.scores) 

data.scores<-cbind(expressed$Block,data.scores)
colnames(data.scores)<-c("Block",colnames(data.scores[,-1]))


ggplot() + 
  geom_text(data=gene.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  geom_line(data=data.scores,aes(group = Block,x=NMDS1,y=NMDS2))+
  scale_colour_manual(values=c("high" = "red", "low" = "blue")) +
  coord_equal() +
  theme_bw()

# lets now split to pine and spruce

expressed_PS<-expressed[,c(1,2,grep("_PS_",colnames(expressed)))]
expressed_PA<-expressed[,c(1,2,grep("_PA_",colnames(expressed)))]
expressed_po<-expressed[,c(1,2,which(colnames(expressed) == "Potra_Potra2n1c3346.3"))]

#now back to long format
expressed_PS_long<- pivot_longer(expressed_PS,cols = c(3:ncol(expressed_PS)))
expressed_PA_long<- pivot_longer(expressed_PA,cols = c(3:ncol(expressed_PA)))
expressed_po_long<- pivot_longer(expressed_po,cols = c(3:ncol(expressed_po)))

#maybe I can also make a filtered set where only genes that are expressed in at least 3 samples 
#but then I need a presence/absence to filter 
expressed_PS_01<-expressed_PS[,-c(1,2)]
expressed_PS_01[expressed_PS_01>0]<-1
expressed_PS_f<-expressed_PS[,which(colSums(expressed_PS_01)>5)+2]
expressed_PS_f<-cbind(expressed_PS[,1:2],expressed_PS_f)

expressed_PA_01<-expressed_PA[,-c(1,2)]
expressed_PA_01[expressed_PA_01>0]<-1
expressed_PA_f<-expressed_PA[,which(colSums(expressed_PA_01)>5)+2]
expressed_PA_f<-cbind(expressed_PA[,1:2],expressed_PA_f)

#check if there is mostly overall in samples with PS or PA expression
plot(rowSums(expressed_PA_01[,-c(1,2)]),rowSums(expressed_PS_01[,-c(1,2)]))

#now back to long format
expressed_PS_long_f<- pivot_longer(expressed_PS_f,cols = c(3:ncol(expressed_PS_f)))
expressed_PA_long_f<- pivot_longer(expressed_PA_f,cols = c(3:ncol(expressed_PA_f)))


#plot

ggplot(expressed_PS_long_f)+
  geom_boxplot(aes(x=name,y=value,colour = level))+
  #ylim(0,200000)+ #just so I can see the other... 
  theme_classic()

ggplot(expressed_PA_long_f)+
  geom_boxplot(aes(x=name,y=value,colour = level))+
  theme_classic()

ggplot(expressed_po_long)+
  geom_boxplot(aes(x=name,y=value,colour = level))+
  theme_classic()


#one is highly expressed
expressed_PS[,36]


ggplot(norm_wide)+
  geom_boxplot(aes(x=level,y=Pinsy_PS_chr08_G035767.mRNA.1,fill = level))+
  geom_line(aes(group = Block,x=level,y=Pinsy_PS_chr08_G035767.mRNA.1))+
  geom_point(aes(fill = level,group = Block,x=level,y=Pinsy_PS_chr08_G035767.mRNA.1))+
  theme_classic()

#now read in cororm data so i can compare the MnP expression to root
norm_wide_coromn<-read_csv("clean_data/norm_wide_clean.csv")

#before binding make sure are in same order 
paste(expressed_PS$level,expressed_PS$Block,sep = "") == norm_wide_coromn$levelblock
expressed_PS$levelblock <-paste(expressed_PS$level,expressed_PS$Block,sep = "")
norm_wide_coromn<-norm_wide_coromn[match(expressed_PS$levelblock,norm_wide_coromn$levelblock),]
paste(expressed_PS$level,expressed_PS$Block,sep = "") == norm_wide_coromn$levelblock #double check

expressed_PS$MnP<-norm_wide_coromn$MnP

ggplot(expressed_PS,aes(x=sqrt(MnP),y=sqrt(Pinsy_PS_chr08_G035767.mRNA.1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

cor.test(sqrt(expressed_PS$MnP),sqrt(expressed_PS$Pinsy_PS_chr08_G035767.mRNA.1))

express_PS_sqrt <- expressed_PS[,-c(1,2,54,55)] |> mutate_all(sqrt)

corsum<-as.data.frame(cor(express_PS_sqrt,sqrt(expressed_PS$MnP),method = "pearson"))

ggplot(expressed_PS,aes(x=sqrt(MnP),y=sqrt(Pinsy_PS_chr07_G028798.mRNA.2)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()


#try ordination again wih split
multi<- metaMDS(expressed[,-c(1,2)],autotransform = FALSE)
multi
data.scores <- as.data.frame(scores(multi,"sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- expressed$level  #  add the grp variable created earlier
head(data.scores)
gene.scores <- as.data.frame(scores(multi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
gene.scores$species <- rownames(gene.scores)  # create a column of species, from the rownames of species.scores
head(gene.scores) 

data.scores<-cbind(expressed$Block,data.scores)
colnames(data.scores)<-c("Block",colnames(data.scores[,-1]))


ggplot() + 
  geom_text(data=gene.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  geom_line(data=data.scores,aes(group = Block,x=NMDS1,y=NMDS2))+
  scale_colour_manual(values=c("high" = "red", "low" = "blue")) +
  coord_equal() +
  theme_bw()

multi<- metaMDS(expressed_PA[,-c(1,2)],autotransform = FALSE)
multi
data.scores <- as.data.frame(scores(multi,"sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- expressed_PA$level  #  add the grp variable created earlier
head(data.scores)
gene.scores <- as.data.frame(scores(multi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
gene.scores$species <- rownames(gene.scores)  # create a column of species, from the rownames of species.scores
head(gene.scores) 

data.scores<-cbind(expressed_PA$Block,data.scores)
colnames(data.scores)<-c("Block",colnames(data.scores[,-1]))


ggplot() + 
  geom_text(data=gene.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  geom_line(data=data.scores,aes(group = Block,x=NMDS1,y=NMDS2))+
  scale_colour_manual(values=c("high" = "red", "low" = "blue")) +
  coord_equal() +
  theme_bw()
