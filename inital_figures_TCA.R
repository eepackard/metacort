library(readr)
library(tidyverse)

#read in ----
express_table_TCA<-read_csv("clean_data/gene_interest_express_TCA.csv")
express_table_GC<-read_csv("clean_data/gene_interest_express_GC.csv")

express_table_TCA$ecNum<-as.factor(express_table_TCA$ecNum)
express_table_TCA$level<-as.factor(express_table_TCA$level)
express_table_TCA$Block<-as.factor(express_table_TCA$Block)

express_table_GC$gene<-as.factor(express_table_GC$gene)
express_table_GC$level<-as.factor(express_table_GC$level)
express_table_GC$Block<-as.factor(express_table_GC$Block)


ggplot(express_table_TCA)+
  geom_boxplot(aes(x=ecNum,y=sum_reads,colour = level))+
  theme_classic()

ggplot(express_table_TCA)+
  geom_boxplot(aes(x=ecNum,y=sum_TPM,colour = level))+
  theme_classic()

express_table_TCA<-express_table_TCA[-which(express_table_TCA$Block == "block19"),] #this block is so variant

ggplot(express_table_TCA)+
  geom_boxplot(aes(x=ecNum,y=sum_TPM,colour = level))+
  theme_classic()

ggplot(express_table_TCA)+
  geom_boxplot(aes(x=ecNum,y=sum_reads,colour = level))+
  theme_classic()

#reduced to top 8 ----
##there is some Block/blocks where the difference is much stronger - lets limit to those where difference in norm_readsalized difference is greater than 2 
## Blocks - 1,4,7,8,11,12,15,20
express_table_TCA$Block<-gsub("block","",express_table_TCA$Block)
express_table.2<-express_table_TCA[which(express_table_TCA$Block %in% c(1,4,7,8,11,12,15,20)),]

express_table_GC$Block<-gsub("block","",express_table_GC$Block)
express_table.3<-express_table_GC[which(express_table_GC$Block %in% c(1,4,7,8,11,12,15,20)),]


ggplot(express_table.2)+
  geom_boxplot(aes(x=ecNum,y=sum_TPM,colour = level))+
  theme_classic()

ggplot(express_table.2)+
  geom_boxplot(aes(x=ecNum,y=norm_reads,colour = level))+
  theme_classic()

ggplot(express_table.3)+
  geom_boxplot(aes(x=gene,y=sum_TPM,colour = level))+
  theme_classic()

ggplot(express_table.3)+
  geom_boxplot(aes(x=gene,y=norm_reads,colour = level))+
  theme_classic()


#scatterplots ----

##make some wide dfs so that it is easier to plot the normalized read data

wide_TCA_sub<-express_table.2[,-c(5,6)] |>  pivot_wider(names_from = "ecNum",values_from = norm_reads)
wide_TCA_sub$allTCA<-rowSums(wide_TCA_sub[,3:14])

wide_TCA<-express_table_TCA[,-c(5,6)] |>  pivot_wider(names_from = "ecNum",values_from = norm_reads)
wide_TCA$allTCA<-rowSums(wide_TCA[,3:14])


ggplot(wide_TCA_sub)+
  geom_boxplot(aes(x=level,y=allTCA,fill = level))+
  geom_line(aes(group = Block,x=level,y=allTCA))+
  geom_point(aes(fill = level,group = Block,x=level,y=allTCA))+
  theme_classic()

ggplot(wide_TCA)+
  geom_boxplot(aes(x=level,y=allTCA,fill = level))+
  geom_line(aes(group = Block,x=level,y=allTCA))+
  geom_point(aes(fill = level,group = Block,x=level,y=allTCA))+
  theme_classic()

wide_sub<-express_table.3[,-c(5,6)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

#resp vs MnP
plot(wide_sub$MnP,wide_TCA_sub$allTCA)+text(wide_sub$MnP,wide_TCA_sub$allTCA,labels = paste(wide_sub$level,wide_sub$Block))+
  abline(lm(wide_TCA_sub$allTCA~wide_sub$MnP))
#growth vs resp
plot(wide_sub$GT48,wide_TCA_sub$allTCA)+text(wide_sub$GT48,wide_TCA_sub$allTCA,labels = paste(wide_sub$level,wide_sub$Block))
plot(wide_sub[-which(wide_sub$GT48 > 550),]$GT48,wide_TCA_sub[-which(wide_sub$GT48 > 550),]$allTCA)+text(wide_sub[-which(wide_sub$GT48 > 550),]$GT48,wide_TCA_sub[-which(wide_sub$GT48 > 550),]$allTCA,labels = paste(wide_sub[-which(wide_sub$GT48 > 550),]$level,wide_sub[-which(wide_sub$GT48 > 550),]$Block))

#CCat and SOD vs resp
plot(wide_sub$SOD+wide_sub$CAT,wide_TCA_sub$allTCA)+
  text(wide_sub$SOD+wide_sub$CAT,wide_TCA_sub$allTCA,labels = paste(wide_sub$level,wide_sub$Block))+
  abline(lm(wide_TCA_sub$allTCA~c(wide_sub$SOD+wide_sub$CAT)))


#CUE vs MnP
wide_TCA_sub$CUE<-wide_sub$GT48/wide_TCA_sub$allTCA
plot(wide_sub$MnP,wide_TCA_sub$CUE)+text(wide_sub$MnP,wide_TCA_sub$CUE,labels = paste(wide_sub$level,wide_sub$Block))+
  abline(lm(wide_TCA_sub$CUE~wide_sub$MnP))

