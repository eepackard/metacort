library(readr)
library(tidyverse)
library(gridExtra)


#read in data
express_table<-read_csv("clean_data/gene_interest_express_GC.csv")


express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)


ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=sum_reads,colour = level))+
  theme_classic()

ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=UBC,colour = level))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=UBC,fill = level))+
  geom_line(aes(group = Block,x=level,y=UBC))+
  geom_point(aes(fill = level,group = Block,x=level,y=UBC))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=BTub_1,fill = level))+
  geom_line(aes(group = Block,x=level,y=BTub_1))+
  geom_point(aes(fill = level,group = Block,x=level,y=BTub_1))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=TEF_1,fill = level))+
  geom_line(aes(group = Block,x=level,y=TEF_1))+
  geom_point(aes(fill = level,group = Block,x=level,y=TEF_1))+
  theme_classic()


levelblock<-as.list(unique(express_table$levelblock))
gene<-as.list(unique(express_table$gene))
plots<-list()
for (i in 1:30){
plots[[i]]<- ggplot(express_table[which(express_table$levelblock == levelblock[[i]]),],aes(x=UBC,y=BTub_1))+geom_point()
}

for (i in 1:22){
  plots[[i]]<- ggplot(express_table[which(express_table$gene == gene[[i]]),],aes(x=UBC,y=sum_reads))+geom_point()
}

grid.arrange(grobs= plots)

