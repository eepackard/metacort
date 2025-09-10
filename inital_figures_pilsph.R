library(readr)
library(tidyverse)

#read in ----
express_table<-read_csv("clean_data/gene_interest_express_pilsph.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

# high vs low ----
ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=sum_TPM,colour = level))+
  theme_classic()


ggplot(express_table[which(express_table$gene == "KGD"),])+
  geom_boxplot(aes(x=level,y=sum_TPM,fill = level))+
  geom_line(aes(group = Block,x=level,y=sum_TPM))+
  geom_point(aes(fill = level,group = Block,x=level,y=sum_TPM))+
  theme_classic()


ggplot(express_table[which(express_table$gene == "NAG"),])+#& express_table$Block < 15
  geom_boxplot(aes(x=level,y=sum_TPM,fill = level))+
  geom_line(aes(group = Block,x=level,y=sum_TPM))+
  geom_point(aes(fill = level,group = Block,x=level,y=sum_TPM))+
  theme_classic()


#scatterplots ----

##make some wide dfs so that it is easier to plot the normalized read data

norm_wide<-express_table[,-5] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

#growth vs MnP
plot(norm_wide$KGD,norm_wide$NAG)+text(norm_wide$KGD,norm_wide$NAG,labels = paste(norm_wide$level,norm_wide$Block))
plot(norm_wide$KGD,norm_wide$CAT)+text(norm_wide$KGD,norm_wide$CAT,labels = paste(norm_wide$level,norm_wide$Block))

