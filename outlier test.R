library(readr)
library(tidyverse)


express_table<-read_csv("clean_data/gene_interest_express_GC.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

#express_table<-express_table[-which(express_table$Block == "block19"),]


##wide dfs ----
norm_wide<-express_table[,-c(5,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)

##CUE ----
#make ratio
norm_wide$CUE<-norm_wide$GT48/norm_wide$KGD

stats_mean <- norm_wide[,-c(1:3)] |> summarise(across(everything(),mean))
stats_sd <- norm_wide[,-c(1:3)] |> summarise(across(everything(),sd))
upper<- stats_mean+stats_sd*1.5
lower<- stats_mean-stats_sd*1.5

all<-rbind(upper,lower,norm_wide[,-c(1:3)])
empty<- data.frame(level=c("upper","lower"),Block=c("upper","lower"),levelblock=c("upper","lower")) 
all<-cbind(rbind(empty,norm_wide[,1:3]),all)
