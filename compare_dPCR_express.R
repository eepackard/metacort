library(readr)
library(tidyverse)

#here I want to compare the measurements made with dPCR to the actual determined expression data 

#read in ----

express_table<-read_csv("clean_data/gene_interest_express_GC.csv")
dPCR<-read_delim("clean_data/dPCR_results_slected.csv",delim =";")

#compare reads to copies/µl
norm_wide<-express_table[,-c(4,5)] |>  pivot_wider(names_from = "gene",values_from = norm_reads)
sum_wide<-express_table[,-c(4,8)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)

#reorder so dPCR is in same order as express table
dPCR<-dPCR[match(paste(norm_wide$Block,norm_wide$level),paste(dPCR$Block,dPCR$level)),]
norm_wide$Block == dPCR$Block
sum_wide$Block == dPCR$Block

plot(norm_wide$MnP,dPCR$MnP_conc_cps_µl)+text(norm_wide$MnP,dPCR$MnP_conc_cps_µl,labels = paste(norm_wide$level,norm_wide$Block))
plot(sum_wide$MnP,dPCR$MnP_conc_cps_µl)+text(sum_wide$MnP,dPCR$MnP_conc_cps_µl,labels = paste(sum_wide$level,sum_wide$Block))

plot(sum_wide$KGD,dPCR$KGD_conc_cps_µl)+text(sum_wide$KGD,dPCR$KGD_conc_cps_µl,labels = paste(sum_wide$level,sum_wide$Block))

plot(sum_wide$MnP/sum_wide$KGD,dPCR$ratio)+text(sum_wide$MnP/sum_wide$KGD,dPCR$ratio,labels = paste(sum_wide$level,sum_wide$Block))
#pretty good linear relationship here