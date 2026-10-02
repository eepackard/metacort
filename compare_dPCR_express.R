library(readr)
library(tidyverse)

#here I want to compare the measurements made with dPCR to the actual determined expression data 

#read in ----

express_table<-read_csv("clean_data/gene_interest_express_GC.csv")
dPCR<-read_delim("clean_data/dPCR_results_slected.csv",delim =";")

#compare reads to copies/µl
norm_wide<-express_table[,-c(5,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_TPM)
sum_wide<-express_table[,-c(4,7:12)] |>  pivot_wider(names_from = "gene",values_from = sum_reads)

#reorder so dPCR is in same order as express table
dPCR<-dPCR[match(paste(norm_wide$Block,norm_wide$level),paste(dPCR$Block,dPCR$level)),]
norm_wide$Block == dPCR$Block
sum_wide$Block == dPCR$Block

plot(sum_wide$MnP/sum_wide$KGD,dPCR$ratio)+text(sum_wide$MnP/sum_wide$KGD,dPCR$ratio,labels = paste(sum_wide$level,sum_wide$Block))
plot(norm_wide$MnP/norm_wide$KGD,dPCR$ratio)+text(norm_wide$MnP/norm_wide$KGD,dPCR$ratio,labels = paste(norm_wide$level,norm_wide$Block))

summary(lm(norm_wide$MnP/norm_wide$KGD~dPCR$ratio))
#pretty good linear relationship here

plot_df<-as.data.frame(cbind(norm_wide$MnP/norm_wide$KGD,dPCR$ratio))
colnames(plot_df)<-c("RNA-seq MnP/KGD gene expression ratio", "dPCR MnP/KGD copy numbers ratio")

tiff("figures/supp_ratios.tiff")
ggplot(aes(x= `RNA-seq MnP/KGD gene expression ratio`, y= `dPCR MnP/KGD copy numbers ratio`),data = plot_df)+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

