library(readr)
library(tidyverse)
library(patchwork)


##read in
HMM_results<- read_csv("clean_data/HMM_results.csv")
norm_wide <- read_csv("clean_data/norm_wide_clean.csv")

HMM_results <- HMM_results[match(norm_wide$levelblock,HMM_results$levelblock),]

HMM_results$levelblock == norm_wide$levelblock

HMM_results$MnP_exp <- norm_wide$MnP


plot(HMM_results$MnP_exp,HMM_results$tot_coromn)
plot(HMM_results$MnP_exp,HMM_results$Cor)
plot(HMM_results$MnP_exp,HMM_results$prop_corom)

HMM_wide <- pivot_longer(HMM_results[,c(1:24,29)],cols=1:24,names_to = "genus")

ggplot(HMM_wide,aes(level,value,fill = genus))+
  geom_bar(position = "stack",stat = "identity")

##read in
blast_results <- read_delim("clean_data/blast_1000_results.csv",delim = ";")
blast_results$block <- paste("block",blast_results$block, sep = "")
blast_results$levelblock <- paste(blast_results$level,blast_results$block,sep = "")

blast_results <- blast_results[match(HMM_results$levelblock,blast_results$levelblock),]

blast_results$prop_cort_tot <- blast_results$Cortinariaceae/blast_results$`total hits`
blast_results$prop_cort_fun <- blast_results$Cortinariaceae/blast_results$Fungi
blast_results$prop_fun_tot <- blast_results$Fungi/blast_results$`total hits`

plot(blast_results$prop_cort_tot,HMM_results$Cor)

plot(blast_results$prop_fun_tot,HMM_results$tot_AA2)

blast_results$levelblock == norm_wide$levelblock

plot(blast_results$prop_cort_fun,norm_wide$MnP)


write_csv(blast_results,"clean_data/blast_1000_results_clean.csv")
