library(readr)
library(tidyverse)

#read in ----
express_table<-read_csv("clean_data/gene_interest_express.csv")

express_table$gene<-as.factor(express_table$gene)
express_table$level<-as.factor(express_table$level)
express_table$Block<-as.factor(express_table$Block)

# high vs low ----
ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=norm_reads,colour = level))+
  theme_classic()

ggplot(express_table)+
  geom_boxplot(aes(x=gene,y=sum_TPM,colour = level))+
  theme_classic()


ggplot(express_table[which(express_table$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "KGD"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()


ggplot(express_table[which(express_table$gene == "GT48"),])+#& express_table$Block < 15
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "SOD"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "GMC"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "GLY"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "SOD"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "ALAS"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table[which(express_table$gene == "CHIT"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

#reduced to top 8 ----
##there is some Block/blocks where the difference is much stronger - lets limit to those where difference in norm_readsalized difference is greater than 2 
## Blocks - 1,2,4,8,9,12,13,14

express_table$Block<-gsub("block","",express_table$Block)
express_table.2<-express_table[which(express_table$Block %in% c(1,2,4,8,9,12,13,14)),]

ggplot(express_table.2)+
  geom_boxplot(aes(x=gene,y=norm_reads,colour = level))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "GMC"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "MnP"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "GLY"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "GT48"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "ALAS"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "SOD"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "CAT"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

ggplot(express_table.2[which(express_table.2$gene == "CHIT"),])+
  geom_boxplot(aes(x=level,y=norm_reads,fill = level))+
  geom_line(aes(group = Block,x=level,y=norm_reads))+
  geom_point(aes(fill = level,group = Block,x=level,y=norm_reads))+
  theme_classic()

#scatterplots ----

##make some wide dfs so that it is easier to plot the normalized read data

norm_wide<-express_table[,-c(4,5)] |>  pivot_wider(names_from = "gene",values_from = norm_reads)

#growth vs MnP
plot(norm_wide$MnP,norm_wide$GT48)+text(norm_wide$MnP,norm_wide$GT48,labels = paste(norm_wide$level,norm_wide$Block))
#resp vs MnP
plot(norm_wide$MnP,norm_wide$KGD)+text(norm_wide$MnP,norm_wide$KGD,labels = paste(norm_wide$level,norm_wide$Block))
#chitin break down vs MnP
plot(norm_wide$MnP,norm_wide$NAG)+text(norm_wide$MnP,norm_wide$NAG,labels = paste(norm_wide$level,norm_wide$Block))
plot(norm_wide$MnP,norm_wide$CHIT)+text(norm_wide$MnP,norm_wide$CHIT,labels = paste(norm_wide$level,norm_wide$Block))
norm_wide$NAGCHIT<-norm_wide$NAG+norm_wide$CHIT
plot(norm_wide$MnP,norm_wide$NAGCHIT)+text(norm_wide$MnP,norm_wide$NAGCHIT,labels = paste(norm_wide$level,norm_wide$Block))
#relationship between NAG and endo chitinase
plot(norm_wide$NAG,norm_wide$CHIT)+text(norm_wide$NAG,norm_wide$CHIT,labels = paste(norm_wide$level,norm_wide$Block))
#growth vs respiration
plot(norm_wide$GT48,norm_wide$KGD)+text(norm_wide$GT48,norm_wide$KGD,labels = paste(norm_wide$level,norm_wide$Block))
#Catalase vs super oxide dismutase
plot(norm_wide$SOD,norm_wide$CAT)+text(norm_wide$SOD,norm_wide$CAT,labels = paste(norm_wide$level,norm_wide$Block))

#CUE vs MnP
norm_wide$CUE<-norm_wide$GT48/norm_wide$KGD
plot(norm_wide$MnP,norm_wide$CUE)+text(norm_wide$MnP,norm_wide$CUE,labels = paste(norm_wide$level,norm_wide$Block))

ggplot(norm_wide)+
  geom_point(aes(x=CUE,y=MnP))+
  theme_classic()
##test if its a unimodel relationship
fit<-lm(MnP~CUE+I(CUE^2),data=norm_wide)
prd <- data.frame(CUE = seq(from = range(norm_wide$CUE)[1], to = range(norm_wide$CUE)[2], length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = CUE, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = norm_wide, aes(x = CUE, y = MnP))


norm_wide$NAGCHIT<-norm_wide$NAG+norm_wide$CHIT
plot(norm_wide$MnP,norm_wide$NAGCHIT)+text(norm_wide$MnP,norm_wide$NAGCHIT,labels = paste(norm_wide$level,norm_wide$Block))

