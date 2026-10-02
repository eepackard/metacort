library(readr)
library(tidyverse)
library(ggrepel)

#read in the expression data and the data from JGI

express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")


Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)


A1_prot_Coromn <- read_delim("raw_data/A1_prot_Coromn.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                                `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                                Location = col_skip(), Scaffold = col_skip(), 
                                                                                                                Start = col_skip(), End = col_skip(), 
                                                                                                                Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

OPT_Coromn <- read_delim("raw_data/OPT_2.A.67.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                        `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                        Location = col_skip(), Scaffold = col_skip(), 
                                                                                                        Start = col_skip(), End = col_skip(), 
                                                                                                        Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

AAAP_Coromn <- read_delim("raw_data/AAAP_2.A.18.csv",delim = ";", escape_double = FALSE, col_types = cols(Score = col_skip(), 
                                                                                                          `Organism Name` = col_skip(), Track = col_skip(), 
                                                                                                          Location = col_skip(), Scaffold = col_skip(), 
                                                                                                          Start = col_skip(), End = col_skip(), 
                                                                                                          Strand = col_skip(), `User Annotations` = col_skip()),trim_ws = TRUE)

transporters <- read_delim("raw_data/transporters.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)

Coromn1_FilteredModels1_go_2025_04_22 <- read_delim("raw_data/Coromn1_FilteredModels1_go_2025-04-22.tab", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)


Coromn1_FilteredModels1_domain_2025_04_22 <- read_delim("raw_data/Coromn1_FilteredModels1_domain_2025-04-22.tab", 
                                                        delim = "\t", escape_double = FALSE, 
                                                        trim_ws = TRUE)

OrgN<-read_csv("clean_data/N_trans_clean.csv")
GMC<-read_csv("clean_data/GMC_clean.csv")

#need to remove block19
express_all<-express_all[-which(express_all$Block == "block19"),]

#want to pull out genes that either have EC, I have pulled out otherwise (GMC, OrgN), or are proteases (maybe have ec)

express_all_EC<- express_all[which(express_all$protID %in% Coromn1_FilteredModels1_ec$proteinId),]
Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% A1_prot_Coromn$`Protein Id`),] #EC 3.4.23. are aspartic proteases
A1_prot_Coromn[-which(A1_prot_Coromn$`Protein Id` %in% Coromn1_FilteredModels1_ec$proteinId),] ## there is five prots that do not have EC but are A1
express_all_EC$ecNum <- Coromn1_FilteredModels1_ec[match(express_all_EC$protID,Coromn1_FilteredModels1_ec$proteinId),]$ecNum

express_all_exprot<-express_all[which(express_all$protID %in% A1_prot_Coromn[-which(A1_prot_Coromn$`Protein Id` %in% Coromn1_FilteredModels1_ec$proteinId),]$`Protein Id`),]
express_all_exprot$ecNum <- rep("3.4.23.-",nrow(express_all_exprot))# not sure exactly which type they fall into

Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% GMC$proteinId),] #all 1.1.3.13 are already there
express_all_GMC<- express_all[which(express_all$protID %in% GMC$proteinId),]
express_all_GMC$ecNum <- GMC[match(express_all_GMC$protID,GMC$proteinId),]$based_on_tree# not sure exactly which type they fall into

express_all_lacc<- express_all[which(express_all$protID %in% c("2079975","1581758","1821865","1205939","1668803")),]#picking out laccases based on ProtID I have confirmed
express_all_lacc$ecNum <- c(rep("Lacc",nrow(express_all_lacc)))

express_all_mnp<- express_all[which(express_all$protID %in% c("1753449","2106982","2121138","1733174","1773443","2004932")),]#picking out mnp based on ProtID because now two most expressed that i change do not have ECNUM
express_all_mnp$ecNum <- c(rep("mnp",nrow(express_all_mnp)))


Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% OrgN$proteinId),] #non have EC
Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% AAAP_Coromn$`Protein Id`),] #non have EC
Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% OPT_Coromn$`Protein Id`),]#non have EC
Coromn1_FilteredModels1_ec[which(Coromn1_FilteredModels1_ec$proteinId %in% transporters$`Protein Id`),]#non have EC
express_all_orgN<- express_all[which(express_all$protID %in% OrgN$proteinId),]
express_all_AAAP<- express_all[which(express_all$protID %in% AAAP_Coromn$`Protein Id`),]
express_all_OPT<- express_all[which(express_all$protID %in% OPT_Coromn$`Protein Id`),]
##first there will be overlap in my full transporters DB so I should remove dups
transporters<-transporters[-which(transporters$`Protein Id` %in% c(express_all_orgN$protID,express_all_AAAP$protID,express_all_OPT$protID)),]
express_all_tother<- express_all[which(express_all$protID %in% transporters$`Protein Id`),]
express_all_tother$ecNum<- transporters[match(express_all_tother$protID,transporters$`Protein Id`),]$type
express_all_orgN$ecNum <- OrgN[match(express_all_orgN$protID,OrgN$proteinId),]$gene# not sure exactly which type they fall into
OPT_Coromn$gene <- c(rep("OPT",nrow(OPT_Coromn)))
express_all_OPT$ecNum <- OPT_Coromn[match(express_all_OPT$protID,OPT_Coromn$`Protein Id`),]$gene# not sure exactly which type they fall into
AAAP_Coromn$gene <- c(rep("AAAP",nrow(AAAP_Coromn)))
express_all_AAAP$ecNum <- AAAP_Coromn[match(express_all_AAAP$protID,AAAP_Coromn$`Protein Id`),]$gene

#bind all
express_select<-rbind(express_all_EC[-which(express_all_EC$ecNum == "1.1.3.13"),],express_all_mnp,express_all_exprot,express_all_lacc,express_all_orgN,express_all_OPT,express_all_AAAP,express_all_tother,express_all_GMC)

#will try to methods one that is prot by prot and one with them combined by EC

#make summed df 
express_select_sum<- express_select |>  group_by(ecNum,level,Block) |> summarise(sum_TPM = sum(TPM),sum_reads = sum(NumReads))

#add levelblock factor
express_all$levelblock <- paste(express_all$level,express_all$Block,sep = "")
express_select_sum$levelblock <- paste(express_select_sum$level,express_select_sum$Block,sep = "")

#make wide dfs
express_wide_sum<- express_select_sum[,-c(2,3,5)] |>  pivot_wider(names_from = "ecNum",values_from = sum_TPM)
express_wide<- express_all[,-c(1,4,5)] |>  pivot_wider(names_from = "protID",values_from = TPM)


#pull out MnP
MnP_sum<-express_wide_sum[,which(colnames(express_wide_sum) == "mnp")]
express_wide_sum<-express_wide_sum[,-which(colnames(express_wide_sum) == "mnp")]

MnP<-as.data.frame(rowSums(express_wide[,which(colnames(express_wide) %in% express_all_mnp$protID)]))
colnames(MnP)<-"MnP"
express_wide<-express_wide[,-which(colnames(express_wide) %in% express_all_mnp$protID)]

#make first column row names - note that the different dataframes are in different orders!
samples_sum<-express_wide_sum$levelblock
express_wide_sum<-express_wide_sum[,-1]
rownames(express_wide_sum)<-samples_sum

samples<- express_wide$levelblock
express_wide<-express_wide[,-1]
rownames(express_wide)<-samples

#remove any columns where the abudnace is less than 100
express_wide_sum <- express_wide_sum[,-which(colSums(express_wide_sum)< 100)]
express_wide <- express_wide[,-which(colSums(express_wide)< 100)]

#should transform all to be safe

MnP$MnP<-sqrt(MnP$MnP)
MnP_sum$MnP<-sqrt(MnP_sum$mnp)

express_wide <- express_wide |> mutate_all(sqrt)
express_wide_sum <- express_wide_sum |> mutate_all(sqrt)
#express_wide_sum[sapply(express_wide_sum, is.infinite)]<-NA

#now can test
cordfsum<-as.data.frame(cor(express_wide_sum,MnP_sum,method = "pearson"))
cordfsum$ec<-rownames(cordfsum)
cordfsum<-cordfsum[order(cordfsum$MnP,decreasing = TRUE),]
cordfsum$rank<- c(1:length(cordfsum$MnP))

cordf<-as.data.frame(cor(express_wide,MnP,method = "pearson"))
cordf$prot<-rownames(cordf)
cordf<-cordf[order(cordf$MnP,decreasing = TRUE),]
cordf$rank<- c(1:length(cordf$MnP))

#lets match some info to make it easier to inspect
#but some have multiple so i will combine
iprID<-aggregate(iprId ~ proteinId, Coromn1_FilteredModels1_domain_2025_04_22,FUN = "paste")
iprDesc<-aggregate(iprDesc ~ proteinId, Coromn1_FilteredModels1_domain_2025_04_22,FUN = "paste")
ipr<-cbind(iprID,iprDesc)
ipr<-ipr[,-3]

goname<-aggregate(go_name ~ proteinId, Coromn1_FilteredModels1_go_2025_04_22,FUN = "paste")
goacc<-aggregate(go_acc ~ proteinId, Coromn1_FilteredModels1_go_2025_04_22,FUN = "paste")
go<-cbind(goname,goacc)
go<-go[,-3]


cordf$goname<-go[match(cordf$prot,go$proteinId),]$go_name
cordf$goacc<-go[match(cordf$prot,go$proteinId),]$go_acc
cordf$iprID<-ipr[match(cordf$prot,ipr$proteinId),]$iprId
cordf$iprDesc<-ipr[match(cordf$prot,ipr$proteinId),]$iprDesc

cordfsum$def<- Coromn1_FilteredModels1_ec[match(cordfsum$ec,Coromn1_FilteredModels1_ec$ecNum),]$definition
cordfsum$path<-Coromn1_FilteredModels1_ec[match(cordfsum$ec,Coromn1_FilteredModels1_ec$ecNum),]$pathway
cordfsum$catact<- Coromn1_FilteredModels1_ec[match(cordfsum$ec,Coromn1_FilteredModels1_ec$ecNum),]$catalyticActivity

##lets try to sort out the ones i looked at a priori
express_interest<-read_csv("clean_data/gene_interest_express_GC_all.csv")
express_interest$ecNum<-Coromn1_FilteredModels1_ec[match(express_interest$protID,Coromn1_FilteredModels1_ec$proteinId),]$ecNum

cordf$apriori<-if_else(cordf$prot %in% express_interest$protID, "red","black")
cordfsum$apriori<-if_else(cordfsum$ec %in% express_interest$ecNum, "red",if_else(cordfsum$ec %in% c("OPT","ACT","AAAP","YAT","LAT","Nitrate","ammonium","DHA","Urea","AAT","Lacc","AAO","GDH","AOx"),"red",if_else(cordfsum$ec %in% c("1.1.3.7","1.1.99.18"),"red","black")))


ggplot(cordf)+
  geom_point(aes(x=rank,y=MnP,colour = apriori),alpha=0.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()

#pdf(file ="figures/figures_s8.pdf")
ggplot(cordfsum)+
  geom_point(aes(x=rank,y=MnP,colour = apriori),alpha=0.5)+
  scale_color_manual(values = c("black","red"))+
  annotate(geom = "segment",x=0,xend = 800,y=-0.38540,yend = -0.38540)+
  annotate(geom = "segment",x=0,xend = 800,y=0.299682,yend = 0.299682,)+
  geom_text_repel(data = cordfsum[c(3,4,6,9,10,11,21,721,745,750,788),],aes(x=rank,y=MnP,label = ec),nudge_x = 100)+
  theme_classic()


#now lets looks at some of the top hits that I didnt look at before
ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.59`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.59`),method = "lm")+
  ylab(label = expression("Glucanase 1,3-"*alpha*"glucan (GH71)"))+
  xlab(label = "MnP ")+
  theme_classic() #GH71

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.28`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.28`),method = "lm")+
  ylab(label = "Trehalase")+
  xlab(label = "MnP ")+
  theme_classic() #trehalose breakdown? 

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.1.3.12`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.1.3.12`),method = "lm")+
  ylab(label = "TTP")+
  xlab(label = "MnP ")+
  theme_classic() #trehalose synthese? 

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`2.4.1.15`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`2.4.1.15`),method = "lm")+
  ylab(label = "TPS")+
  xlab(label = "MnP ")+
  theme_classic() #trehalose syn? 

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.28`/express_wide_sum$`2.4.1.15`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.28`/express_wide_sum$`2.4.1.15`),method = "lm")+
  ylab(label = "ratio Trehalaase/TPS")+
  xlab(label = "MnP ")+
  theme_classic() #trehalose syn? 

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.28`/express_wide_sum$`2.4.1.15`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.2.1.28`/express_wide_sum$`2.4.1.15`),method = "lm")+
  ylab(label = "ratio Trehalaase/TPS")+
  xlab(label = "MnP ")+
  theme_classic() #trehalose syn? 



ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`2.3.2.2`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`2.3.2.2`),method = "lm")+
  ylab(label = "")+
  xlab(label = "MnP ")+
  theme_classic() #GGT

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.5.1.49`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.5.1.49`),method = "lm")+
  ylab(label = "")+
  xlab(label = "MnP ")+
  theme_classic() #Forminidase

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.4.14.9`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.4.14.9`),method = "lm")+
  ylab(label = "Peptidase")+
  xlab(label = "MnP ")+
  theme_classic() #tripeptidtlydase?

#tiff(filename="figures/figures_s6.tiff")
ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.1.3.2`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`3.1.3.2`),method = "lm")+
  ylab(label = "Acid phosphatase")+
  xlab(label = "MnP")+
  theme_classic() #acid phosphatase

ggplot()+
  geom_point(aes(y=express_wide_sum$`3.2.1.52`,x=MnP_sum$MnP))+
  geom_smooth(aes(y=express_wide_sum$`3.2.1.52`,x=MnP_sum$MnP),method = "lm")+
  ylab(label = "NAG")+
  xlab(label = "MnP")+
  theme_classic() #NAG

#tiff(filename="figures/figures_s5.tiff")
ggplot()+
  geom_point(aes(y=express_wide_sum$`3.2.1.14`,x=MnP_sum$MnP))+
  geom_smooth(aes(y=express_wide_sum$`3.2.1.14`,x=MnP_sum$MnP),method = "lm")+
  ylab(label = "Chitinase")+
  xlab(label = "MnP")+
  theme_classic() #chitinase

ggplot()+
  geom_point(aes(x=express_wide_sum$`3.2.1.52`,y=express_wide_sum$`3.1.3.2`))+
  geom_smooth(method = "lm",se=FALSE,linetype=2)+
  ylab(label = "aP")+
  xlab(label = "NAG ")+
  theme_classic() #acid phosphatase v NAG

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`6.3.1.2`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`6.3.1.2`),method = "lm")+
  ylab(label = "GS")+
  xlab(label = "MnP")+
  theme_classic() #gluamate-ammonia ligase

ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$`1.4.1.4`))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$`1.4.1.4`),method = "lm")+
  ylab(label = "GDH")+
  xlab(label = "MnP")+
  theme_classic() #gluamate-ammonia ligase

ggplot()+
  geom_point(aes(x=express_wide_sum$`1.2.4.2`,y=express_wide_sum$`1.4.1.4`))+
  geom_smooth(aes(x=express_wide_sum$`1.2.4.2`,y=express_wide_sum$`1.4.1.4`),method = "lm")+
  ylab(label = "GDH")+
  xlab(label = "KGD")+
  theme_classic() #gluamate-ammonia ligase

ggplot()+
  geom_point(aes(x=express_wide_sum$`6.3.1.2`,y=express_wide_sum$`1.4.1.4`))+
  geom_smooth(aes(x=express_wide_sum$`6.3.1.2`,y=express_wide_sum$`1.4.1.4`),method = "lm")+
  ylab(label = "GDH")+
  xlab(label = "GS")+
  theme_classic()

ggplot()+
  geom_point(aes(x=express_wide_sum$ammonium,y=express_wide_sum$`6.3.1.2`))+
  geom_smooth(aes(x=express_wide_sum$ammonium,y=express_wide_sum$`6.3.1.2`),method = "lm")+
  ylab(label = "AmT")+
  xlab(label = "GS")+
  theme_classic() #gluamate-ammonia ligase v ammonium transport

ggplot()+
  geom_point(aes(x=express_wide_sum$ammonium,y=express_wide_sum$`1.4.1.4`))+
  geom_smooth(aes(x=express_wide_sum$ammonium,y=express_wide_sum$`1.4.1.4`),method = "lm")+
  ylab(label = "AmT")+
  xlab(label = "GDH")+
  theme_classic()

ggplot()+
  geom_point(aes(x=express_wide_sum$ammonium,y=express_wide_sum$OPT))+
  geom_smooth(aes(x=express_wide_sum$ammonium,y=express_wide_sum$OPT),method = "lm")+
  ylab(label = "AmT")+
  xlab(label = "OPT")+
  theme_classic() #gluamate-ammonia ligase v ammonium transport
cor.test(express_wide_sum$ammonium,express_wide_sum$OPT)

# just to take a closer look at the EC classes
cordfsum <-separate(cordfsum,col = ec, into = c("class","subclass","subsubclass","spe"))

#check out acid phosphotases
#based on jgi 2106987, 1584111, 1753485, 1816830, 1822257
prot_acid<- c("2106987", "1584111", "1753485", "1816830", "1822257")
express_acidpho<- express_all[which(express_all$protID %in% prot_acid),]
#only 1822257 and 1753485 have sigP above 80 %

ggplot(express_acidpho,aes(x=level,y=TPM,fill = as.factor(protID)))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly one expressed

#tiff(filename="figures/figures_s6.tiff")
ggplot()+
  geom_point(aes(x=MnP$MnP,y=express_wide$`1753485`))+
  geom_smooth(aes(x=MnP$MnP,y=express_wide$`1753485`),method = "lm")+
  ylab(label = "Acid phosphatase")+
  xlab(label = "MnP ")+
  theme_classic() #acid phosphatase

#check out alpha-glucosidase
#based on jgi 
alphag<- c("2014010", "1748086", "1770501", "1958061", "1979606","1591610", "1749380", "1981828")
express_alphag<- express_all[which(express_all$protID %in% alphag),]

ggplot(express_alphag,aes(x=level,y=TPM,fill = as.factor(protID)))+
  geom_boxplot(position = position_dodge(0.8))+
  theme_classic()# mostly two expressed

ggplot()+
  geom_point(aes(x=MnP$MnP,y=express_wide$`1749380`))+
  geom_smooth(aes(x=MnP$MnP,y=express_wide$`1749380`),method = "lm")+
  ylab(label = "")+
  xlab(label = "MnP ")+
  theme_classic() #alpha-glucosidase

#check ammonium transport
#tiff(filename="figures/figures_s7.tiff")
ggplot()+
  geom_point(aes(x=MnP_sum$MnP,y=express_wide_sum$ammonium))+
  geom_smooth(aes(x=MnP_sum$MnP,y=express_wide_sum$ammonium),method = "lm")+
  ylab(label = "Ammonium transporters")+
  xlab(label = "MnP ")+
  theme_classic() #all in 1.A.11.3

ggplot()+
  geom_point(aes(x=MnP$MnP,y=express_wide$`1685047`))+
  geom_smooth(aes(x=MnP$MnP,y=express_wide$`1685047`),method = "lm")+
  ylab(label = "Ammonium transporter (1.A.11.3.3)")+
  xlab(label = "MnP ")+
  theme_classic() #1.A.11.3.3 - high affinity

ggplot()+
  geom_point(aes(x=MnP$MnP,y=express_wide$`1870932`))+
  geom_smooth(aes(x=MnP$MnP,y=express_wide$`1870932`),method = "lm")+
  ylab(label = "")+
  xlab(label = "MnP ")+
  theme_classic() #1.A.11.3.3 - high affinity

ggplot()+
  geom_point(aes(x=MnP$MnP,y=express_wide$`1204656`))+
  geom_smooth(aes(x=MnP$MnP,y=express_wide$`1204656`),method = "lm")+
  ylab(label = "")+
  xlab(label = "MnP ")+
  theme_classic() #1.A.11.2.2

ggplot()+
  geom_point(aes(x=MnP$MnP,y=express_wide$`1795000`))+
  geom_smooth(aes(x=MnP$MnP,y=express_wide$`1795000`),method = "lm")+
  ylab(label = "")+
  xlab(label = "MnP ")+
  theme_classic() #1.A.11.3.4


write_csv(cordfsum,"clean_data/explore_cor.csv")
