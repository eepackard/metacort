library(readr)
library(tidyverse)
library(DESeq2)
library(apeglm)
library(pheatmap)

#read in the expression data and the data from JGI

express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")

express_all$levelblock <-paste(express_all$level,express_all$Block,sep = "")
#need to remove 19
express_all<-express_all[-which(express_all$Block == "block19"),]

##reduce ----
##there is some Block/blocks where the difference is much stronger - lets limit to those where difference is greatest - these blocks can be pick out from the express table file 
express_table<-read_csv("clean_data/gene_interest_express_GC.csv")
select<-which(express_table[which(express_table$gene == "MnP" & express_table$level == "high"),]$sum_TPM-express_table[which(express_table$gene == "MnP" & express_table$level == "low"),]$sum_TPM > 1000)
block<-as.data.frame(express_table[which(express_table$gene == "MnP" & express_table$level == "high"),]$Block)
block[c(select),]

express_all<-express_all[which(express_all$Block %in% block[c(select),]),]

##need to get quant file in format for DEseq
##should have a counts matrix and a metadata matrix - which the rows of meta data are in same order as the columns of count
express_wide<- express_all[,-c(1,3,5)] |>  pivot_wider(names_from = "levelblock",values_from = NumReads)
prots<-express_wide$protID
express_wide<- express_wide[,-1]
rownames(express_wide)<-prots
express_meta<-as.data.frame(colnames(express_wide))
express_meta<-separate(express_meta,col = `colnames(express_wide)`,into = c("level","block"),sep="block")
rownames(express_meta)<-colnames(express_wide)

#change to interger
express_wide <- express_wide |> mutate_if(is.numeric,as.integer)
rownames(express_wide)<-prots

#now make deseq item
dds<-DESeqDataSetFromMatrix(countData = express_wide,colData = express_meta,design = ~level)

#prefilter a little
keep<-rowSums(counts(dds)) >=10
dds<-dds[keep,]

dds$level <- relevel(dds$level, ref= "low")

#run

dds<-DESeq(dds)
res<-results(dds)
#shrink
#resLFC<-lfcShrink(dds,coef = "level_high_vs_low",type = "apeglm")
#resLFC


resOrdered <- res[order(res$pvalue),]
summary(res)

res05<- results(dds,alpha = 0.05)
summary(res05)

ressig <-subset(resOrdered,padj < 0.1)
ressig
ressig.df<-as.data.frame(ressig)

#try to match some info to protID
Coromn1_FilteredModels1_go_2025_04_22 <- read_delim("raw_data/Coromn1_FilteredModels1_go_2025-04-22.tab", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)

Coromn1_FilteredModels1_ec <- read_delim("raw_data/Coromn1_FilteredModels1_ec_2025-04-22.tab", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

goname<-aggregate(go_name ~ proteinId, Coromn1_FilteredModels1_go_2025_04_22,FUN = "paste")
goacc<-aggregate(go_acc ~ proteinId, Coromn1_FilteredModels1_go_2025_04_22,FUN = "paste")
go<-cbind(goname,goacc)
go<-go[,-3]
ressig.df$goname<-go[match(rownames(ressig.df),go$proteinId),]$go_name
ressig.df$goacc<-go[match(rownames(ressig.df),go$proteinId),]$go_acc
ressig.df$ec<-Coromn1_FilteredModels1_ec[match(rownames(ressig.df),Coromn1_FilteredModels1_ec$proteinId),]$definition
ressig.df$ecNum<-Coromn1_FilteredModels1_ec[match(rownames(ressig.df),Coromn1_FilteredModels1_ec$proteinId),]$ecNum


# Add significance column
res$significant <- res$padj < 0.1 & abs(res$log2FoldChange) > 1

# Volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value")


vsd <- vst(dds, blind = TRUE)  # Variance-stabilizing transformation
plotPCA(vsd, intgroup = "level")

# Heatmap of Sample Distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, main = "Sample Distance Heatmap")
