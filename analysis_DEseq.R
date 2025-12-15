library(readr)
library(tidyverse)
library(DESeq2)
library(apeglm)
library(pheatmap)

#read in the expression data and the data from JGI

express_all<-read_csv("clean_data/expression_salmon_quant_all.csv")


##need to get quant file in format for DEseq
##should have a counts matrix and a metadata matrix - which the rows of meta data are in same order as the columns of count

express_all$levelblock <-paste(express_all$level,express_all$Block,sep = "")
#need to remove low19
express_all<-express_all[-which(express_all$Block == "block19"),]

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
