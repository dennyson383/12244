library(tximportData)
library(tximport)
library(readr)
library(ggplot2)
library (dplyr)
library (DESeq2)
library(plyr)
library(reshape2)
library(tidyverse)
library(biomaRt)
library(reshape2)
library(RColorBrewer) 
library(rhdf5)

setwd('~/Desktop/ZIKVmeta/Results/GSE125558/')
getwd()

# Analise de expressao difenrencial/ importar dataset e criar fatores
condition <- factor(c(rep("infected",3), rep("control", 1)))
sampleTable <- data.frame(condition = as.factor(condition))
deseq=read.table("CountsGSE125558.txt", sep = "\t", header = TRUE, row.names = 1)
rownames(sampleTable) <- colnames(deseq)
sampleTable
deseq <- DESeqDataSetFromMatrix(countData = deseq,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
DGE.results <- results(d.deseq,
                       independentFiltering = TRUE,
                       alpha = 0.05)
table(DGE.results$padj < 0.05)
#organizar tabela
resOrdered <- DGE.results[order(DGE.results$padj),]
head(resOrdered)
#padj
resSig <- subset(resOrdered, padj < 0.05)
head(resSig)
#up e down
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 1.0) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -1.0)
#export
#genes down
resOrdeger <- resl2fc_down[order(resl2fc_down$padj),]
head(resOrdeger)
ressig <- subset(resl2fc_down, padj < 0,05)
write.csv(as.data.frame(resl2fc_down), '~/Desktop/ZIKVmeta/Results/GSE125558/Comparation_down')
#export
#genesUp
resOrdeger <- resl2fc_up[order(resl2fc_up$padj),]
head(resOrdeger)
ressig <- subset(resl2fc_up, padj < 0,05)
write.csv(as.data.frame(resl2fc_up), '~/Desktop/ZIKVmeta/Results/GSE125558/Comparation_up')




