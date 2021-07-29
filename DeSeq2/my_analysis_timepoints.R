library(tximportData)
library(tximport)
BiocManager::install("tximport")
library(readr)
library(tximportData)
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


#############################################################
#############################################################
#Deseq2

#analisa em ordem alfabética os valores de ref ouuu:
#dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
#ouuuu
#dds$condition <- relevel(dds$condition, ref = "untreated")

library("ggplot2")
library ("DESeq2")
condition <- factor(c(rep("Th0",27), rep("Th17",27)))
timepoints <- factor(c(rep("0.5h",3), rep("1h", 3), rep("2h", 3), rep("4h", 3), rep("6h", 3), rep("12h", 3), rep("24h", 3), rep("48h", 3), rep("72h", 3)))
sampleTable <- data.frame(condition = as.factor(condition), timepoints = as.factor(timepoints))

sampleTable$timepoints <- factor(sampleTable$timepoints, levels = c("0.5h","1h","2h","4h","6h","12h","24h","48h","72h"))

countdata=read.table("GSE52260_rc_mouse_samples_to_test.txt", sep="", header = TRUE, row.names= 1)
rownames(sampleTable) <- colnames(countdata)
sampleTable$condition
sampleTable$timepoints
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~ condition + timepoints + condition:timepoints)

deseq
ddsTP <- DESeq(deseq, test="LRT", reduced = ~ condition + timepoints)
resTP <- results(ddsTP)
resTPbackup <- resTP

write.csv(as.data.frame(resTPbackup),
          file="transcripts_Timepoints_1_teste.csv")

resTP$X <- mcols(ddsTP)$X
head(resTP[order(resTP$padj),],4)

data <- plotCounts(ddsTP, gene = c("Lnx2","Ccdc134"),
                  intgroup=c("timepoints","condition"), returnData=TRUE) #which.min(resTP$padj),
ggplot(data, aes(x=timepoints, y=count, color=condition, group=condition)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()



#resTP$X <- mcols(ddsTP)$X
#head(resTC[order(resTC$padj),],4)
#organizar tabela
resOrdered <- resTPbackup[order(resTPbackup$padj),]
head(resOrdered)
#padj
resSig <- subset(resOrdered, padj < 0.05)
head(resSig)
#up e down
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 1.0) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -1.0)
#export
write.csv(as.data.frame(resSig),
          file="significat_transcripts.csv")
write.csv(as.data.frame(resl2fc_up),
          file="transcripts_up.csv")
write.csv(as.data.frame(resl2fc_down),
          file="transcripts_down.csv")


###################

d.deseq<-DESeq(deseq)
d.deseq
res_deseq<- results(d.deseq)
res_deseq

#organizar tabela
resOrdered <- res_deseq[order(res_deseq$padj),]
head(resOrdered)
#padj
resSig <- subset(resOrdered, padj < 0.05)
head(resSig)
#up e down
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 1.0) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -1.0)
#export
write.csv(as.data.frame(resSig),
          file="significat_transcripts_k048.csv")
write.csv(as.data.frame(resl2fc_up),
          file="transcripts_up_k048.csv")
write.csv(as.data.frame(resl2fc_down),
          file="transcripts_down_k048.csv")

#Convert gene id
library(MAGeCKFlute)
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#up
table=read.table("transcripts_up_k048.csv", sep=",", header = TRUE)
transcript_ids <- table$X
transcript_ids
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "X", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[!duplicated(e$ensembl_gene_id),,drop=FALSE]

#Exportar
write.csv(as.data.frame(f),
          file="K048_up.csv")

#down
table=read.table("transcripts_down_k048.csv", sep=",", header = TRUE)
transcript_ids <- table$X
transcript_ids
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "X", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[!duplicated(e$ensembl_gene_id),,drop=FALSE]

#Exportar
write.csv(as.data.frame(f),
          file="K048_down.csv")


options(error=recover)

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

logcounts <- log2(counts(d.deseq, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts ) )
plot(hclust(dist(t(logcounts))), labels=colData(d.deseq)$condition)
plot(logcounts[,1], logcounts[,2], cex=.1)
plotMA(d.deseq, ylim=c(-5,5))
plotMA(res_deseq, ylim=c(-5,5))

var_genes <- apply(logcounts, 1,)
head(var_genes)

dds_zkifnpos_cd3 <- DESeqDataSetFromTximport(txi_kallisto_Syl_Vac, sample_info_Syl_Vac, ~treatment)
dds_zkifnpos_cd3 <- read.delim("zkifnpos_cd3.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dds_zkifnpos_cd3 <- DESeq(dds_zkifnpos_cd3)
res_zkifnpos_cd3 <- results(dds_zkifnpos_cd3)

#PCA plot#
vsdB <- varianceStabilizingTransformation(d.deseq)

plotPCA(vsdB, intgroup=c("condition"))
#plotPCA(vsdB, intgroup=c("contidion","timepoints"))

#heatmap
library(gplots)
vsdB_table<-as.data.frame(assay(vsdB))
vsdB_table_rowsum<-transform(vsdB_table, sum=rowSums(vsdB_table))
colnames(vsdB_table_rowsum)
selected<-order(vsdB_table_rowsum$sum, decreasing = TRUE) [100:160]
vsdB_table[selected,]
heatmap.2(as.matrix(vsdB_table[selected,]), Rowv=F, dendrogram = "column", col = "heat.colors", 
          density.info ="none", trace="none", cexCol=0.9, labRow=NA) 

mypalette <- brewer.pal(11, 'RdYlBu')
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple", "orange")[sampleTable$condition]
heatmap.2(vsdB_table, col=rev(morecols(50)), trace = 'none',
          main = "Top 500 most variable genes across samples",
          ColSideColors =col.cell,
          scale = 'row')
#





DGE.results <- results(dds_zkifnpos_cd3,
                       independentFiltering = TRUE,
                       alpha = 0.05)
table(DGE.results$padj < 0.05)



#EXTRA
#dds_Syl_Vac <- estimateSizeFactors( dds_Syl_Vac )
#sizeFactors(dds_Syl_Vac)
#colSums(counts(dds_Syl_Vac))
logcounts <- log2( counts(dds_zkifnpos_cd3, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts ) )
plot(hclust(dist(t(logcounts))), labels=colData(dds_Syl_Vac)$protocol)
plot(hclust(dist(t(logcounts))), labels=colData(dds_Syl_Vac)$treatment)
plot(logcounts[,1], logcounts[,2], cex=.1)
plotMA(dds_Syl_Vac, ylim=c(-5,5))
plotMA(res_Syl_Vac, ylim=c(-5,5))