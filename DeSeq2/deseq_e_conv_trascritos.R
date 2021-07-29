# Selecione as linhas que quer comentar e aperte ctrl+shift+c para coment?-las ou descoment?-las
# Os c?digos a baixo s?o para mostrar a diferen?a entre a abordagem tradicional e a com shrinkage
# Dados shrink tendem a representar melhor o que ter?amos caso nossas quantidades de amostras fossem maiores
# resNORM <- lfcShrink(res, coef = 2, type="normal")
# resAPE <- lfcShrink(res, coef= 2, type="apeglm")
# resASH <- lfcShrink(res, coef= 2, type="ashr")
# res_deseq <- resAPE
# 
# plotMA(res, ylim=c(-10,10))
# plotMA(res, ylim=c(-2,2))
# 
# plotMA(resNORM, ylim=c(-10,10))
# plotMA(resNORM, ylim=c(-2,2))
# 
# plotMA(resAPE, ylim=c(-10,10))
# plotMA(resAPE, ylim=c(-2,2))
# 
# plotMA(resASH, ylim=c(-10,10))
# plotMA(resASH, ylim=c(-2,2))


setwd("E:/Metan?lise/Novo/GSE139181")
library (DESeq2)
library(biomaRt)
library(MAGeCKFlute)
BiocManager::install("apeglm")
# install.packages("ashr")
library(apeglm)
library (ashr)
#library(writexl)
library(openxlsx)
#install.packages("readxl")
library(readxl)
library(tximport)
library(tximportData)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#############################################
################ni_control_t1################
#############################################

countdata=read.table("gse99081.txt", sep = "\t", header = T, row.names = 1)
condition <- factor(c(rep("control",8), rep("infected",8)))
head(countdata)
sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 1.0)
resl2fc_down <- subset(resl2fc, log2FoldChange < -1.0)
write.xlsx(as.data.frame(resSig), "significat_transcripts.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up.xlsx")
table
transcript_ids <- table$...1

res <- getBM(attributes = c('ensembl_transcript_id_version',
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)

c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]

write.xlsx(as.data.frame(g),"transcripts_converted_up.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),
          "transcripts_converted_down.xlsx", row.names=TRUE)



#############################################
################ni_control_t2################
#############################################


countdata=read.table("controlvaccine.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("vaccine",3)))

sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig),
          "significat_transcripts_peripheralNeurons.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_peripheralNeurons.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_peripheralNeurons.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_peripheralNeurons.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),
          "transcripts_conv_up_ni_control_t2.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_ni_control_t2.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),
          "transcripts_conv_down_ni_control_t2.xlsx", row.names=TRUE)




#############################################
##################ni_t1_t2###################
#############################################


countdata=read.table("Macrophag.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",2), rep("infected",2)))

sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_macrophag.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_macrophag.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_macrophag.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_up_macrophag.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_macrophag.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)

c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_down_macrophag.xlsx", row.names=TRUE)


#############################################
##############immune_control_t1##############
#############################################


countdata=read.table("Microglia .txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",2), rep("infected",2)))

sampleTable <- data.frame(condition = as.factor(condition))
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95)
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_microglia.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_microglia.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_microglia.xlsx", row.names=TRUE)

#up
table=read_excel("transcripts_up_microglia.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_converted_up_microglia.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_microglia.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_converted_down_microglia.xlsx", row.names=TRUE)


#############################################
##############immune_control_t2##############
#############################################


countdata=read.table("B cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_Bcells.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_Bcells.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_Bcells.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_Bcells.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_converted_up_Bcells.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_Bcells.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_converted_down_Bcells.xlsx", row.names=TRUE)




#############################################
################immune_t1_t2#################
#############################################


countdata=read.table("CD4 T cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_CD4.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_CD4.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_CD4.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_CD4.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),"transcripts_converted_up_CD4.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_CD4.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_down_CD4.xlsx", row.names=TRUE)





#############################################
##############fetal_control_t1###############
#############################################


countdata=read.table("CD8 T cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_CD8.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_CD8.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_CD8.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_CD8.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_converted_up_CD8.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_CD8.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_down_CD8.xlsx", row.names=TRUE)



#############################################
##############fetal_control_t2###############
#############################################


countdata=read.table("NK cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_NKcells.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_NKcells.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_NKcells.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_NKcells.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_converted_up_NKcells.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_NKcells.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_down_NKcells.xlsx", row.names=TRUE)




#############################################
################fetal_t1_t2##################
#############################################


countdata=read.table("Monocyte cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_monocyte.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_monocyte.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_monocyte.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_monocyte.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
write.xlsx(as.data.frame(g),"transcripts_conv_up_monocyte.xlsx", row.names=TRUE)
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_up_monocyte.xlsx_external_and_entrez", row.names=TRUE)

#down
table=read_excel("transcripts_down_monocyte.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),
          "transcripts_conv_down_fetal_t1_t2.xlsx", row.names=TRUE)




#############################################
################fetal_t1_t2##################
#############################################


countdata=read.table("mDC cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_mDC.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_mDC.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_mDC.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_mDC.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_up_mDC.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_mDC.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),
           "transcripts_conv_down_fetal_t1_t2.xlsx", row.names=TRUE)



#############################################
################fetal_t1_t2##################
#############################################


countdata=read.table("pDC Cells.txt", sep="\t", header = TRUE, row.names= 1)
condition <- factor(c(rep("control",3), rep("infected",3)))

sampleTable <- data.frame(condition = as.factor(condition))
sampleTable$condition <- relevel(sampleTable$condition, ref = "control")
rownames(sampleTable) <- colnames(countdata)
deseq <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = sampleTable,
                                design = ~condition)
d.deseq<-DESeq(deseq)
res_deseq <- lfcShrink(d.deseq, coef= 2, type="apeglm")

resOrdered <- res_deseq[order(res_deseq$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_up <- subset(resl2fc, log2FoldChange > 0.95) 
resl2fc_down <- subset(resl2fc, log2FoldChange < -0.95)

write.xlsx(as.data.frame(resSig), "significat_transcripts_pDC.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_up), "transcripts_up_pDC.xlsx", row.names=TRUE)
write.xlsx(as.data.frame(resl2fc_down), "transcripts_down_pDC.xlsx", row.names=TRUE)


#up
table=read_excel("transcripts_up_pDC.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g), "transcripts_conv_up_pDC.xlsx", row.names=TRUE)

#down
table=read_excel("transcripts_down_pDC.xlsx")
transcript_ids <- table$...1
res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)


c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]
write.xlsx(as.data.frame(g),
           "transcripts_conv_down_pDC.xlsx", row.names=TRUE)
