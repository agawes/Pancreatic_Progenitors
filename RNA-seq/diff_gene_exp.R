#### Analysis 3: Figure 4 #######
#### differential gene expression #####

## load libraries:
library(edgeR)
library("DESeq2")
require("VennDiagram")
### all labels in Arial! - use library(extrafont)
library(extrafont)
font_import()

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_reanalysis_October18/")

counts = read.table("21.07.2017.SmartSeq2.gene.counts.tsv",h=T,row.names=1)
gene_ann = counts[,1,drop=F]
counts = counts[,20:46]
colnames(counts) = c("A_GFPneg_d1", "A_GFPneg_d2", "A_GFPneg_d3","A_GFPpos_d1","A_GFPpos_d2","A_GFPpos_d3","A_presort_d1","A_presort_d2","A_presort_d3","B_GFPneg_d1","B_GFPneg_d2","B_GFPneg_d3","B_GFPpos_d1","B_GFPpos_d2","B_GFPpos_d3","B_presort_d1","B_presort_d2","B_presort_d3","C_GFPneg_d1", "C_GFPneg_d2", "C_GFPneg_d3","C_GFPpos_d1","C_GFPpos_d2","C_GFPpos_d3","C_presort_d1","C_presort_d2","C_presort_d3")
tpms = read.table("21.07.2017.SmartSeq2.gene.tpm.tsv",h=T,row.names=1)
tpm_ann = tpms[,1,drop=F]
tpms = tpms[,20:46]
colnames(tpms)=c("A_GFPneg_d1", "A_GFPneg_d2", "A_GFPneg_d3","A_GFPpos_d1","A_GFPpos_d2","A_GFPpos_d3","A_presort_d1","A_presort_d2","A_presort_d3","B_GFPneg_d1","B_GFPneg_d2","B_GFPneg_d3","B_GFPpos_d1","B_GFPpos_d2","B_GFPpos_d3","B_presort_d1","B_presort_d2","B_presort_d3","C_GFPneg_d1", "C_GFPneg_d2", "C_GFPneg_d3","C_GFPpos_d1","C_GFPpos_d2","C_GFPpos_d3","C_presort_d1","C_presort_d2","C_presort_d3")

#### remove outlier: A_presort_d3
counts = counts[,-9]
tpms = tpms[,-9]

sort=sapply(strsplit(colnames(counts), split="_"), function(x) x[2])
diff = sapply(strsplit(colnames(counts), split="_"), function(x) x[3])
protocol=sapply(strsplit(colnames(counts), split="_"), function(x) x[1])

design = data.frame(row.names = colnames(counts), sort=sort, diff = diff, protocol=protocol)

############ differential expression between methods #############
############### each method vs remaining two - separately for each sort fraction #################
### in presort:
## Protocol A vs rest
design_A_pre = design[grep("presort",design$sort),]
design_A_pre$A = rep(0,nrow(design_A_pre))
design_A_pre$A[grep("A",design_A_pre$protocol)]<- 1
design_A_pre$A=factor(design_A_pre$A, levels=c("0","1"))

dds_A_pre <- DESeqDataSetFromMatrix(
  countData = counts[,grep("presort",design$sort)],
  colData = design_A_pre,
  design = ~diff+A)

dds_A_pre <- DESeq(dds_A_pre)
res_A_pre = results(dds_A_pre)
length(which(res_A_pre$padj < 0.05))	### 876

res_A_pre = cbind(gene_ann, res_A_pre)
head(data.frame(res_A_pre[order(res_A_pre$padj),]),n=20)
write.table(res_A_pre, file = "A_vs_rest.presort.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## B vs rest
design_B_pre = design[grep("presort",design$sort),]
design_B_pre$B = rep(0,nrow(design_B_pre))
design_B_pre$B[grep("B",design_B_pre$protocol)]<- 1
design_B_pre$B=factor(design_B_pre$B, levels=c("0","1"))

dds_B_pre <- DESeqDataSetFromMatrix(
  countData = counts[,grep("presort",design$sort)],
  colData = design_B_pre,
  design = ~diff+B)

dds_B_pre <- DESeq(dds_B_pre)
res_B_pre = results(dds_B_pre)
length(which(res_B_pre$padj < 0.05))	### 2218

res_B_pre = cbind(gene_ann, res_B_pre)
head(data.frame(res_B_pre[order(res_B_pre$padj),]),n=20)
write.table(res_B_pre, file = "B_vs_rest.presort.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## C vs rest
design_C_pre = design[grep("presort",design$sort),]
design_C_pre$C = rep(0,nrow(design_C_pre))
design_C_pre$C[grep("C",design_C_pre$protocol)]<- 1
design_C_pre$C=factor(design_C_pre$C, levels=c("0","1"))

dds_C_pre <- DESeqDataSetFromMatrix(
  countData = counts[,grep("presort",design$sort)],
  colData = design_C_pre,
  design = ~diff+C)

dds_C_pre <- DESeq(dds_C_pre)
res_C_pre = results(dds_C_pre)
length(which(res_C_pre$padj < 0.05))	### 650

res_C_pre = cbind(gene_ann, res_C_pre)
head(data.frame(res_C_pre[order(res_C_pre$padj),]),n=20)
write.table(res_C_pre, file = "C_vs_rest.presort.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

############### in GFP+: ##################
## A vs rest
design_A_plus = design[grep("GFPpos",design$sort),]
design_A_plus$A = rep(0,nrow(design_A_plus))
design_A_plus$A[grep("A",design_A_plus$protocol)]<- 1
design_A_plus$A=factor(design_A_plus$A, levels=c("0","1"))

dds_A_plus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFPpos",design$sort)],
  colData = design_A_plus,
  design = ~diff+A)

dds_A_plus <- DESeq(dds_A_plus)
res_A_plus = results(dds_A_plus)
length(which(res_A_plus$padj < 0.05))	### 1206

res_A_plus = cbind(gene_ann, res_A_plus)
head(data.frame(res_A_plus[order(res_A_plus$padj),]),n=20)
write.table(res_A_plus, file = "A_vs_rest.GFPpos.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## B vs rest
design_B_plus = design[grep("GFPpos",design$sort),]
design_B_plus$B = rep(0,nrow(design_B_plus))
design_B_plus$B[grep("B",design_B_plus$protocol)]<- 1
design_B_plus$B=factor(design_B_plus$B, levels=c("0","1"))

dds_B_plus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFPpos",design$sort)],
  colData = design_B_plus,
  design = ~diff+B)

dds_B_plus <- DESeq(dds_B_plus)
res_B_plus = results(dds_B_plus)
length(which(res_B_plus$padj < 0.05))	### 2312

res_B_plus = cbind(gene_ann, res_B_plus)
head(data.frame(res_B_plus[order(res_B_plus$padj),]),n=20)
write.table(res_B_plus, file = "B_vs_rest.GFPpos.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## C vs rest
design_C_plus = design[grep("GFPpos",design$sort),]
design_C_plus$C = rep(0,nrow(design_C_plus))
design_C_plus$C[grep("C",design_C_plus$protocol)]<- 1
design_C_plus$C=factor(design_C_plus$C, levels=c("0","1"))

dds_C_plus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFPpos",design$sort)],
  colData = design_C_plus,
  design = ~diff+C)

dds_C_plus <- DESeq(dds_C_plus)
res_C_plus = results(dds_C_plus)
length(which(res_C_plus$padj < 0.05))	### 924

res_C_plus = cbind(gene_ann, res_C_plus)
head(data.frame(res_C_plus[order(res_C_plus$padj),]),n=20)
write.table(res_C_plus, file = "C_vs_rest.GFP_plus.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

### in GFPneg:
## A vs rest
design_A_minus = design[grep("GFPneg",design$sort),]
design_A_minus$A = rep(0,nrow(design_A_minus))
design_A_minus$A[grep("A",design_A_minus$protocol)]<- 1
design_A_minus$A=factor(design_A_minus$A, levels=c("0","1"))

dds_A_minus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFPneg",design$sort)],
  colData = design_A_minus,
  design = ~diff+A)

dds_A_minus <- DESeq(dds_A_minus)
res_A_minus = results(dds_A_minus)
length(which(res_A_minus$padj < 0.05))	### 3126

res_A_minus = cbind(gene_ann, res_A_minus)
head(data.frame(res_A_minus[order(res_A_minus$padj),]),n=20)
write.table(res_A_minus, file = "A_vs_rest.GFPneg.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## B vs rest
design_B_minus = design[grep("GFPneg",design$sort),]
design_B_minus$B = rep(0,nrow(design_B_minus))
design_B_minus$B[grep("B",design_B_minus$protocol)]<- 1
design_B_minus$B=factor(design_B_minus$B, levels=c("0","1"))

dds_B_minus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFPneg",design$sort)],
  colData = design_B_minus,
  design = ~diff+B)

dds_B_minus <- DESeq(dds_B_minus)
res_B_minus = results(dds_B_minus)
length(which(res_B_minus$padj < 0.05))	### 3751

res_B_minus = cbind(gene_ann, res_B_minus)
head(data.frame(res_B_minus[order(res_B_minus$padj),]),n=20)
write.table(res_B_minus, file = "B_vs_rest.GFPneg.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## C vs rest
design_C_minus = design[grep("GFPneg",design$sort),]
design_C_minus$C = rep(0,nrow(design_C_minus))
design_C_minus$C[grep("C",design_C_minus$protocol)]<- 1
design_C_minus$C=factor(design_C_minus$C, levels=c("0","1"))

dds_C_minus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFPneg",design$sort)],
  colData = design_C_minus,
  design = ~diff+C)

dds_C_minus <- DESeq(dds_C_minus)
res_C_minus = results(dds_C_minus)
length(which(res_C_minus$padj < 0.05))	### 476

res_C_minus = cbind(gene_ann, res_C_minus)
head(data.frame(res_C_minus[order(res_C_minus$padj),]),n=20)
write.table(res_C_minus, file = "C_vs_rest.GFPneg.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)


##### GFP+ vs GFP- - for each method separately

## A - GFP+ vs GFP-
design_A_GFP = design[intersect(grep("A",design$protocol),grep("GFP",design$sort)),]

dds_A_GFP <- DESeqDataSetFromMatrix(
  countData = counts[,intersect(grep("A",design$protocol),grep("GFP",design$sort))],
  colData = design_A_GFP,
  design = ~diff+sort)

dds_A_GFP <- DESeq(dds_A_GFP)
res_A_GFP = results(dds_A_GFP)
length(which(res_A_GFP$padj < 0.05))	### 5355

res_A_GFP = cbind(gene_ann, res_A_GFP)
head(data.frame(res_A_GFP[order(res_A_GFP$padj),]),n=20)
write.table(res_A_GFP, file = "A.GFPpos_v_neg.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## B - GFP+ vs GFP-
design_B_GFP = design[intersect(grep("B",design$protocol),grep("GFP",design$sort)),]

dds_B_GFP <- DESeqDataSetFromMatrix(
  countData = counts[,intersect(grep("B",design$protocol),grep("GFP",design$sort))],
  colData = design_B_GFP,
  design = ~diff+sort)

dds_B_GFP <- DESeq(dds_B_GFP)
res_B_GFP = results(dds_B_GFP)
length(which(res_B_GFP$padj < 0.05))	### 718

res_B_GFP = cbind(gene_ann, res_B_GFP)
head(data.frame(res_B_GFP[order(res_B_GFP$padj),]),n=20)
write.table(res_B_GFP, file = "B.GFPpos_v_neg.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## C - GFP+ vs GFP-
design_C_GFP = design[intersect(grep("C",design$protocol),grep("GFP",design$sort)),]

dds_C_GFP <- DESeqDataSetFromMatrix(
  countData = counts[,intersect(grep("C",design$protocol),grep("GFP",design$sort))],
  colData = design_C_GFP,
  design = ~diff+sort)

dds_C_GFP <- DESeq(dds_C_GFP)
res_C_GFP = results(dds_C_GFP)
length(which(res_C_GFP$padj < 0.05))	### 1707

res_C_GFP = cbind(gene_ann, res_C_GFP)
head(data.frame(res_C_GFP[order(res_C_GFP$padj),]),n=20)
write.table(res_C_GFP, file = "C.GFPpos_v_neg.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

