#### differentially open chromatin #####
library(edgeR)
library("DESeq2")
library(sva)

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18/")
load("endo_comp.atac.Rdata")


design$diff=factor(design$diff)
design$pool=factor(design$pool)

design$A=factor(ifelse(design$protocol=="A",1,0))
design$B=factor(ifelse(design$protocol=="B",1,0))
design$C=factor(ifelse(design$protocol=="C",1,0))

design_presort = design[grep("presort",design$GFP),]
design_GFP_plus = design[grep("GFP\\+",design$GFP),]
design_GFP_minus = design[grep("GFP\\-",design$GFP),]

############### each method -vs rest-  for each sort fraction #################
################ Protocol A  ################
## presort

dds_A_pre <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$GFP=="presort")],
  colData = design_presort,
  design = ~pool+diff+A)

dds_A_pre <- DESeq(dds_A_pre)
res_A_pre = results(dds_A_pre)
length(which(res_A_pre$padj < 0.05))	### 6013
length(which(res_A_pre$padj < 0.05 & res_A_pre$log2FoldChange>0))	### 5999
length(which(res_A_pre$padj < 0.05 & res_A_pre$log2FoldChange<0))	### 14

res_A_pre = cbind(peak_ann, res_A_pre)
head(data.frame(res_A_pre[order(res_A_pre$padj),]),n=20)
write.table(res_A_pre, file = "A_vs_all.presort.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## GFP+
dds_A_plus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFP\\+",design$GFP)],
  colData = design_GFP_plus,
  design = ~pool+diff+A)

dds_A_plus <- DESeq(dds_A_plus)
res_A_plus = results(dds_A_plus)
length(which(res_A_plus$padj < 0.05))	### 10445
length(which(res_A_plus$padj < 0.05 & res_A_plus$log2FoldChange>0))	### 10307
length(which(res_A_plus$padj < 0.05 & res_A_plus$log2FoldChange<0))	### 138

res_A_plus = cbind(peak_ann, res_A_plus)
head(data.frame(res_A_plus[order(res_A_plus$padj),]),n=20)
write.table(res_A_plus, file = "A_vs_all.GFP_plus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## GFP-
dds_A_min <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFP\\-",design$GFP)],
  colData = design_GFP_minus,
  design = ~pool+diff+A)

dds_A_min <- DESeq(dds_A_min)
res_A_min = results(dds_A_min)
length(which(res_A_min$padj < 0.05))	### 11223
length(which(res_A_min$padj < 0.05 & res_A_min$log2FoldChange>0))	### 10829
length(which(res_A_min$padj < 0.05 & res_A_min$log2FoldChange<0))	### 394

res_A_min = cbind(peak_ann, res_A_min)
head(data.frame(res_A_min[order(res_A_min$padj),]),n=20)
write.table(res_A_min, file = "A_vs_all.GFP_minus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)


######################## Protocol B ########################
## presort

dds_B_pre <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$GFP=="presort")],
  colData = design_presort,
  design = ~pool+diff+B)

dds_B_pre <- DESeq(dds_B_pre)
res_B_pre = results(dds_B_pre)
length(which(res_B_pre$padj < 0.05))	### 0
length(which(res_B_pre$padj < 0.05 & res_B_pre$log2FoldChange>0))	### 0
length(which(res_B_pre$padj < 0.05 & res_B_pre$log2FoldChange<0))	### 0

res_B_pre = cbind(peak_ann, res_B_pre)
head(data.frame(res_B_pre[order(res_B_pre$padj),]),n=20)
write.table(res_B_pre, file = "B_vs_all.presort.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## GFP+

dds_B_plus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFP\\+",design$GFP)],
  colData = design_GFP_plus,
  design = ~pool+diff+B)

dds_B_plus <- DESeq(dds_B_plus)
res_B_plus = results(dds_B_plus)
length(which(res_B_plus$padj < 0.05))	### 83
length(which(res_B_plus$padj < 0.05 & res_B_plus$log2FoldChange>0))	### 83
length(which(res_B_plus$padj < 0.05 & res_B_plus$log2FoldChange<0))	### 0

res_B_plus = cbind(peak_ann, res_B_plus)
head(data.frame(res_B_plus[order(res_B_plus$padj),]),n=20)
write.table(res_B_plus, file = "B_vs_all.GFP_plus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## GFP-
dds_B_min <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFP\\-",design$GFP)],
  colData = design_GFP_minus,
  design = ~pool+diff+B)

dds_B_min <- DESeq(dds_B_min)
res_B_min = results(dds_B_min)
length(which(res_B_min$padj < 0.05))	### 271
length(which(res_B_min$padj < 0.05 & res_B_min$log2FoldChange>0))	### 246
length(which(res_B_min$padj < 0.05 & res_B_min$log2FoldChange<0))	### 25

res_B_min = cbind(peak_ann, res_B_min)
head(data.frame(res_B_min[order(res_B_min$padj),]),n=20)
write.table(res_B_min, file = "B_vs_all.GFP_minus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)


#################### Protocol C ########################
## presort

dds_C_pre <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$GFP=="presort")],
  colData = design_presort,
  design = ~pool+diff+C)

dds_C_pre <- DESeq(dds_C_pre)
res_C_pre = results(dds_C_pre)
length(which(res_C_pre$padj < 0.05))	### 117
length(which(res_C_pre$padj < 0.05 & res_C_pre$log2FoldChange>0))	### 95
length(which(res_C_pre$padj < 0.05 & res_C_pre$log2FoldChange<0))	### 22

res_C_pre = cbind(peak_ann, res_C_pre)
head(data.frame(res_C_pre[order(res_C_pre$padj),]),n=20)
write.table(res_C_pre, file = "C_vs_all.presort.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## GFP+

dds_C_plus <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFP\\+",design$GFP)],
  colData = design_GFP_plus,
  design = ~pool+diff+C)

dds_C_plus <- DESeq(dds_C_plus)
res_C_plus = results(dds_C_plus)
length(which(res_C_plus$padj < 0.05))	### 18
length(which(res_C_plus$padj < 0.05 & res_C_plus$log2FoldChange>0))	### 17
length(which(res_C_plus$padj < 0.05 & res_C_plus$log2FoldChange<0))	### 1

res_C_plus = cbind(peak_ann, res_C_plus)
head(data.frame(res_C_plus[order(res_C_plus$padj),]),n=20)
write.table(res_C_plus, file = "C_vs_all.GFP_plus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## in GFP-
dds_C_min <- DESeqDataSetFromMatrix(
  countData = counts[,grep("GFP\\-",design$GFP)],
  colData = design_GFP_minus,
  design = ~pool+diff+C)

dds_C_min <- DESeq(dds_C_min)
res_C_min = results(dds_C_min)
length(which(res_C_min$padj < 0.05))	### 106
length(which(res_C_min$padj < 0.05 & res_C_min$log2FoldChange>0))	### 106
length(which(res_C_min$padj < 0.05 & res_C_min$log2FoldChange<0))	### 0

res_C_min = cbind(peak_ann, res_C_min)
head(data.frame(res_C_min[order(res_C_min$padj),]),n=20)
write.table(res_C_min, file = "C_vs_all.GFP_minus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)


##### GFP+ vs GFP- - for each method separately
design$sort=as.character(design$GFP)
design$sort[design$sort=="GFP+"] <- "GFP_plus"
design$sort[design$sort=="GFP-"] <- "GFP_minus"
design$sort=factor(design$sort)

## A - GFP+ vs GFP-
design_A_GFP = design[intersect(grep("A",design$protocol),grep("GFP",design$GFP)),]

dds_A_GFP <- DESeqDataSetFromMatrix(
  countData = counts[,intersect(grep("A",design$protocol),grep("GFP",design$GFP))],
  colData = design_A_GFP,
  design = ~pool+diff+sort)

dds_A_GFP <- DESeq(dds_A_GFP)
res_A_GFP = results(dds_A_GFP)
length(which(res_A_GFP$padj < 0.05))	### 4083
length(which(res_A_GFP$padj < 0.05 & res_A_GFP$log2FoldChange>0))	### 1165
length(which(res_A_GFP$padj < 0.05 & res_A_GFP$log2FoldChange<0))	### 2918

res_A_GFP = cbind(peak_ann, res_A_GFP)
head(data.frame(res_A_GFP[order(res_A_GFP$padj),]),n=20)
write.table(res_A_GFP, file = "A.GFP_plus_vs_minus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## B - GFP+ vs GFP-
design_B_GFP = design[intersect(grep("B",design$protocol),grep("GFP",design$GFP)),]

dds_B_GFP <- DESeqDataSetFromMatrix(
  countData = counts[,intersect(grep("B",design$protocol),grep("GFP",design$GFP))],
  colData = design_B_GFP,
  design = ~pool+diff+sort)

dds_B_GFP <- DESeq(dds_B_GFP)
res_B_GFP = results(dds_B_GFP)
length(which(res_B_GFP$padj < 0.05))	### 0

res_B_GFP = cbind(peak_ann, res_B_GFP)
head(data.frame(res_B_GFP[order(res_B_GFP$padj),]),n=20)
write.table(res_B_GFP, file = "B.GFP_plus_vs_minus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## C - GFP+ vs GFP-
design_C_GFP = design[intersect(grep("C",design$protocol),grep("GFP",design$GFP)),]

dds_C_GFP <- DESeqDataSetFromMatrix(
  countData = counts[,intersect(grep("C",design$protocol),grep("GFP",design$GFP))],
  colData = design_C_GFP,
  design = ~pool+diff+sort)

dds_C_GFP <- DESeq(dds_C_GFP)
res_C_GFP = results(dds_C_GFP)
length(which(res_C_GFP$padj < 0.05))	### 0

res_C_GFP = cbind(peak_ann, res_C_GFP)
head(data.frame(res_C_GFP[order(res_C_GFP$padj),]),n=20)
write.table(res_C_GFP, file = "C.GFP_plus_vs_minus.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)



