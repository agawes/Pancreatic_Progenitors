#### Jennings data - differential expression analysis #####
library(DESeq2)

jen=read.table("external_data/Jennings_2017/Jennings.count_table.txt",h=T,sep="\t",row.names=1)
jen_ann = jen[,c(1,8)]
jen=jen[,c(2:7)]

design_jen = data.frame(row.names = colnames(jen), tissue=c("DP","DP","HBP","HBP","HC","HC"),dp_vs_rest=factor(c(rep(1,2),rep(0,4))), hbp_vs_rest=factor(c(rep(0,2),rep(1,2),rep(0,2))), hc_vs_rest=factor(c(rep(0,4),rep(1,2))))

## DP vs rest
dds_DP <- DESeqDataSetFromMatrix(
  countData = round(jen),
  colData = design_jen,
  design = ~dp_vs_rest)

dds_DP <- DESeq(dds_DP)
res_DP = results(dds_DP)
length(which(res_DP$padj < 0.05))	### 482

res_DP = cbind(jen_ann, res_DP)
write.table(res_DP, file = "external_data/Jennings_2017/DP_vs_rest.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)
write.table(res_DP[which(res_DP$padj < 0.05 & res_DP$log2FoldChange>0),]$gene.name, file="external_data/Jennings_2017/DP.txt",quote=F, col.names=F, row.names=F)

# HBP vs rest
dds_HBP <- DESeqDataSetFromMatrix(
  countData = round(jen),
  colData = design_jen,
  design = ~hbp_vs_rest)

dds_HBP <- DESeq(dds_HBP)
res_HBP = results(dds_HBP)
length(which(res_HBP$padj < 0.05))	### 8

res_HBP = cbind(jen_ann, res_HBP)
write.table(res_HBP, file = "external_data/Jennings_2017/HBP_vs_rest.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)
write.table(res_HBP[which(res_HBP$padj < 0.05 & res_HBP$log2FoldChange>0),]$gene.name, file="external_data/Jennings_2017/HBP.txt",quote=F, col.names=F, row.names=F)

# HC vs rest
dds_HC <- DESeqDataSetFromMatrix(
  countData = round(jen),
  colData = design_jen,
  design = ~hc_vs_rest)

dds_HC <- DESeq(dds_HC)
res_HC = results(dds_HC)
length(which(res_HC$padj < 0.05))	### 903

res_HC = cbind(jen_ann, res_HC)
write.table(res_HC, file = "external_data/Jennings_2017/HC_vs_rest.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)
write.table(res_HC[which(res_HC$padj < 0.05 & res_HC$log2FoldChange>0),]$gene.name, file="external_data/Jennings_2017/HC.txt",quote=F, col.names=F, row.names=F)

