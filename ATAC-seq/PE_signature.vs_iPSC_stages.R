#### Analysis 2 - Figure 3 #####
#### PE signature #####

## load libraries:
library(edgeR)
library("DESeq2")
require("VennDiagram")
library(annotables)
library(RUVSeq)
### all labels in Arial! - use library(extrafont)
library(extrafont)
font_import()
y

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18/")
load("endo_comp.atac.Rdata")

## read in merged ATAC-seq counts - for PE samples, compared with iPSC samples
merged_counts=read.table("external_data/iPSC_endoderm.merged_peaks.count_table.031218.txt",h=T,sep="\t",row.names=1)
merged_peak_ann = merged_counts[,1:5]
merged_counts=merged_counts[,-c(1:5)]
## run RUV setting stages as groups in RUVg
filter <- apply(merged_counts, 1, function(x) length(x[x>10])>=3)
filtered <- merged_counts[filter,]
genes <- rownames(filtered)

ipsc_stage=gsub("_.+","",colnames(merged_counts[28:51]))
design_merged=data.frame(sample=colnames(filtered), stage = c(rep("PE",27), ipsc_stage ),
                         group=c(gsub("_d.","",colnames(filtered))[1:27],ipsc_stage))
design_merged$group[which(design_merged$group=="PE")]="A_presort"

set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(design_merged, row.names=colnames(filtered)))
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4))
plotPCA(set,  cex=1.2)

differences <- makeGroups(design_merged$group)
differences

set_ruv <- RUVs(set, genes, k=1, differences)
plotRLE(set_ruv, outline=FALSE, ylim=c(-4, 4))

col=c("darkred","darkblue","cyan","darkolivegreen1","green","darkred","coral","cyan","red","darkgreen","green","coral",
      "red","darkblue","darkgreen","darkolivegreen1","coral","red","cornflowerblue","darkgreen","red","darkblue",
      "darkolivegreen1","darkred","cyan","cornflowerblue","green","#6DA567","#6DA567","#6DA567","#7883BA","#7883BA",
      "#7883BA","#96665A","#96665A","#96665A","#F4B8B0","#F4B8B0","#F4B8B0","#755A91","#755A91","#755A91","#C15858",
      "#C15858","#C15858","#CC85B1","#CC85B1","#CC85B1","#CADAE8", "#CADAE8","#CADAE8")
pch = c(rep(16,5),rep(17,6),rep(15,5), rep(16,4),rep(17,3),rep(15,4), rep(18,24))

pdf("PCA.ATAC.iPSC_controls.pdf", family="Arial")
plotPCA(set_ruv, col=col, pch=pch, labels=F, cex=2, cex.axis=1.5, cex.lab=1.5)
legend("bottomright", c("A GFP+","A presort", "A GFP-","B GFP+","B presort", "B GFP-", "C GFP+", "C presort",
 "C GFP-","iPSC","DE","GT","PF","PE","EP","EN","BLC"), col=c("darkred","red","coral","darkblue","cornflowerblue",
 "cyan","darkgreen","green","darkolivegreen1","#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567"),
       pch=c(rep(16,9),rep(18,8)), bty="n")
dev.off()

norm_counts=normCounts(set_ruv)
cpm=cpm(norm_counts)

### calculate mean CPM - for each stage from control iPSCs and all PE samples by cell population:
##   exclude PE here 
ipsc_mean_cpm=data.frame(iPSC=apply(cpm[,grepl("iPSC",colnames(cpm))],1,mean),
                         DE=apply(cpm[,grepl("DE",colnames(cpm))],1,mean),
                         PGT=apply(cpm[,grepl("GT",colnames(cpm))],1,mean),
                         PFG=apply(cpm[,grepl("PF",colnames(cpm))],1,mean),
                         EP=apply(cpm[,grepl("EP",colnames(cpm))],1,mean),
                         EN=apply(cpm[,grepl("EN",colnames(cpm))],1,mean),
                         BLC=apply(cpm[,grepl("BLC",colnames(cpm))],1,mean)
)

PE_mean_cpm=data.frame(A_presort=apply(cpm[,grepl("A_presort",colnames(cpm))],1,mean),
                       B_presort=apply(cpm[,grepl("B_presort",colnames(cpm))],1,mean),
                       C_presort=apply(cpm[,grepl("C_presort",colnames(cpm))],1,mean),
                       A_GFPpos=apply(cpm[,grepl("A_GFPpos",colnames(cpm))],1,mean),
                       B_GFPpos=apply(cpm[,grepl("B_GFPpos",colnames(cpm))],1,mean),
                       C_GFPpos=apply(cpm[,grepl("C_GFPpos",colnames(cpm))],1,mean),
                       A_GFPneg=apply(cpm[,grepl("A_GFPneg",colnames(cpm))],1,mean),
                       B_GFPneg=apply(cpm[,grepl("B_GFPneg",colnames(cpm))],1,mean),
                       C_GFPneg=apply(cpm[,grepl("C_GFPneg",colnames(cpm))],1,mean)
)

merged_CPMs = merge(PE_mean_cpm, ipsc_mean_cpm, by="row.names")
rownames(merged_CPMs) = merged_CPMs$Row.names
merged_CPMs=merged_CPMs[,-1]
merged_CPMs=merged_CPMs[rownames(merged_counts),]

#### calculate PE signatures for each method & fraction; report overlap
### calculate CV for presort, plus, minus

### for being a PE signature gene, require:
### merged_TPMs$CV>1 --> tissue selectivity
### merged_TPMs$Z_A_presort >1 &  merged_TPMs$Z_B_presort >1 & merged_TPMs$Z_C_presort >1
# ### merged_TPMs$A_presort >= 1 & merged_TPMs$B_presort >= 1 & merged_TPMs$C_presort >=1

PE_signatures=data.frame(peak=rownames(merged_CPMs))
PE_signatures=cbind(PE_signatures, PE_mean_cpm[as.character(PE_signatures$peak),] )

library(pbapply)
for (i in 1:9){  ### for the 9 comparisons
  CV=pbapply(merged_CPMs[,c(i,10:16)], 1, function(x) sd(x)/mean(x))
  Z = pbsapply(rownames(merged_CPMs), function(x) scale(as.numeric(merged_CPMs[x,c(i,10:16)]), center=T, scale=T)[1,1])

  signature = rownames(merged_CPMs[intersect(which(merged_CPMs[,i]>=1), intersect(which(CV>1), which(Z>1)) ),])
  in_signature=rownames(merged_CPMs) %in% signature
  print(paste0(colnames(merged_CPMs)[i],": ",length(signature)))
  PE_signatures=cbind(PE_signatures, in_signature)
  colnames(PE_signatures)[ncol(PE_signatures)] = paste0(colnames(merged_CPMs)[i],"_signature")
}

# # [1] "A_presort: 1951"
# [1] "B_presort: 2731"
# [1] "C_presort: 1242"
# [1] "A_GFPpos: 2146"
# [1] "B_GFPpos: 2029"

# # [1] "C_GFPpos: 1544"
# # [1] "A_GFPneg: 1887"
# # [1] "B_GFPneg: 2530"
# # [1] "C_GFPneg: 978"

presort=rownames(PE_signatures[which(PE_signatures$A_presort_signature & PE_signatures$B_presort_signature & PE_signatures$C_presort_signature),])
write.table(merged_peak_ann[presort,1:3], file="PE_iPSC_signature.presort.bed",sep="\t",quote=F,row.names=F, col.names=F)
GFPpos=rownames(PE_signatures[which(PE_signatures$A_GFPpos_signature & PE_signatures$B_GFPpos_signature & PE_signatures$C_GFPpos_signature),])
write.table(merged_peak_ann[GFPpos,1:3], file="PE_iPSC_signature.GFPpos.bed",sep="\t",quote=F,row.names=F, col.names=F)
GFPneg=rownames(PE_signatures[which(PE_signatures$A_GFPneg_signature & PE_signatures$B_GFPneg_signature & PE_signatures$C_GFPneg_signature),])
write.table(merged_peak_ann[GFPneg,1:3], file="PE_iPSC_signature.GFPneg.bed",sep="\t",quote=F,row.names=F, col.names=F)

write.table(PE_signatures, file="PE_signatures.vs_other_diff_stages.txt",sep="\t",quote=F,row.names=F)

# 
# __END__
# ### try diff exp - PE samples against all other stages - adjusting for RUV covariate
# 
# design_merged$RUV_cov= pData(set_ruv)$W_1
# 
# PE_signatures=data.frame(peak=rownames(merged_CPMs))
# PE_signatures=cbind(PE_signatures, merged_CPMs[,1:9] )
# 
# for (i in 1:9){  ### for the 9 comparisons
#   gr=levels(design_merged$group)[i]
#   select=c(grep(gr,design$sample),28:42,46:51)
#   design_gr = design_merged[select,]
#   design_gr$Gr = rep(0,nrow(design_gr))
#   design_gr$Gr[grep(gr,design_gr$group)]<- 1
#   design_gr$Gr=factor(design_gr$Gr, levels=c("0","1"))
#   
#   dds_gr <- DESeqDataSetFromMatrix(
#     countData = merged_counts[,select],
#     colData = design_gr,
#     design = ~RUV_cov+Gr)
#   
#   dds_gr <- DESeq(dds_gr)
#   res_gr = results(dds_gr)
#   length(which(res_gr$padj < 0.01 & res_gr$log2FoldChange>=1 & merged_CPMs[,gr]>=1))	### 1942
#   write.table(res_gr, file = paste0(gr,".DOCs.DESeq2.tsv"),sep="\t",quote=F,col.names=NA)
#   
#   signature = rownames(res_gr[which(res_gr$padj < 0.01 & res_gr$log2FoldChange>=1 & merged_CPMs[,gr]>=1),])
#   in_signature=rownames(merged_CPMs) %in% signature
#   print(paste0(gr,": ",length(signature)))
#   PE_signatures=cbind(PE_signatures, in_signature)
#   colnames(PE_signatures)[ncol(PE_signatures)] = paste0(gr,"_signature")
# }
# 
# presort_degs=rownames(PE_signatures[which(PE_signatures$A_presort_signature & PE_signatures$B_presort_signature & PE_signatures$C_presort_signature),])
# write.table(merged_peak_ann[presort_degs,1:3], file="PE_iPSC_signature.DEG.presort.bed",sep="\t",quote=F,row.names=F, col.names=F)
# GFPpos_degs=rownames(PE_signatures[which(PE_signatures$A_GFPpos_signature & PE_signatures$B_GFPpos_signature & PE_signatures$C_GFPpos_signature),])
# write.table(merged_peak_ann[GFPpos_degs,1:3], file="PE_iPSC_signature.DEG.GFPpos.bed",sep="\t",quote=F,row.names=F, col.names=F)
# GFPneg_degs=rownames(PE_signatures[which(PE_signatures$A_GFPneg_signature & PE_signatures$B_GFPneg_signature & PE_signatures$C_GFPneg_signature),])
# write.table(merged_peak_ann[GFPneg_degs,1:3], file="PE_iPSC_signature.DEG.GFPneg.bed",sep="\t",quote=F,row.names=F, col.names=F)
# 
# write.table(PE_signatures, file="PE_signatures.DEG.vs_other_diff_stages.txt",sep="\t",quote=F,row.names=F)

######## Venn Diagram of the PE signatures:
PE_GFPpos_signatures <- list(A=rownames(PE_signatures[PE_signatures$A_GFPpos_signature,]), 
                             B=rownames(PE_signatures[PE_signatures$B_GFPpos_signature,]), 
                             C=rownames(PE_signatures[PE_signatures$C_GFPpos_signature,]) )
require(VennDiagram)
pdf("venn.PE_iPSC_signature.GFPpos.pdf")
venn.plot <- venn.diagram(PE_GFPpos_signatures , NULL, fill=c("red", "blue","green"), alpha=c(0.5,0.5,0.5), cex = 3, cat.fontface=4, category.names=c("Protocol A", "Protocol B", "Protocol C"), main="PE signature - presort",  scaled=T)
grid.draw(venn.plot)
dev.off()

###  it doesn't work great --> I used the http://eulerr.co/ app inputing the numbers

