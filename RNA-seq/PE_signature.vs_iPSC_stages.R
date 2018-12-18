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

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_reanalysis_October18/")
load("endocomp.Rdata")

## read in RNA-seq counts from control cell lines
ipsc_counts= read.table("external_data/31.01.2017.Differentiation_v2.gene.counts.tsv",h=T,sep="\t",row.names=1)
ipsc_gene_ann=ipsc_counts[,1,drop=F]
ipsc_counts=ipsc_counts[,-1]
ipsc_stage = gsub(".+\\..","",colnames(ipsc_counts))

which(rownames(ipsc_counts) != rownames(counts))
merged_counts=cbind(counts, ipsc_counts)

## run RUV setting stages as groups in RUVg
filter <- apply(merged_counts, 1, function(x) length(x[x>5])>=3)
filtered <- merged_counts[filter,]
genes <- rownames(filtered)
design_merged=data.frame(sample=colnames(filtered), stage = c(rep("PE",26), ipsc_stage ),
                         group=c(gsub("_d.","",colnames(filtered))[1:26],ipsc_stage))
design_merged$group[which(design_merged$group=="PE")]="A_presort"

set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(design_merged, row.names=colnames(filtered)))
set <- betweenLaneNormalization(set, which="upper")
# plotRLE(set, outline=FALSE, ylim=c(-4, 4))
# plotPCA(set,  cex=1.2)

differences <- makeGroups(design_merged$group)
# differences

### better: set only PE and A GFP presort as group, rest to each own.

set_ruv <- RUVs(set, genes, k=1, differences)
# plotRLE(set_ruv, outline=FALSE, ylim=c(-4, 4))
# plotPCA(set_ruv,  cex=1.2)

col=c(rep("coral",3),rep("darkred",3),rep("red",2), rep("cyan",3),rep("darkblue",3),rep("cornflowerblue",3),rep("darkolivegreen1",3),
      rep("darkgreen",3), rep("green",3), "#7883BA","#96665A","#6DA567","#F4B8B0","#CADAE8","#C15858","#CC85B1","#755A91","#7883BA",
      "#96665A","#6DA567","#F4B8B0","#CADAE8","#C15858","#CC85B1","#755A91","#7883BA","#96665A","#6DA567","#F4B8B0","#CADAE8","#C15858",
      "#CC85B1","#755A91")

pch = c(16,17,15,16,17,15,16,17,16,17,15,16,17,15,16,17,15,16,17,15,16,17,15,16,17,15, rep(18,24))

pdf("PCA.RNA.iPSC_controls.pdf", family="Arial")
par(mar=c(5.1,5.1,1,1))
plotPCA(set_ruv, col=col, pch=pch, labels=F, cex=2, cex.axis=1.5, cex.lab=1.5)
legend("bottomleft", c("A GFP+","A presort", "A GFP-","B GFP+","B presort", "B GFP-", "C GFP+", "C presort",
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
                         PGT=apply(cpm[,grepl("PGT",colnames(cpm))],1,mean),
                         PFG=apply(cpm[,grepl("PFG",colnames(cpm))],1,mean),
                         EP=apply(cpm[,grepl("EP",colnames(cpm))],1,mean),
                         EN6=apply(cpm[,grepl("ENstage6",colnames(cpm))],1,mean),
                         EN7=apply(cpm[,grepl("ENstage7",colnames(cpm))],1,mean)
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

merged_CPMs$gene = sapply(rownames(merged_CPMs), function(x) as.character(gene_ann[x,]))

# ## filter down to only protein coding genes
prot_coding=sapply(gsub("\\..+","",rownames(merged_CPMs)), function(x) grch37[which(grch37$ensgene==x)[1],"biotype"]=="protein_coding")
merged_CPMs=merged_CPMs[prot_coding,]


#### calculate PE signatures for each method & fraction; report overlap
### calculate CV for presort, plus, minus

### for being a PE signature gene, require:
### merged_TPMs$CV>1 --> tissue selectivity
### merged_TPMs$Z_A_presort >1 &  merged_TPMs$Z_B_presort >1 & merged_TPMs$Z_C_presort >1
### merged_TPMs$A_presort >= 1 & merged_TPMs$B_presort >= 1 & merged_TPMs$C_presort >=1

PE_signatures=data.frame(EnsGene=rownames(merged_CPMs), GeneName=merged_CPMs$gene)
PE_signatures=cbind(PE_signatures, PE_mean_cpm[as.character(PE_signatures$EnsGene),] )

for (i in 1:9){  ### for the 9 comparisons
  CV=apply(merged_CPMs[,c(i,10:16)], 1, function(x) sd(x)/mean(x))
  Z = sapply(rownames(merged_CPMs), function(x)
    (merged_CPMs[x,i]-mean(as.numeric(merged_CPMs[x,c(i,10:16)])))/sd(as.numeric(merged_CPMs[x,c(i,10:16)])))
  signature = merged_CPMs[intersect(which(merged_CPMs[,i]>=1), intersect(which(CV>1), which(Z>1)) ),"gene"]
  in_signature=merged_CPMs$gene %in% signature
  print(paste0(colnames(merged_CPMs)[i],": ",length(signature)))
  PE_signatures=cbind(PE_signatures, in_signature)
  colnames(PE_signatures)[ncol(PE_signatures)] = paste0(colnames(merged_CPMs)[i],"_signature")
}

# [1] "A_presort: 281"
# [1] "B_presort: 727"
# [1] "C_presort: 629"
# [1] "A_GFPpos: 405"
# [1] "B_GFPpos: 634"
# [1] "C_GFPpos: 498"
# [1] "A_GFPneg: 571"
# [1] "B_GFPneg: 769"
# [1] "C_GFPneg: 651" 
presort=as.character(PE_signatures[which(PE_signatures$A_presort_signature & PE_signatures$B_presort_signature & PE_signatures$C_presort_signature),]$GeneName)
write.table(presort, file="PE_iPSC_signature.presort.txt",sep="\t",quote=F,row.names=F, col.names=F)
GFPpos=as.character(PE_signatures[which(PE_signatures$A_GFPpos_signature & PE_signatures$B_GFPpos_signature & PE_signatures$C_GFPpos_signature),]$GeneName)
write.table(GFPpos, file="PE_iPSC_signature.GFPpos.txt",sep="\t",quote=F,row.names=F, col.names=F)
GFPneg=as.character(PE_signatures[which(PE_signatures$A_GFPneg_signature & PE_signatures$B_GFPneg_signature & PE_signatures$C_GFPneg_signature),]$GeneName)
write.table(GFPneg, file="PE_iPSC_signature.GFPneg.txt",sep="\t",quote=F,row.names=F, col.names=F)

write.table(PE_signatures, file="PE_signatures.vs_other_diff_stages.txt",sep="\t",quote=F,row.names=F)
PE_signatures=read.table("PE_signatures.vs_other_diff_stages.txt",sep="\t",h=T)

######## Venn Diagram of the PE signatures:
PE_GFPpos_signatures <- list(A=as.character(PE_signatures[PE_signatures$A_GFPpos_signature,]$GeneName),
                             B=as.character(PE_signatures[PE_signatures$B_GFPpos_signature,]$GeneName),
                             C=as.character(PE_signatures[PE_signatures$C_GFPpos_signature,]$GeneName) )
require(VennDiagram)
pdf("venn.PE_iPSC_signature.GFPpos.pdf")
venn.plot <- venn.diagram(PE_GFPpos_signatures , NULL, fill=c("red", "blue","green"), alpha=c(0.5,0.5,0.5), cex = 3, cat.fontface=4, category.names=c("Protocol A", "Protocol B", "Protocol C"), main="PE signature - presort",  scaled=T)
grid.draw(venn.plot)
dev.off()

# ###  it doesn't work great --> I used the http://eulerr.co/ app inputing the numbers


__END__
### this was another attempt, but doesn't work as well, e.g. ONECUT1 is not DE, even though it picks up at PE


design_merged$RUV_cov= pData(set_ruv)$W_1

PE_signatures=data.frame(peak=rownames(merged_CPMs))
PE_signatures=cbind(PE_signatures, merged_CPMs[,1:9] )

for (i in 1:9){  ### for the 9 comparisons
  gr=levels(design_merged$group)[i]
  select=c(grep(gr,design_merged$sample),27:31,33:39, 41:47,49:50)
  design_gr = design_merged[select,]
  design_gr$Gr = rep(0,nrow(design_gr))
  design_gr$Gr[grep(gr,design_gr$group)]<- 1
  design_gr$Gr=factor(design_gr$Gr, levels=c("0","1"))
  
  dds_gr <- DESeqDataSetFromMatrix(
    countData = merged_counts[,select],
    colData = design_gr,
    design = ~RUV_cov+Gr)
  
  dds_gr <- DESeq(dds_gr)
  res_gr = results(dds_gr)
  length(which(res_gr$padj < 0.01 & res_gr$log2FoldChange>=1 & merged_CPMs[,gr]>=1))	### 1942
  write.table(res_gr, file = paste0(gr,".DOCs.DESeq2.tsv"),sep="\t",quote=F,col.names=NA)
  
  signature = rownames(res_gr[which(res_gr$padj < 0.01 & res_gr$log2FoldChange>=1 & merged_CPMs[,gr]>=1),])
  in_signature=rownames(merged_CPMs) %in% signature
  print(paste0(gr,": ",length(signature)))
  PE_signatures=cbind(PE_signatures, in_signature)
  colnames(PE_signatures)[ncol(PE_signatures)] = paste0(gr,"_signature")
}

presort_degs=rownames(PE_signatures[which(PE_signatures$A_presort_signature & PE_signatures$B_presort_signature & PE_signatures$C_presort_signature),])
write.table(gene_ann[presort_degs,], file="PE_iPSC_signature.DEG.presort.txt",sep="\t",quote=F,row.names=F, col.names=F)
GFPpos_degs=rownames(PE_signatures[which(PE_signatures$A_GFPpos_signature & PE_signatures$B_GFPpos_signature & PE_signatures$C_GFPpos_signature),])
write.table(gene_ann[GFPpos_degs,], file="PE_iPSC_signature.DEG.GFPpos.txt",sep="\t",quote=F, col.names=F,row.names=F)
GFPneg_degs=rownames(PE_signatures[which(PE_signatures$A_GFPneg_signature & PE_signatures$B_GFPneg_signature & PE_signatures$C_GFPneg_signature),])
write.table(gene_ann[GFPneg_degs,], file="PE_iPSC_signature.DEG.GFPneg.txt",sep="\t",quote=F, col.names=F,row.names=F)

write.table(PE_signatures, file="PE_signatures.DEG.vs_other_diff_stages.txt",sep="\t",quote=F,row.names=F)

