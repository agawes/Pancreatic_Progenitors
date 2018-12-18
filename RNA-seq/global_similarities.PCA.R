#### Analysis 1: Figure 2 #######
#### global PCA analysis ####
#### genes expressed at 1 TPMs ##### 

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

############################ similarities between protocols  #######################################

############################## Global view of similarities and differences ##############################
#####  MDS and PCA plots
sort=sapply(strsplit(colnames(counts), split="_"), function(x) x[2])
diff = sapply(strsplit(colnames(counts), split="_"), function(x) x[3])
protocol=sapply(strsplit(colnames(counts), split="_"), function(x) x[1])

design = data.frame(row.names = colnames(counts), sort=sort, diff = diff, protocol=protocol)
design$protocol_sort = factor(paste(design$protocol, design$sort, sep="_"))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = design,
  design = ~1)

vsd <- varianceStabilizingTransformation(dds) 
vstMat = assay(vsd)

### col by method: A - red, B - blue, C - green
### col intensity - sort: darkest - GFP+, lightest - GFP-, med - presort

# col by sort (captures both subject and line), different reps will be the same col
# GFP-  red, GFP+  blue, presort - green
# coral, red, darkred
# cyan, cornflowerblue, darkblue
# chartreuse, green, darkgreen

design$protocol_sort = factor(paste(design$protocol, design$sort, sep="_"))

col_sort = c("coral","darkred","red","cyan","darkblue","cornflowerblue","darkolivegreen1","darkgreen","green")
col = sapply(design$protocol_sort, function(x) col_sort[which(levels(design$protocol_sort) == x)])

pch_diff=c(16,17,15)
pch = sapply(design$diff, function(x) pch_diff[which(levels(design$diff) == x)])

labels=c("A: GFP+","A: GFP-", "A: presort", "B: GFP+","B: GFP-","B: presort", "C: GFP+","C: GFP-","C: presort")

save.image("endocomp.Rdata")


pdf("Fig2B.RNA_PCA.pdf")
par(mar=c(5.1,5.2,1,1), xpd=T)
plotMDS(vstMat, col = col, pch = pch, cex=2, gene.selection="common",cex.lab=1.5, cex.axis=1.5)

par(mar=c(5.1,5.1,1,10), xpd=T)
plotMDS(vstMat, col = col, pch = pch, cex=2, gene.selection="common",cex.lab=1.5, cex.axis=1.5)
legend(1.75,1.1,labels,col=col_sort[c(2,1,3,5,4,6,8,7,9)], pch=16,bty="n", title="Protocol", cex=1.5)
legend(1.75,-0.5,c(1,2,3),pch=pch_diff,bty="n",title="Diff", cex=1.5)
dev.off()


##### Venn Diagrams of genes expressed at 1TPM

A_presort_mean = apply(tpms[,grep("A_presort",colnames(tpms))], 1, mean)
B_presort_mean = apply(tpms[,grep("B_presort",colnames(tpms))], 1, mean) 
C_presort_mean = apply(tpms[,grep("C_presort",colnames(tpms))], 1, mean) 

### Venn diagrams of genes expressed at >=1 TPM - separately for presort, GFP+, GFP-

######## method vs rest, by fraction
presort_tpm1 <- list(A=as.character(gene_ann[which(A_presort_mean>=1),]), 
                     C=as.character(gene_ann[which(C_presort_mean>=1),]), 
                     B=as.character(gene_ann[which(B_presort_mean>=1),]))
pdf("venn.presort.tpm1.pdf", family="Arial")
venn.plot <- venn.diagram(presort_tpm1 , NULL, fill=c("darkgreen", "darkblue","yellow"), alpha=c(0.5,0.5,0.5), 
                          cex = 2, cat.fontface=4, category.names=c("A", "C", "B"), main="Presort: tpm>=1",  scaled=T)
grid.draw(venn.plot)
dev.off()

A_plus_mean = apply(tpms[,grep("A_GFPpos",colnames(tpms))], 1, mean)
B_plus_mean = apply(tpms[,grep("B_GFPpos",colnames(tpms))], 1, mean)
C_plus_mean = apply(tpms[,grep("C_GFPpos",colnames(tpms))], 1, mean)

plus_tpm1 <- list(A=as.character(gene_ann[which(A_plus_mean>=1),]), 
                  C=as.character(gene_ann[which(C_plus_mean>=1),]), 
                  B=as.character(gene_ann[which(B_plus_mean>=1),]))
pdf("venn.GFP_plus.tpm1.pdf")
venn.plot <- venn.diagram(plus_tpm1 , NULL, fill=c("darkgreen", "darkblue","yellow"), alpha=c(0.5,0.5,0.5), cex = 2, 
                          cat.fontface=4, category.names=c("A", "C", "B"), main="GFP+ - tpm>=10",  scaled=T)
grid.draw(venn.plot)
dev.off()

A_minus_mean = apply(tpms[,grep("A_GFPneg",colnames(tpms))], 1, mean)
B_minus_mean = apply(tpms[,grep("B_GFPneg",colnames(tpms))], 1, mean)
C_minus_mean = apply(tpms[,grep("C_GFPneg",colnames(tpms))], 1, mean)

minus_tpm1 <- list(A=as.character(gene_ann[which(A_minus_mean>=1),]), 
                   C=as.character(gene_ann[which(C_minus_mean>=1),]), 
                   B=as.character(gene_ann[which(B_minus_mean>=1),]))
pdf("venn.GFP_minus.tpm1.pdf")
venn.plot <- venn.diagram(minus_tpm1 , NULL, fill=c("darkgreen", "darkblue","yellow"), alpha=c(0.5,0.5,0.5), cex = 2,
                          cat.fontface=4, category.names=c("A", "C", "B"), main="GFP- TPM>=1",  scaled=T)
grid.draw(venn.plot)
dev.off()

#### PCA plot together with data from Ramond

counts_Ramond=read.table("external_data/11.01.2017.sorted_cells.gene.counts.tsv",h=T,sep="\t",row.names=1)
gene_ann_Ramond = counts_Ramond[,1,drop=F]
counts_Ramond=counts_Ramond[,-c(1,14:ncol(counts_Ramond))]

counts_merged = merge(counts, counts_Ramond, by="row.names")
rownames(counts_merged) = counts_merged$Row.names
counts_merged=counts_merged[,-1]

## run vstMat together and plotMDS:
design_merged = data.frame(row.names = colnames(counts_merged), 
                        sort=c(as.character(design$sort),rep(c("A","B","C","D"),3)), 
                        diff = c(as.character(design$diff),rep(c("A","B","C","D"),3)), 
                        protocol=c(as.character(design$protocol),rep("Ramond",12)),
                        protocol_sort = c(as.character(design$protocol_sort),rep(c("A","B","C","D"),3)))

dds_merged<- DESeqDataSetFromMatrix(
  countData = counts_merged,
  colData = design_merged,
  design = ~1)

vsd_merged <- varianceStabilizingTransformation(dds_merged) 
vstMat_merged = assay(vsd_merged)

### rm batch effects:
library(sva)
batch<-c(rep("gfp",26),rep("ramond",12))

modcombat<-model.matrix(~1, data=design_merged)
cv=apply(vstMat_merged, 1, function(x) sd(x)/mean(x))
combat_mydata= ComBat(dat=vstMat_merged[which(cv>0),], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# Ramond A: yellow	Ramond B: goldenrod2 	Ramond C: indianred2	Ramond D: blue

col_sort_merged = c("yellow","coral","darkred","red","goldenrod2","cyan","darkblue","cornflowerblue","indianred2","darkolivegreen1","darkgreen","green",
                 "blue")
col_merged = sapply(design_merged$protocol_sort, function(x) col_sort_merged[which(levels(design_merged$protocol_sort) == x)])

pch_diff_merged=c( rep(7,4), 16,17,15)
pch_merged = sapply(design_merged$diff, function(x) pch_diff_merged[which(levels(design_merged$diff) == x)])

labels=c("A GFP+","A presort","A GFP-","B GFP+","B presort","B GFP-","C GFP+","C presort","C GFP-",
         "popA", "popB", "popC", "popD")

pdf("MDS.Ramond.pdf")
par(mar=c(5.1,4.1,4.1,10), xpd=T)
plotMDS(combat_mydata, col = col_merged, pch = pch_merged, cex=2)
legend("topright",inset=c(-0.35,0),labels,col=col_sort_merged[c(3,4,2,7,8,6,11,12,10,1,5,9,13)], pch=c(rep(16,9),rep(7,4)),bty="n", title="Protocol /\nCell population")
par(mar=c(5.1,5.1,4.1,1))
plotMDS(combat_mydata, col = col_merged, pch = pch_merged, cex=2, gene.selection="common", cex.axis=1.5, cex.lab=1.5)

dev.off()



