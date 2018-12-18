### ATAC analysis - global similarities #####

## load libraries:
library(edgeR)
library("DESeq2")
library(sva)

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18/")
counts = read.table("25.07.18.merged_peaks.clean.txt",h=T,row.names=1)
peak_ann = counts[,1:5,drop=F]
counts = counts[,-c(1:5)]

design=read.table("design_matrix",h=T)

## rename samples
colnames(counts) = as.character(sapply(colnames(counts), function(x) design$sample[which(design$seq_id == x)]))

### MDS plot
rownames(design) = design$sample
design = design[colnames(counts),]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = design,
  design = ~1)

dds=estimateSizeFactors(dds)
norm_counts = counts(dds, norm=T)

vsd <- varianceStabilizingTransformation(dds) 
vstMat = assay(vsd)

### col by method: Beta - red, Hebrok - blue, Nostro - green
### col intensity - sort: darkest - GFP+, lightest - GFP-, med - presort

# col by sort (captures both subject and line), different reps will be the same col
# GFP-  red, GFP+  blue, presort - green
# coral, red, darkred
# cyan, cornflowerblue, darkblue
# chartreuse, green, darkgreen

design$protocol_sort = factor(paste(design$protocol, design$GFP, sep="_"))

col_sort = c("coral","darkred","red","cyan","darkblue","cornflowerblue","darkolivegreen1","darkgreen","green")
col = sapply(design$protocol_sort, function(x) col_sort[which(levels(design$protocol_sort) == x)])

design$diff = factor(design$diff)
pch_diff=c(16,17,15)
pch = sapply(design$diff, function(x) pch_diff[which(levels(design$diff) == x)])

labels=gsub("_"," ",levels(design$protocol_sort))

### batch correction for pooling:
batch<-design$pool
modcombat<-model.matrix(~1, data=design)
combat_mydata= ComBat(dat=vstMat, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

labels=c("A: GFP+","A: presort","A: GFP-",  "B: GFP+","B: presort","B: GFP-", "C: GFP+","C: presort","C: GFP-")

save.image("endo_comp.atac.Rdata")

pdf("Fig2C.ATAC_PCA.pdf")
par(mar=c(5.1,5.1,1,10), xpd=T)
plotMDS(combat_mydata, col = col, pch = pch, cex=2, gene.selection="common", main="", cex.lab=1.5, cex.axis=1.5)
legend(1.1,0.6,labels,col=col_sort[c(2,3,1,5,6,4,8,9,7)], pch=16,bty="n", title="Protocol", cex=1.5)
legend(1.1,0.05,c("d1","d2","d3"),pch=pch_diff,bty="n",title="Diff", cex=1.5)

par(mar=c(5.1,5.1,1,2))
plotMDS(combat_mydata, col = col, pch = pch, cex=2, gene.selection="common", main="", cex.lab=1.5, cex.axis=1.5)
dev.off()

