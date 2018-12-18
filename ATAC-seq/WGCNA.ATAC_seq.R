library(edgeR)
library("DESeq2")
library(sva)

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18/")
load("endo_comp.atac.Rdata")


#### run WGCNA
library("WGCNA")
allowWGCNAThreads()

### only peaks with at least 10x in >= 3 samples & with CV > 0.5
peak_CV = apply(norm_counts, 1, function(x) sd(x)/mean(x))
length(intersect(which(peak_CV>0.5),which(rowSums(norm_counts>=10)>=3)))
## [1] 19739

subset=intersect(which(peak_CV>0.5),which(rowSums(norm_counts>=10)>=3))
counts_subset=counts[subset,]

use=rownames(counts[subset,])
wgcna_vst=vstMat[which(rownames(vstMat) %in% use),]
dim(wgcna_vst)
# [1] 19739    27

peak_ann$peak=rownames(peak_ann)
write.table(peak_ann[use,c(1:3,6)], file="WGCNA.input_peaks.bed",sep="\t",quote=F, row.names = F, col.names = F)
peak_gene_ann = read.table("20181031-public-3.0.0-gRgKAm-hg19-all-region.txt",skip=1, sep="\t")
degData = t(wgcna_vst)
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

### pick power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
pdf(file = "WGCNA_softThreshold.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


##### run WGCNA module detection

net = blockwiseModules(degData, power = 3, maxBlockSize = 21000,
                       TOMType = "signed", networkType="signed",
                       numericLabels = TRUE, pamStage = T, pamRespectsDendro = T,
                       saveTOMs=T, verbose = 3, deepSplit=4, minModuleSize=50)
table(net$colors)

# 0    1    2    3    4    5    6    7    8    9   10   11 
# 5231 7598 1773 1742  957  818  588  466  216  198   83   69 

moduleLabels = net$colors
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

lsdatME=moduleEigengenes(degData,moduleLabels)$eigengenes
datME=moduleEigengenes(degData,moduleLabels)$eigengenes
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

MEs = moduleEigengenes(degData, moduleLabels)$eigengenes

order=c(1,6,24,18,21,13,17,7, 12, 2, 22, 14, 19, 9, 26, 3, 8, 25, 20, 10, 15, 5, 11, 27, 4, 23, 16 )

pdf("Module_eigengenes.barplots.pdf",width=10, height=4)
par(las=2, mar=c(8,5,4,1))
for (i in 1:length(unique(moduleLabels))){
  which.module = unique(moduleLabels)[i]
  ME=datME[, paste("ME",which.module, sep="")]
  barplot(ME[order], col=col[order], ylab="eigengene expression",xlab="", main=which.module, names.arg=design$sample[order], cex.names=0.75)
}
dev.off()

peak2module = data.frame(peak=colnames(degData), module=paste0("M",moduleLabels), peak_genes=peak_gene_ann$V2)
peak2module_ann=merge(peak2module, peak_ann, by="peak")

write.table(peak2module_ann, file = "WGCNA.peak2module.txt", sep = "\t", row.names=F, quote=F)

rownames(datME) = colnames(wgcna_vst)
write.table(datME[order,], file="ATAC.module_eigengenes.endo_comp.txt",sep="\t",quote=F)

library(RColorBrewer)
pdf("datME.heatmap.pdf")
colnames(datME)=gsub("ME","A_M",colnames(datME))
heatmap(data.matrix(datME[rev(order),]), col=brewer.pal(9,"Reds"), Rowv=NA)
dev.off()

save.image("atac.wgcna.Rdata")
