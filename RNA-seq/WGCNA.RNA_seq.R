### WGCNA on RNA-seq #####

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18/")
library("WGCNA")
allowWGCNAThreads()
library(dplyr)
library(annotables)
library(extrafont)
font_import()
y
load("endocomp.Rdata")


### only only protein coding genes expressed at >= 1TPM in >= 3 samples:
### filter the results to only protein-coding genes

use=rownames(tpms)[which(rowSums(tpms>=1)>=3)] # 30393
wgcna_vst=vstMat[which(rownames(vstMat) %in% use),]

prot_coding = grch37$ensgene[grch37$biotype == "protein_coding"]
ensrownames=sapply(strsplit(rownames(wgcna_vst),split="\\."), function(x) x[1])
coding_vst = wgcna_vst[which(ensrownames %in% prot_coding),]
coding_gene_ann = gene_ann[rownames(coding_vst),,drop=F]
wgcna_vst=coding_vst
dim(wgcna_vst)
# [1] 15555    26

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

##### run WGCNA module detection ### 

degData = t(wgcna_vst)
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

net = blockwiseModules(degData, power = 6, maxBlockSize = 20000,
                       TOMType = "signed", networkType="signed",
                       numericLabels = TRUE, pamStage = T, pamRespectsDendro = T,
                       verbose = 3, minModuleSize=50, deepSplit=1)
table(net$colors)

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
# 1757 2238 1583 1363 1221 1079 1020  712  585  531  484  467  436  417  378  326  238  225  179  178  138 

save(net, file="net.20M.Rda")

pdf("WGCNA_dendrogram.pdf", width = 12, height = 9)
moduleColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

lsdatME=moduleEigengenes(degData,moduleColors)$eigengenes
#datME=moduleEigengenes(degData,moduleColors)$eigengenes
### use numbers instead of colours:
datME=moduleEigengenes(degData,moduleLabels)$eigengenes

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

#MEs = moduleEigengenes(degData, moduleColors)$eigengenes
MEs = moduleEigengenes(degData, moduleLabels)$eigengenes

order=c(4:8,1:3,12:17,9:11,21:26,18:20)

pdf("Module_eigengenes.barplots.pdf",width=10, height=4, family="Arial")
par(las=2, mar=c(8,5,4,1))
for (i in 1:length(unique(moduleLabels))){
  which.module = unique(moduleLabels)[i]
  ME=datME[, paste("ME",which.module, sep="")]
  barplot(ME[order], col=col[order], ylab="eigengene expression",xlab="",main=which.module, names.arg=rownames(design)[order], cex.names=0.75)
}
dev.off()

 gene2module = data.frame(gene=colnames(degData), geneName=sapply(colnames(degData), function(x) as.character(gene_ann[x,])), module=paste0("M",moduleLabels))
 write.table(gene2module, file = "WGCNA.gene2module.20M.txt", sep = "\t", row.names=F, quote=F)

datME=moduleEigengenes(degData,moduleLabels)$eigengenes
rownames(datME)=colnames(vstMat)
write.table(datME[order,-1], file="RNA.module_eigengenes.endo_comp.txt",sep="\t",quote=F)

library(RColorBrewer)
pdf("datME.heatmap.pdf")
colnames(datME)=gsub("ME","R_M",colnames(datME))
heatmap(data.matrix(datME[rev(order),]), col=brewer.pal(9,"Reds"), Rowv=NA)
dev.off()

save.image("endocomp.WGCNA.Rdata")

load("endocomp.WGCNA.Rdata")

#######  enrichment analyses for WGCNA modules ########
gene2module[gene2module$geneName=="NKX6-1",]
# gene geneName module
# ENSG00000163623.5 ENSG00000163623.5   NKX6-1     M5
gene2module[gene2module$geneName=="PDX1",]
# gene geneName module
# ENSG00000139515.5 ENSG00000139515.5     PDX1     M5

#### which module is enriched for the PE signature? 
Cebola_PE_signature=scan("external_data/TEAD_YAP.500genes.txt",what="character")

### signatures of different cell populations from Ramond: 
popA=scan("external_data/Ramond2018.popA.txt",what="character")
popB=scan("external_data/Ramond2018.popB.txt",what="character")
popC=scan("external_data/Ramond2018.popC.txt",what="character")
popD=scan("external_data/Ramond2018.popD.txt",what="character")

### enrichment for known developmental genes from GO:
go_digestive=scan("external_data/GO_0048565.digestive_tract_development.txt",what="character")
go_intestinal=scan("external_data/GO_0060575.intestinal_epithelial_cell_differentiation.txt",what="character")
go_pancreas=scan("external_data/GO_0031016.pancreas_development.txt",what="character")
go_liver=scan("external_data/GO_0001889.liver_development.txt",what="character")

##### three lists of genes differentially and highly expressed in ISC derived from duodenum, jejunum and ileum 
##### Wang 2016, Nature, Fig2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4853906/figure/F8/
duodenum=scan("external_data/Wang2016.duodenum_genes.txt",what="character")
jejunum=scan("external_data/Wang2016.jejunum_genes.txt",what="character")
ileum=scan("external_data/Wang2016.ileum_genes.txt",what="character")

#### Jennings 2017 - DEGs - markers of dorsal pancreas, hepatic cords, hepatobiliary primordium
dp=scan("external_data/Jennings_2017/DP.txt",what="character")
hbp=scan("external_data/Jennings_2017/HBP.txt",what="character")
hc=scan("external_data/Jennings_2017/HC.txt",what="character")

gene_lists = list(PE_GFPpos=scan("PE_iPSC_signature.GFPpos.txt",what="character"),
                  PE_pre=scan("PE_iPSC_signature.presort.txt",what="character"),
                  PE_GFPneg=scan("PE_iPSC_signature.GFPneg.txt",what="character"),
                  Cebola_PE_signature=Cebola_PE_signature, popA=popA, popB=popB, popC=popC, popD=popD,
                  go_digestive=go_digestive, go_intestinal=go_intestinal, go_pancreas=go_pancreas, go_liver=go_liver,
                  duodenum=duodenum, jejunum=jejunum, ileum=ileum, dp=dp, hbp=hbp, hc=hc
                  
)
####################


hypergeo<-function(genes, gene2module){
  enrichment = rep(1, length(levels(gene2module$module)))
  for(i in 1:length(unique(gene2module$module))) {
    A=length(intersect(genes, gene2module$geneName[gene2module$module==levels(gene2module$module)[i]] )) ## overlap
    B=length(gene2module$geneName[gene2module$module==levels(gene2module$module)[i]])## module size
    C=length(gene2module$geneName)
    D=length(which(genes %in% gene2module$geneName))
    enrichment[i]= phyper(A-1,B,C-B,D,lower.tail=F)
  }
  
  names(enrichment) = levels(gene2module$module)
  return(enrichment)
}

enrichments = data.frame(lapply(gene_lists, function(x) hypergeo(x, gene2module)))
write.table(enrichments, file="WGCNA.module_enrichments.txt",sep="\t",quote=F)
