### expression plots ###### 
setwd("~/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18/")
load("endocomp.Rdata")

### plot NKX6.1 expression
library("beeswarm")

design$protocol_sort = factor(design$protocol_sort, levels=c("A_presort","A_GFPpos","A_GFPneg",
                                                             "B_presort","B_GFPpos","B_GFPneg",
                                                             "C_presort","C_GFPpos","C_GFPneg"))
design$ABC = factor(c(rep("A: GFP-",3),rep("A: GFP+",3),rep("A: presort",2), rep("B: GFP-",3), rep("B: GFP+",3), rep("B: presort",3),rep("C: GFP-",3),rep("C: GFP+",3),rep("C: presort",3)), levels=c("A: GFP+","A: presort","A: GFP-","NA1","B: GFP+","B: presort","B: GFP-","NA2","C: GFP+","C: presort","C: GFP-"))


gene="CDX2"
pdf(paste0(gene,".pdf"),width=7, height=4, family="Arial")
par(las=2)
ensgene=rownames(gene_ann[which(gene_ann$GeneName == gene),,drop=F])
ensgene=gsub("\\..","",ensgene)
#beeswarm(as.numeric(tpms[ensgene,]) ~ design$ABC, pwpch=pch, pwcol=col,cex=2.5, main=gene,ylab= paste(gene,"TPMs"), xlab="", cex.main=2)

beeswarm(as.numeric(tpms[ensgene,]) ~ design$ABC, pwpch=pch, pwcol=col,cex=1.5, main=gene,
         ylab= paste(gene,"TPMs"), xlab="", xaxt="n", cex.main=1.5, bty="n", 
         ylim=c(0,max(as.numeric(tpms[ensgene,]))+0.1*max(as.numeric(tpms[ensgene,]))), yaxs="i")

for (i in c(1:3,5:7,9:11)){
  y=mean(as.numeric(tpms[ensgene,which(design$ABC==levels(design$ABC)[i])]))
  segments(i-0.4, y, i+0.4,y,lty=1, lwd=2 ) 
}

axis(1, labels = FALSE, tick=T, at=c(0,12), lty=1)
labels = c("A: GFP+","A: presort","A: GFP-","","B: GFP+","B: presort","B: GFP-","","C: GFP+","C: presort","C: GFP-")
text(1:11, par("usr")[3] - 0.25, srt = 45, adj = 1.2,
     labels = labels, xpd = TRUE)
#abline(v=c(3.5, 6.5), lty=2)
dev.off()

#### for NKX6-1 we want to label it NKX6.1 so it's consistent with the rest of the paper
gene="NKX6-1"
pdf(paste0(gene,".pdf"),width=7, height=4, family="Arial")
par(las=2)
ensgene=rownames(gene_ann[which(gene_ann$GeneName == gene),,drop=F])
ensgene=gsub("\\..","",ensgene)
#beeswarm(as.numeric(tpms[ensgene,]) ~ design$ABC, pwpch=pch, pwcol=col,cex=2.5, main=gene,ylab= paste(gene,"TPMs"), xlab="", cex.main=2)

beeswarm(as.numeric(tpms[ensgene,]) ~ design$ABC, pwpch=pch, pwcol=col,cex=1.5, main="NKX6.1",
         ylab= paste("NKX6.1 TPMs"), xlab="", xaxt="n", cex.main=1.5, bty="n", 
         ylim=c(0,max(as.numeric(tpms[ensgene,]))+0.1*max(as.numeric(tpms[ensgene,]))), yaxs="i")

for (i in c(1:3,5:7,9:11)){
  y=mean(as.numeric(tpms[ensgene,which(design$ABC==levels(design$ABC)[i])]))
  segments(i-0.4, y, i+0.4,y,lty=1, lwd=2 ) 
}

axis(1, labels = FALSE, tick=T, at=c(0,12), lty=1)
labels = c("A: GFP+","A: presort","A: GFP-","","B: GFP+","B: presort","B: GFP-","","C: GFP+","C: presort","C: GFP-")
text(1:11, par("usr")[3] - 0.25, srt = 45, adj = 1.2,
     labels = labels, xpd = TRUE)
#abline(v=c(3.5, 6.5), lty=2)
dev.off()



####### gene expression plots for MODY and neonatal diabetes genes #############

neonatal=scan("external_data/neonatal_diabetes_genes.txt",what="character")
MODY=scan("external_data/MODY_genes.txt",what="character")
genes = sort(unique(c(neonatal,MODY)))
genes=gsub("\\.","-",genes)

pdf("MODY_and_neonatal.gene_exp.pdf", paper="a4", family="Arial", width=8.3, height=11.7)
par(mfrow=c(6,2), las=2, mar=c(5.1,4,3,1))

for (gene in genes){
  ensgene=rownames(gene_ann[which(gene_ann$GeneName == gene),,drop=F])
  ensgene=gsub("\\..+","",ensgene)

  beeswarm(as.numeric(tpms[ensgene,]) ~ design$ABC, pwpch=pch, pwcol=col,cex=1.5, main=gsub("-","\\.",gene),
           ylab= paste0(gsub("-","\\.",gene), " TPMs"), xlab="", xaxt="n", cex.main=1.5, bty="n", 
           ylim=c(0,1.1*max(as.numeric(tpms[ensgene,]))), yaxs="i")
  
  for (i in c(1:3,5:7,9:11)){
    y=mean(as.numeric(tpms[ensgene,which(design$ABC==levels(design$ABC)[i])]))
    segments(i-0.4, y, i+0.4,y,lty=1, lwd=2 ) 
  }
  
  axis(1, labels = FALSE, tick=T, at=c(0,12), lty=1)
  labels = c("A: GFP+","A: presort","A: GFP-","","B: GFP+","B: presort","B: GFP-","","C: GFP+","C: presort","C: GFP-")
  text(1:11, par("usr")[3] - 0.25, srt = 45, adj = 1.2,
       labels = labels, xpd = TRUE)
}
dev.off()




