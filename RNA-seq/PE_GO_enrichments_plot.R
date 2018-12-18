## Fig3C. PE signatures GO enrichment plot

library(ggplot2)
library(reshape)
library("gridExtra")

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18//")
enr_RNA=read.table("GO_BP.panc_signature.txt",h=T,sep="\t")
enr_RNA_melt=melt(enr_RNA)
names(enr_RNA_melt)[2]="Population"
levels(enr_RNA_melt$Population)= c("GFP+","Presort","GFP-")
enr_RNA_melt$Population=factor(enr_RNA_melt$Population, levels=c("GFP-","Presort","GFP+"))
levs=levels(enr_RNA_melt$GO_BP)[c(8,3,10, 6,5,4,7,9,1,2)]
enr_RNA_melt$GO_BP = factor(enr_RNA_melt$GO_BP, levels=levs)

pdf("PE_signature.RNA_enrichments.pdf", family="Arial")
ggplot(enr_RNA_melt, aes(fill=Population, y=-log10(value), x=GO_BP)) +
  geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
  scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 14) + coord_flip() +ylab("-log10(p-value)") + 
  xlab("") + guides(fill = guide_legend(reverse=TRUE)) 
dev.off()

