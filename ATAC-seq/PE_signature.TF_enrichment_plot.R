## TF enrichment plots - PE signature

library(ggplot2)
library(reshape)
library("gridExtra")

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18//")
enr_ATAC=read.table("PE_signature.TFs_enriched.txt",h=T)
enr_ATAC_melt=melt(enr_ATAC)
names(enr_ATAC_melt)[2]="Population"
levels(enr_ATAC_melt$Population)= c("GFP-","GFP+","Presort")
enr_ATAC_melt$Population=factor(enr_ATAC_melt$Population, levels=c("GFP-","Presort","GFP+"))
  enr_ATAC_melt$tf_motif = factor(enr_ATAC_melt$tf_motif, levels=c("Cux2","Hnf6","Hnf1b","Fosl2","Foxa2","Hoxb4","Pdx1",
                         "Tead4","Gata1","Sox9")[10:1])

pdf("PE_signature.ATAC_enrichments.pdf",family="Arial")
ggplot(enr_ATAC_melt, aes(fill=Population, y=-log10(value), x=tf_motif)) +
  geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
  scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 14) + coord_flip() +ylab("-log10(p-value)") + 
  xlab("") + guides(fill = guide_legend(reverse=TRUE)) 
dev.off()

