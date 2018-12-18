library(ggplot2)
library(reshape)
library("gridExtra")

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18//")
enr=read.table("RNA_modules.enrichments.plot_input.txt",h=T,sep="\t")
enr_GO = enr[,c(1:8)]
enr_pops= enr[c(1,9:12)]

### separate the GO enrichments and pops enrichment plots:
enr_GO_melt=melt(enr_GO)
enr_GO_melt$Module=factor(enr_GO_melt$Module, levels=rev(c("R_M1","R_M2","R_M3","R_M5","R_M7","R_M8","R_M11","R_M13","R_M16","R_M20")))
levels(enr_GO_melt$variable) = c("PP signature\nGFP+","PP signature\n(Cebola 2015)","GO:0031016\npancreas development",
                                 "GO:0048565 digestive\ntract development",
                              "Duodenum signature\n(Wang 2015)", "Ileum signature\n(Wang 2015)", "Hepatic Cords\n(Jennings 2017)")
enr_GO_melt$variable = factor(enr_GO_melt$variable,rev(levels(enr_GO_melt$variable)))

enr_pops_melt=melt(enr_pops)
enr_pops_melt$Module=factor(enr_pops_melt$Module, levels=rev(c("R_M1","R_M2","R_M3","R_M5","R_M7","R_M8","R_M11","R_M13","R_M16","R_M20")))
levels(enr_pops_melt$variable) =  c(
                                    "Fetal panc.\nprogenitors 1", "Fetal panc.\nprogenitors 2", 
                                    "Fetal endocrine\nprogenitors", "Fetal endocrine\ncells")
enr_pops_melt$variable = factor(enr_pops_melt$variable,rev(levels(enr_pops_melt$variable)))

pdf("Fig4.GExp_enrichments.v3.pdf",height=10,width=5, family="Arial")
p1=ggplot(enr_GO_melt, aes(fill=Module, y=-log10(value), x=variable)) +
  geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
  scale_fill_brewer(palette = "Paired") + theme_bw(base_size = 14) + coord_flip() +ylab("-log10(p-value)") + xlab("") + guides(fill = guide_legend(reverse=TRUE)) 
p2=ggplot(enr_pops_melt, aes(fill=Module, y=-log10(value), x=variable)) +
  geom_bar( stat="identity", position=position_dodge(width=0.5)) + coord_cartesian(ylim=c(0,50)) +
  scale_fill_brewer(palette = "Paired") + theme_bw(base_size = 14) + coord_flip() +ylab("-log10(p-value)") + xlab("") + guides(fill = guide_legend(reverse=TRUE)) 

grid.arrange(p1, p2, nrow = 2, layout_matrix=matrix(c(1,1,1,2), nrow=4))
dev.off()

setwd("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18/")
enr_ATAC=read.table("ATAC_enrichments.txt",h=T)
  enr_ATAC_melt=melt(enr_ATAC)
  enr_ATAC_melt$Module=factor(enr_ATAC_melt$Module, levels=rev(c("A_M1","A_M2","A_M3","A_M4","A_M5","A_M6")))
  enr_ATAC_melt$variable = factor(enr_ATAC_melt$variable,rev(levels(enr_ATAC_melt$variable)))
  
pdf("Fig4E.ATAC_enrichments.v2.pdf", height=10,width=4)
  ggplot(enr_ATAC_melt, aes(fill=Module, y=-log10(value), x=variable)) +
    geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
    scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 14) + coord_flip() +ylab("-log10(p-value)") + 
    xlab("") + guides(fill = guide_legend(reverse=TRUE)) 
  dev.off()
  
  