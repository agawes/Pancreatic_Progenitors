### RNA and ATAC WGCNA module eigengenes #### 

### all labels in Arial! - use library(extrafont)
library(extrafont)
font_import()
y

RNA_ME=read.table("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18/RNA.module_eigengenes.endo_comp.txt",
                  sep="\t",h=T,row.names=1)

ATAC_ME=read.table("/Users/agata/Desktop/StemBANCC/endodermal competence/FINAL_ATAC_Reanalysis_October18/ATAC.module_eigengenes.endo_comp.txt",
                   sep="\t",h=T,row.names=1)
## rm the sample missing from RNA-seq
ATAC_ME=ATAC_ME[-6,]
## rm the M0 module
ATAC_ME=ATAC_ME[,-1]
rownames(RNA_ME)=rownames(ATAC_ME)

write.table(cor(ATAC_ME, RNA_ME), file="RNA_vs_ATAC.module_eigengene_cor.txt",sep="\t",quote=F)

cormat=round(cor(ATAC_ME, RNA_ME),2)
colnames(cormat)=gsub("ME","R_M",colnames(cormat))
rownames(cormat)=gsub("ME","A_M",rownames(cormat))

#cluster cor's
hc=heatmap(cormat)
cormat_hc=cormat[hc$rowInd,hc$colInd]

library(reshape2)
library(ggplot2)
melted_cormat <- melt(cormat)

melted_cormat_hc <- melt(cormat_hc)

# Heatmap
pdf("ATAC_RNA_cor_matrix.pdf",width=8,height=4, family="Arial")
# ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   theme_minimal()+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1))+
#   coord_fixed() +xlab("RNA-seq") + ylab("ATAC-seq")

ggplot(data = melted_cormat_hc, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(size=12))+
  coord_fixed() +xlab("RNA-seq") + ylab("ATAC-seq")

dev.off()
