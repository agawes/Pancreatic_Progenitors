#### Analysis 2 - Figure 3 #####
#### PE signature #####

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

sort=sapply(strsplit(colnames(counts), split="_"), function(x) x[2])
diff = sapply(strsplit(colnames(counts), split="_"), function(x) x[3])
protocol=sapply(strsplit(colnames(counts), split="_"), function(x) x[1])

design = data.frame(row.names = colnames(counts), sort=sort, diff = diff, protocol=protocol)
design$protocol_sort = factor(paste(design$protocol, design$sort, sep="_"))

####### PE signature - across the methods #########
# compare to FPKMs from MSB-12-862-s001_human_tissue_FPKMs (on Desktop) - from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4848759/#msb155865-sup-0001

gtex_fpkms=read.table("external_data/GTEx_tissues_fpkm.txt",h=T,row.names=1,sep="\t")
gtex_ann = gtex_fpkms[,1,drop=F]
gtex_fpkms = gtex_fpkms[,-c(1,32,33)]

### first need to convert the FPKM's to TPM's:
gtex_tpms=matrix(, ncol=ncol(gtex_fpkms), nrow=nrow(gtex_fpkms))
colnames(gtex_tpms) = gsub("fpkm","tpm",colnames(gtex_fpkms))
rownames(gtex_tpms) = rownames(gtex_fpkms)
for (i in 1:ncol(gtex_fpkms)){
  gtex_tpms[,i] = sapply(gtex_fpkms[,i], function(x) x/sum(gtex_fpkms[,i])*10^6)
}

# merge with TPM table:
rownames(tpms) = sapply(strsplit(rownames(tpms),split="\\."),function(x) x[1])
merged_TPMs = merge(tpms, gtex_tpms, by="row.names")
rownames(merged_TPMs) = merged_TPMs$Row.names
merged_TPMs=merged_TPMs[,-1]

#merged_TPMs=merged_TPMs[,-9] ## rm outlier
# rescale the TPMs for endocomp to total 1e06:
for (i in 1:26){
  merged_TPMs[,i] = sapply(merged_TPMs[,i], function(x) x/sum(merged_TPMs[,i])*10^6)
}
merged_TPMs$gene = sapply(rownames(merged_TPMs), function(x) gtex_ann[x,])


#### calculate PE signatures for each method & fraction; report overlap
#### mean TPMs per method:
merged_TPMs$A_presort = apply(merged_TPMs[,grep("A_presort",colnames(merged_TPMs))], 1, mean)
merged_TPMs$A_GFP_plus = apply(merged_TPMs[,grep("A_GFPpos",colnames(merged_TPMs))], 1, mean)
merged_TPMs$A_GFP_minus = apply(merged_TPMs[,grep("A_GFPneg",colnames(merged_TPMs))], 1, mean)

merged_TPMs$B_presort = apply(merged_TPMs[,grep("B_presort",colnames(merged_TPMs))], 1, mean)
merged_TPMs$B_GFP_plus = apply(merged_TPMs[,grep("B_GFPpos",colnames(merged_TPMs))], 1, mean)
merged_TPMs$B_GFP_minus = apply(merged_TPMs[,grep("B_GFPneg",colnames(merged_TPMs))], 1, mean)

merged_TPMs$C_presort = apply(merged_TPMs[,grep("C_presort",colnames(merged_TPMs))], 1, mean)
merged_TPMs$C_GFP_plus = apply(merged_TPMs[,grep("C_GFPpos",colnames(merged_TPMs))], 1, mean)
merged_TPMs$C_GFP_minus = apply(merged_TPMs[,grep("C_GFPneg",colnames(merged_TPMs))], 1, mean)

merged_TPMs = merged_TPMs[,c(27:66)]
### calculate CV for presort, plus, minus

merged_TPMs$CV =apply(merged_TPMs[,c(grep("GTEx",colnames(merged_TPMs)), grep("A_|B_|C_",colnames(merged_TPMs)))], 1, function(x) sd(x)/mean(x))
#merged_TPMs$CV_presort = apply(merged_TPMs[,grep("presort$",colnames(merged_TPMs))], 1, function(x) sd(x)/mean(x))
#merged_TPMs$CV_plus = apply(merged_TPMs[,grep("plus",colnames(merged_TPMs))], 1, function(x) sd(x)/mean(x))
#merged_TPMs$CV_minus = apply(merged_TPMs[,grep("minus",colnames(merged_TPMs))], 1, function(x) sd(x)/mean(x))

### for being a PE signature gene, require:
### merged_TPMs$CV>1 --> tissue selectivity
### NOT: merged_TPMs$CV_presort <= 1 --> stable expression between methods
### merged_TPMs$Z_A_presort >1 &  merged_TPMs$Z_B_presort >1 & merged_TPMs$Z_C_presort >1
### merged_TPMs$A_presort >= 1 & merged_TPMs$B_presort >= 1 & merged_TPMs$C_presort >=1

merged_TPMs$Z_A_presort = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$A_presort-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("A_presort",colnames(merged_TPMs))) ])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("A_presort",colnames(merged_TPMs))) ])))
merged_TPMs$Z_B_presort = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$B_presort-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("B_presort",colnames(merged_TPMs))) ])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("B_presort",colnames(merged_TPMs))) ])))
merged_TPMs$Z_C_presort = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$C_presort-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("C_presort",colnames(merged_TPMs))) ])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("C_presort",colnames(merged_TPMs))) ])))

PE_presort_signature.A = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_A_presort >1 & merged_TPMs$A_presort >= 1, select=gene)
PE_presort_signature.B = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_B_presort >1 & merged_TPMs$B_presort >= 1 , select=gene)
PE_presort_signature.C = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_C_presort >1 & merged_TPMs$C_presort >=1, select=gene)
PE_presort_signature=subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$CV_presort <= 1 & merged_TPMs$Z_A_presort >1 
                            & merged_TPMs$Z_B_presort >1 & merged_TPMs$Z_C_presort >1 & merged_TPMs$A_presort >= 1 & merged_TPMs$B_presort >= 1 
                            & merged_TPMs$C_presort >=1, select=gene)

######## Venn Diagram of the PE signatures:
PE_presort_signatures <- list(A=as.character(PE_presort_signature.A$gene), B=as.character(PE_presort_signature.B$gene), 
                              C=as.character(PE_presort_signature.C$gene))
require(VennDiagram)
pdf("venn.PE_signature.presort.pdf")
venn.plot <- venn.diagram(PE_presort_signatures , NULL, fill=c("red", "blue","green"), alpha=c(0.5,0.5,0.5), cex = 3, cat.fontface=4, category.names=c("Protocol A", "Protocol B", "Protocol C"), main="PE signature - presort",  scaled=T)
grid.draw(venn.plot)
dev.off()

###  it doesn't work --> I used the http://eulerr.co/ app inputing the numbers

PE_presort_signature=intersect(intersect(PE_presort_signatures$A, PE_presort_signatures$B),PE_presort_signatures$C)

select=unlist(sapply(PE_presort_signature, function(x) which(merged_TPMs$gene ==x)[1]))
PE_signature_SupTable = cbind(PE_presort_signature, merged_TPMs[select,c(31:44)])
write.table(PE_signature_SupTable, file="STable.PE_presort_signature.txt",sep="\t",quote=F)


## for GFP+ and GFP- specific genes, keep both GFP+ and GFP-
merged_TPMs$Z_A_plus = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$A_GFP_plus-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("A_GFP_plus",colnames(merged_TPMs)))])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("A_GFP_plus",colnames(merged_TPMs)))])))
merged_TPMs$Z_B_plus = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$B_GFP_plus-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("B_GFP_plus",colnames(merged_TPMs)))])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("B_GFP_plus",colnames(merged_TPMs)))])))
merged_TPMs$Z_C_plus = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$C_GFP_plus-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("C_GFP_plus",colnames(merged_TPMs)))])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("C_GFP_plus",colnames(merged_TPMs)))])))

PE_GFP_plus_signature.A = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_A_plus >1 & merged_TPMs$A_GFP_plus >= 1, select=gene)
PE_GFP_plus_signature.B = subset(merged_TPMs, merged_TPMs$CV>1 &  merged_TPMs$Z_B_plus >1 & merged_TPMs$B_GFP_plus >= 1 , select=gene)
PE_GFP_plus_signature.C = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_C_plus >1 & merged_TPMs$C_GFP_plus >=1, select=gene)

######## Venn Diagram of the PE signatures:
PE_GFP_plus_signatures <- list(A=as.character(PE_GFP_plus_signature.A$gene), B=as.character(PE_GFP_plus_signature.B$gene), 
                               C=as.character(PE_GFP_plus_signature.C$gene))
pdf("venn.PE_signature.GFP_plus.pdf")
venn.plot <- venn.diagram(PE_GFP_plus_signatures , NULL, fill=c("red", "blue","green"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("Protocol A", "Protocol B", "Protocol C"), main="PE signature - GFP+",  scaled=T)
grid.draw(venn.plot)
dev.off()

PE_GFP_plus_signature=intersect(intersect(PE_GFP_plus_signatures$A, PE_GFP_plus_signatures$B),PE_GFP_plus_signatures$C)

select=unlist(sapply(PE_GFP_plus_signature, function(x) which(merged_TPMs$gene ==x)[1]))
PE_signature_plus_SupTable = cbind(PE_GFP_plus_signature, merged_TPMs[select,c(31:41,45:47)])
write.table(PE_signature_plus_SupTable, file="STable.PE_presort_signature.GFPpos.txt",sep="\t",quote=F)

merged_TPMs$Z_A_minus = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$A_GFP_minus-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("A_GFP_minus",colnames(merged_TPMs)))])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("A_GFP_minus",colnames(merged_TPMs)))])))
merged_TPMs$Z_B_minus = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$B_GFP_minus-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("B_GFP_minus",colnames(merged_TPMs)))])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("B_GFP_minus",colnames(merged_TPMs)))])))
merged_TPMs$Z_C_minus = sapply(rownames(merged_TPMs), function(x) 
  (merged_TPMs[x,]$C_GFP_minus-mean(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("C_GFP_minus",colnames(merged_TPMs)))])))
  /sd(as.numeric(merged_TPMs[x,c(grep("GTEx",colnames(merged_TPMs)), grep("C_GFP_minus",colnames(merged_TPMs)))])))


PE_GFP_minus_signature.A = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_A_minus >1 & merged_TPMs$A_GFP_minus >= 1, select=gene)
PE_GFP_minus_signature.B = subset(merged_TPMs, merged_TPMs$CV>1 &  merged_TPMs$Z_B_minus >1 & merged_TPMs$B_GFP_minus >= 1 , select=gene)
PE_GFP_minus_signature.C = subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$Z_C_minus >1 & merged_TPMs$C_GFP_minus >=1, select=gene)
PE_GFP_minus_signature=subset(merged_TPMs, merged_TPMs$CV>1 & merged_TPMs$CV_minus <= 1 & merged_TPMs$Z_A_minus >1 &  merged_TPMs$Z_B_minus >1 
                              & merged_TPMs$Z_C_minus >1 & merged_TPMs$A_GFP_minus >= 1 & merged_TPMs$B_GFP_minus >= 1 & merged_TPMs$C_GFP_minus >=1, select=gene)

######## Venn Diagram of the PE signatures:
PE_GFP_minus_signatures <- list(A=as.character(PE_GFP_minus_signature.A$gene), B=as.character(PE_GFP_minus_signature.B$gene), 
                                C=as.character(PE_GFP_minus_signature.C$gene))
pdf("venn.PE_signature.GFP_minus.pdf")
venn.plot <- venn.diagram(PE_GFP_minus_signatures , NULL, fill=c("red", "blue","green"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("Protocol A", "Protocol B", "Protocol C"), main="PE signature - GFP-",  scaled=T)
grid.draw(venn.plot)
dev.off()

PE_GFP_minus_signature=intersect(intersect(PE_GFP_minus_signatures$A, PE_GFP_minus_signatures$B),PE_GFP_minus_signatures$C)
select=unlist(sapply(PE_GFP_minus_signature, function(x) which(merged_TPMs$gene ==x)[1]))
PE_signature_minus_SupTable = cbind(PE_GFP_minus_signature, merged_TPMs[select,c(31:41,48:50 )])
write.table(PE_signature_minus_SupTable, file="STable.PE_signature.GFPneg.txt",sep="\t",quote=F)

save(merged_TPMs, file="GTEx_TPMs.PE_Z_scores.Rdata")
