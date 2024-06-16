# 制作三阴性乳腺癌和普通乳腺癌的分类原始数据
suppressPackageStartupMessages(library(ggplot2))#画图
suppressPackageStartupMessages(library(ggpubr))#画图辅助
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Oxidative_stress_new"
setwd(workingdir)
memory.limit(400000)
load(file = "paper.module.genes.and.pca.RData")
load(file ="tcga33_gtex18_samples11348v5620_and_genes56189.RData")
TNBC.samples = read.csv(file = 'TNBC.samples.csv',header = FALSE)$V1

all.tpm.final = all.tpm.final[,all.colinfo$Study == "BRCA"]
all.colinfo = all.colinfo[all.colinfo$Study == "BRCA",]
all.tpm.final = all.tpm.final[,all.colinfo$Sample_type != "Solid Tissue Normal"]
all.colinfo = all.colinfo[all.colinfo$Sample_type != "Solid Tissue Normal",]

genemap =read.csv("C:/Users/Dell/Documents/R/Metabolic Network/genes_recon1_map.csv")
genenames = unique(genemap$name)
# 得到表达矩阵和label
sample.expr = t(all.tpm.final[intersect(genenames,rownames(sample.expr)),])
sample.label = ifelse(substr(all.colinfo$Sample_ID,1,12) %in% TNBC.samples,"TNBC","non-TNBC")
write.csv(sample.expr,file = "RECON1_TNBC_expr.csv")
write.csv(sample.label,file = "RECON1_TNBC_label.csv")