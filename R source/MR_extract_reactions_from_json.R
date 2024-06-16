# 这个文档用来读入json源文件
# 输入： Recon3D.json
# 输出文件：reactions_recon3D.RData
          #reactions_recon3D.csv
          # metabolites_recon3D.RData
          # metabolites_recon3D.csv
          # genes_recon3D.RData
          # genes_recon3D.csv

library("rjson")
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Metabolic Network"
setwd(workingdir)
memory.limit(400000)

# 读入原始文档
result <- fromJSON(file = "Recon3D.json")
################################################################################
# reactions信息读取
reactions = result$reactions  # notes与annotation是不同平台的命名对应，没什么用
reactions_id = sapply(reactions,function(x) x$id) # 化学反应的唯一标识：例如EX_5adtststerone_e
reactions_name = sapply(reactions,function(x) x$name) # 化学反应的描述性名称
reactions_lower = sapply(reactions,function(x) x$lower_bound) #一般用不到，只有求FBA时会使用
reactions_upper = sapply(reactions,function(x) x$upper_bound) # 一般用不到，只有求FBA时会使用
reactions_gene = sapply(reactions,function(x) x$gene_reaction_rule) # 包含了各个反应所需的所有酶
reactions_sub = sapply(reactions,function(x) x$subsystem) #正常代谢中的子系统标注
reactions_metabolites = sapply(reactions,function(x) x$metabolites)# 在计量矩阵中，正数表示该化学反应生成这种物质，负数表示该化学反应消耗这种物质
  # 求各个反应的底物和产物
products = c()
substrates = c()
for(i in c(1:length(reactions_metabolites)))
{
  metaForThisReaction = reactions_metabolites[[i]]
  coefs = sapply(metaForThisReaction, function(x) x[[1]]) # 取出系数
  products = c(products,paste(names(metaForThisReaction)[coefs>0],collapse = ",")) #根据系数取出产物
  substrates = c(substrates,paste(names(metaForThisReaction)[coefs<0],collapse = ","))#根据系数取出底物
}
# reactions_biggID = c()
# for(i in c(1:length(reactions)))
# {
#   tmp = reactions[[i]]$annotation[[1]]
#   reactions_biggID = append(reactions_biggID,tmp)
# }

reactions_recon3 = data.frame(id = reactions_id,
                              name = reactions_name,
                              lower_bound = reactions_lower,
                              upper_bound = reactions_upper,
                              gene = reactions_gene,
                              subsystem = reactions_sub,
                              substrates = substrates,
                              products = products)

save(reactions_recon3,file = "reactions_recon3D.RData")
write.csv(reactions_recon3,file = "reactions_recon3D.csv")
# 至此得到了reactions当中所有的有效信息
################################################################################

metabolites = result$metabolites
metabolites_id = sapply(metabolites,function(x) x$id) #5835
metabolites_name = sapply(metabolites,function(x) x$name) #5835
metabolites_compartment = sapply(metabolites,function(x) x$compartment) #5835
metabolites_charge = sapply(metabolites,function(x) x$charge) # 5835
metabolites_formula = sapply(metabolites,function(x) x$formula) # 5835
metabolites_formula[sapply(metabolites_formula,is.null)] <- NA
metabolites_formula = unlist(metabolites_formula)
metabolites_notes = c()
for(i in c(1:length(metabolites)))
{
  tmp = metabolites[[i]]$notes[[1]]
  metabolites_notes = append(metabolites_notes,tmp)
} #5835
metabolites_biocyc = c()
for(i in c(1:length(metabolites)))
{
  tmp = metabolites[[i]]$annotation$biocyc
  if(is.null(tmp)){
    metabolites_biocyc = append(metabolites_biocyc,NA)
  }else{metabolites_biocyc = append(metabolites_biocyc,paste(tmp,collapse = ","))}
} #5835
metabolites_hmdb = c()
for(i in c(1:length(metabolites)))
{
  tmp = metabolites[[i]]$annotation$hmdb
  if(is.null(tmp)){
    metabolites_hmdb = append(metabolites_hmdb,NA)
  }else{metabolites_hmdb = append(metabolites_hmdb,paste(tmp,collapse = ","))}
} #5835
metabolites_kegg = c()
for(i in c(1:length(metabolites)))
{
  tmp = metabolites[[i]]$annotation$kegg.compound
  if(is.null(tmp)){
    metabolites_kegg = append(metabolites_kegg,NA)
  }else{metabolites_kegg = append(metabolites_kegg,paste(tmp,collapse = ","))}
} #5835
# 其他几个数据库的名称没有做，方法同上
metabolites_recon3D = data.frame(id = metabolites_id,
                                name = metabolites_name,
                                compartment = metabolites_compartment,
                                charge = metabolites_charge,
                                formula = metabolites_formula,
                                notes = metabolites_notes,
                                biocyc = metabolites_biocyc,
                                hmdb = metabolites_hmdb,
                                kegg = metabolites_kegg)
save(metabolites_recon3D,file = "metabolites_recon3D.RData")
write.csv(metabolites_recon3D,file = "metabolites_recon3D.csv")
# 得到了metabolites 所有有效信息

# 符号除了空格有2247个，但是不重复的基因有1866个
genes = result$genes
genes_id = sapply(genes,function(x) x$id) #2248
genes_name = sapply(genes,function(x) x$name) #2248
genes_ncbigene = sapply(genes,function(x) return(ifelse((!is.null(x$annotation$ncbigene)),x$annotation$ncbigene,'unknown')))

genes_recon3D = data.frame(id = genes_id,
                          name = genes_name,
                          ncbigene = genes_ncbigene)


save(genes_recon3D,file = "genes_recon3D.RData")
write.csv(genes_recon3D,file = "genes_recon3D.csv")
# 得到了genes 所有有效信息
