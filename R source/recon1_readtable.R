





# 第一步 用fbar生成reaction stoich metabolite三个表
# 第二步 进行本程序得到边表 和节点表
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Metabolic Network"
setwd(workingdir)
memory.limit(400000)

#stoich
stoich = read.csv(file = "stoich.csv")



reactions_df = read.csv(file = "reactions.csv" )
fix(stoich)
# numR = length(unique(stoich$abbreviation))
# numM = length(unique(stoich$met))

# 初筛这个矩阵,去掉currency metabolites
#proton, water, oxygen molecule, NADP+ , NADPH, ATP, 
# diphosphate(ppi), carbon dioxide, phosphate, ADP, coA,
# UDP, NAD+ , NADH, AMP, ammonia, hydrogen peroxide,
# oxidized electron acceptor, reduced electron acceptor, 
# 3-5-ADP, GDP, carbon monoxide, GTP, and FAD
currency = c("h","h2o","o2","nadp","nadph","atp","adp","pi","ppi","co2","coa","udp","nad","nadh","amp",
             "nh3","nh4","h2o2","gdp","fad","gtp","co")
compartment = c("c","e","g","l","m","n","r","x")
currency_m = c()
for(i in currency)
{
  currency_m = append(currency_m,paste0(i,"_",compartment))
}
stoich = stoich[!stoich$met %in% currency_m,]

# 初筛这个矩阵，去掉不需要酶的反应
# 酶促反应有哪些呢
en_re = reactions_df$abbreviation[reactions_df$geneAssociation != '']
stoich = stoich[stoich$abbreviation %in% en_re,]
table(data.frame(table(stoich$abbreviation))$Freq)

# 现在开始做met to met 的表格
reactions = unique(stoich$abbreviation)
edge_data = data.frame()
for(i in c(1:length(reactions)))
{
  reaction = reactions[i]
  submatrix = stoich[stoich$abbreviation == reaction,]
  met_out = submatrix$met[submatrix$stoich<0]
  met_in = submatrix$met[submatrix$stoich>0]
  from = c()
  to = c()
  for(m in met_out)
  {
    for(n in met_in)
    {
      from = append(from,m)
      to = append(to,n)
    }
  }
  sub_data = data.frame(from = from,to = to,reaction = reaction)
  if(dim(edge_data)[1]==0)
  {
    edge_data = sub_data
  }else{
    edge_data = rbind(edge_data,sub_data)
  }
}

table(data.frame(table(edge_data$from))$Freq)
table(data.frame(table(edge_data$to))$Freq)

met_data = unique(c(edge_data$from,edge_data$to))
met_data = data.frame(id = c(0:(length(met_data)-1)),name = met_data)
src = met_data$id[match(edge_data$from,met_data$name)]
dst = met_data$id[match(edge_data$to,met_data$name)]

edge_data = cbind(edge_data,src = src,dst = dst)

# 加一个步骤，将起点终点都是currency的边去掉
critiria1 = !edge_data$from %in% currency_m
critiria2 = !edge_data$to %in% currency_m

edge_data = edge_data[critiria1|critiria2,]

save(edge_data,file = "edge_data.RData")
write.csv(edge_data,file = "edge_data.csv")
save(met_data,file = "met_data.RData")
write.csv(met_data,file = "met_data.csv")


# 处理reactions_df的过程
reactions_df = reactions_df[!reactions_df$geneAssociation %in% c("","0"),]
reactions_df$geneAssociation = gsub('\\(', '( ', reactions_df$geneAssociation)
reactions_df$geneAssociation = gsub(')',' )',reactions_df$geneAssociation)
special_index = which(nchar(reactions_df$geneAssociation)>1000)
for(i in special_index)
{
  token = reactions_df$geneAssociation[i]
  tokens = strsplit(token,' ')[[1]]
  tokens = unique(tokens[!tokens %in% c('(',')','or','and')])
  tokens = str_c(tokens,collapse = " or ")
  reactions_df$geneAssociation[i] = tokens
}
write.csv(reactions_df,file = "enzymatic_reactions.csv")

genesid = unique(unlist(strsplit(reactions_df$geneAssociation,' '))) #2248个基因
genesid = genesid[grepl('_AT',genesid)] # 2244个基因 去掉了"or"  "("   "and" ")" 
entrez = sapply(strsplit(genesid,"_"),function(x) x[1])
library(org.Hs.eg.db)
library(annotate)
symbols = getSYMBOL(entrez, data='org.Hs.eg')
genesmap = data.frame(genesid = genesid,id = entrez,name = symbols)

write.csv(genesmap,file = "genes_recon1_map.csv")
a = read.csv(file = "genes_recon1_map.csv")
# ND3 MT-ND3
# ND6 MT-ND6
# ND4L MT-ND4L
# ND4 MT-ND4
# ND2 MT-ND2
# ND1 MT-ND1
# ND5 MT-ND5
# COX1 MT-CO1
# COX2 MT-CO2 
# COX3 MT-CO3

"GFUS" "LOC102724197" NA "LOC344967"  NA  "LOC286297" 