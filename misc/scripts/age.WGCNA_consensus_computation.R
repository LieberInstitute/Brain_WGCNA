
enriched.modules = readRDS("age.WGCNA.paper/consensus/new_consensus_raw_data_metanets.rds")
colnames(enriched.modules$risk.module.negative.corrected)[1:2] = c("new.network","module")
perinatal = c("DLPFC Perinatal",
              # "META Perinatal",
              "Werling2020",
              "Walker2019",
              "Li2018")
postnatal = c("DLPFC Juvenile",
              # "META Juvenile",
              "DLPFC Adult",
              "DLPFC Older Adult",
              "Fromer2016_case",
              "Fromer2016_control",
              "Gandal2018a",
              "Pergola2017",
              "Pergola2019",
              "Pergola2020",
              "Radulescu2020",
              "Hartl2021")

################ positive list
### prenatal
net.mod = enriched.modules$risk.module.positive
net = enriched.modules$gene_module_list_PC
net.mod = net.mod[net.mod$new.network %in% perinatal,]
net = net[perinatal]
consensus <- unlist(sapply(names(net), function(x){
  mod = net.mod$module[net.mod$new.network %in% x]
  net[[x]][mod]
}))
data <- as.data.frame(sort(table(factor(consensus))[table(factor(consensus))>= 3], decreasing = TRUE))
colnames(data) = c("ensembl","hits")
data$ensembl = as.character(data$ensembl)
data$tot <- sapply(data$ensembl,function(x){
  net.check = lapply(net,unlist)
  sum(sapply(net.check,function(y) x %in% y))
})
data$ratio <- data$hits/data$tot
consensus.prenatal = data[data$tot >= 3 & signif(data$ratio,2) >= 0.67,]
consensus.prenatal[,perinatal] = t(sapply(consensus.prenatal$ensembl,function(x){
  sapply(perinatal,function(y){ names(net[[y]])[sapply(net[[y]],function(z) x %in% z)] })
}))
### postnatal
net.mod = enriched.modules$risk.module.positive
net.mod$new.network = gsub("BRNCTX","Hartl2021",net.mod$new.network)
net = enriched.modules$gene_module_list_PC
names(net)[grep("BRNCTX",names(net))] = "Hartl2021"
net.mod = net.mod[net.mod$new.network %in% postnatal,]
net = net[postnatal]
consensus <- unlist(sapply(names(net), function(x){
  mod = net.mod$module[net.mod$new.network %in% x]
  net[[x]][mod]
}))
data <- as.data.frame(sort(table(factor(consensus))[table(factor(consensus))>= 3], decreasing = TRUE))
colnames(data) = c("ensembl","hits")
data$ensembl = as.character(data$ensembl)
data$tot <- sapply(data$ensembl,function(x){
  net.check = lapply(net,unlist)
  sum(sapply(net.check,function(y) x %in% y))
})
data$ratio <- data$hits/data$tot
consensus.postnatal = data[data$tot >= 6 & signif(data$ratio,2) >= 0.67,]
consensus.postnatal[,postnatal] = t(sapply(consensus.postnatal$ensembl,function(x){
  sapply(postnatal,function(y){ 
    result = names(net[[y]])[sapply(net[[y]],function(z) x %in% z)] 
    if(length(result)) result else ""
  })
}))

#################### negative list
### prenatal
net.mod = enriched.modules$risk.module.negative.corrected
net = enriched.modules$gene_module_list_PC
net.mod = net.mod[net.mod$new.network %in% perinatal,]
net = net[perinatal]
consensus <- unlist(sapply(names(net), function(x){
  mod = net.mod$module[net.mod$new.network %in% x]
  net[[x]][mod]
}))
data <- as.data.frame(sort(table(factor(consensus))[table(factor(consensus))>= 3], decreasing = TRUE))
colnames(data) = c("ensembl","hits")
data$ensembl = as.character(data$ensembl)
data$tot <- sapply(data$ensembl,function(x){
  net.check = lapply(net,unlist)
  sum(sapply(net.check,function(y) x %in% y))
})
data$ratio <- data$hits/data$tot
consensus.prenatal.negative = data[data$tot >= 3 & signif(data$ratio,2) >= 0.67,]
consensus.prenatal.negative[,perinatal] = t(sapply(consensus.prenatal.negative$ensembl,function(x){
  sapply(perinatal,function(y){ 
    result = names(net[[y]])[sapply(net[[y]],function(z) x %in% z)] 
    if(length(result)) result else ""
  })
}))
### postnatal
net.mod = enriched.modules$risk.module.negative.corrected
net.mod$new.network = gsub("BRNCTX","Hartl2021",net.mod$new.network)
net = enriched.modules$gene_module_list_PC
names(net)[grep("BRNCTX",names(net))] = "Hartl2021"
net.mod = net.mod[net.mod$new.network %in% postnatal,]
net = net[postnatal]
consensus <- unlist(sapply(names(net), function(x){
  mod = net.mod$module[net.mod$new.network %in% x]
  net[[x]][mod]
}))
data <- as.data.frame(sort(table(factor(consensus))[table(factor(consensus))>= 3], decreasing = TRUE))
colnames(data) = c("ensembl","hits")
data$ensembl = as.character(data$ensembl)
data$tot <- sapply(data$ensembl,function(x){
  net.check = lapply(net,unlist)
  sum(sapply(net.check,function(y) x %in% y))
})
data$ratio <- data$hits/data$tot
consensus.postnatal.negative = data[data$tot >= 6 & signif(data$ratio,2) >= 0.67,]
consensus.postnatal.negative[,postnatal] = t(sapply(consensus.postnatal.negative$ensembl,function(x){
  sapply(postnatal,function(y){ 
    result = names(net[[y]])[sapply(net[[y]],function(z) x %in% z)] 
    if(length(result)) result else ""
  })
}))

save(consensus.prenatal,consensus.prenatal.negative,consensus.postnatal,consensus.postnatal.negative,
     file = "age.WGCNA.paper/consensus/consensus_withWalker_0.6thrershold_PC.RData")



