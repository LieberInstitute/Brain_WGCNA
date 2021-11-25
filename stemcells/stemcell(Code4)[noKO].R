library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(SummarizedExperiment)
library(purrr)
#options(stringsAsFactors = FALSE)
set.seed(123)
library(dendextend)
library(limma)
library(dplyr)
library(pheatmap)
library(gtools)       #Optional. Just to sort path of files read
library(doParallel)
CPU = 3
registerDoParallel(cores=CPU)
getDoParWorkers()
library(furrr)
options(mc.cores = 3, future.globals.maxSize= 3565158400)


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/")

rse_sva_cleaned_full      = readRDS("redone/results/rse_sva_cleaned_protectDx(redone).rds")


pwr.df = readRDS("redone/results/soft-thresholding-pwrs_fulltable_protectDx(redone).rds")

n.samples = unlist(map(rse_sva_cleaned_full,~dim(.)[2])) #Number of samples per set
n.genes   = unlist(map(rse_sva_cleaned_full,~dim(.)[1])) #Number of genes   per set


get_networks = function(rse,mod, net.type, assay.name, beta.pwr,pwr.type, nname){
  if(mod == "pearson"){
    li_more = list(nThreads = CPU, corType = "pearson")
    #Not specifying corFn defaults to pearson
  }
  
  WGCNA::allowWGCNAThreads(nThreads = CPU)
  WGCNAnThreads()
  
  pwr = unlist(beta.pwr[beta.pwr$rse_name == nname & beta.pwr$correlation == mod & beta.pwr$network == net.type, paste0(pwr.type,".pwr")])

  input_matrix = t(assays(rse)[[assay.name]]) #Genes in columns, subjects in rows
  
  li.args = list(datExpr = input_matrix, randomSeed = 123, verbose = 5, maxBlockSize = 3000,
                 TOMType = "signed", saveTOMs = FALSE, saveTOMFileBase = (nname), power = pwr,
                 minModuleSize = 40,
                 pamStage = TRUE, pamRespectsDendro = TRUE,  deepSplit = 4,                     
                 reassignThreshold = 1e-06,                                                     
                 mergeCutHeight    = 0.15,                                                      
                 numericLabels     = TRUE)
  
  print(nname)
  print(pwr)
  cor_default    = cor
  cor            = WGCNA::cor      #Overwrite default cor method as WGCNA::cor instead of stats::cor #Alternative is to use session package

  net            = do.call("blockwiseModules", c(li.args, li_more, networkType = net.type))
  

  li.args.IMC    = list(datExpr =  input_matrix, colors = net$colors, power = pwr,	
                        getWholeNetworkConnectivity = TRUE)	
  IMC            = do.call("intramodularConnectivity.fromExpr", c(li.args.IMC, networkType = net.type))	
  IMC$modules    = net$colors	
  rownames(IMC)  = names(net$colors)	
  net$IMC        = IMC
  
  net$KME        = WGCNA::signedKME(input_matrix, net$MEs)
  cor            = cor_default     #Reassign default cor method


  print(paste0(nname,"...sft.",mod,".",net.type,"|pwr:",pwr,"|pwr.type: ",pwr.type))
  saveRDS(net, paste0(nname,"...sft.",mod,".",net.type," & pwr-",pwr," & pwr.type-",pwr.type, "_protectDx(redone).rds"))
  #gc()
  return(net)
}


setwd("redone/results/")

assay.name = "ranknorm"
pwr.type   = "by_all"
mod        = "pearson"
net.pearson.shybrid   =  imap(rse_sva_cleaned_full,~ get_networks(.x, nname = .y, mod = mod , assay.name = assay.name, beta.pwr = pwr.df, pwr.type = pwr.type, net.type = "signed hybrid"))

 
 
#######################
dir      = getwd()
ls       = gtools::mixedsort(list.files(dir, pattern = "sft.pearson.signed hybrid & pwr-.*&.*._protectDx\\(redone\\).rds"))
obj.name = substring(text = ls,first = 1,last = nchar(ls)-4)
fp       = file.path(dir,ls)

li     = list()
li_rse = list()
for (f in seq_along(fp)){
  li[[obj.name[f]]]     = readRDS(fp[f])
}


###Save blockwise module output for functional analysis
gene_module_list = imap(li,~ {
  all_modules = data.frame(module_numbers = .$colors, module_colors = WGCNA::labels2colors(.$colors), gene = names(.$colors), gene.clean = strsplit2(names(.$colors),"\\.")[,1], stringsAsFactors = F)
  module_li   = split(all_modules$gene.clean,all_modules$module_colors)
  return(module_li)
})
gene_module_list1 = gene_module_list
names(gene_module_list1)[grepl("pwr-15",names(gene_module_list))] = "rse.pwr15.protectDx"
names(gene_module_list1)[grepl("pwr-7" ,names(gene_module_list))]  = "rse.pwr7.protectDx"
saveRDS(gene_module_list1,"gene-module list(stem cells)_protectDx(redone).rds")

###Save Intramodular connectivity as list	
IMC_list = imap(li,"IMC")

###Save KME as list	
KME_list = imap(li,"KME")


############Get metainfo from list li
metadf = tibble(li.name = names(li))
metadf = metadf %>% tidyr::separate(li.name,into = c("full.name", "pwr.used", "match.used"),sep = " & ", remove = F)
metadf$match.used = gsub("pwr.type|_protectDx\\(redone\\)","",metadf$match.used)
metadf$full.name = paste0(metadf$full.name, "_protectDx(redone)")
metadf = merge(metadf,pwr.df, all.x = T, by = "full.name")  #Changes the order of x


metadf = metadf %>% mutate(
  genes         = n.genes[rse_name],
  grey          = map_int(li[match(li.name,names(li))],~table(.$colors)["0"]),
  N.mm          = map_int(li[match(li.name,names(li))],~length(table(.$colors))),
  modules       = map(li[match(li.name,names(li))],~table(.$colors)),                                     #This is a list-column
  modules_flat  = map_chr(modules,~paste0(.,collapse = "__")),                  #Flattened to character
  #modules       = NULL,
  pct.grey      = grey*100/genes,
  MEs           = map(li[match(li.name,names(li))], "MEs"),
  hist_approx   = sapply(modules_flat, function(i) {
    i      = as.numeric(unlist(strsplit(i,"__")))
    bars   = vapply(9601:9608, intToUtf8, character(1))[-c(4,8)]
    factor = cut(i[-1]/1000, breaks = seq(0, 1, length.out = length(bars) + 1), labels = bars, include.lowest = TRUE) ## Absolute scaling, not relative to any group
    chars  = as.character(factor)
    chars[is.na(chars)] = bars[length(bars)]
    char_hist = paste0(paste0(chars, collapse = " "),"||<<<")
    return(char_hist)
  })
)

gene_module_metadf = metadf %>% mutate(test = map(gene_module_list[match(.$li.name,names(gene_module_list))],~.),
                                       IMC  = IMC_list[match(metadf$li.name, names(IMC_list))],
                                       KME  = KME_list[match(metadf$li.name, names(KME_list))])

saveRDS(gene_module_metadf, "metadf_networks_protectDx(redone).rds")
