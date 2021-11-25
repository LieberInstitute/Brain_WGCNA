####################################################################################################
#                   Script to erform WGCNA to obtain samples match HP-DG networks                  #
#                               "Cell-population enrichment" studies                               #
#                                        [QSVA pipeline]                                           #
####################################################################################################


library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(SummarizedExperiment)
library(purrr)

set.seed(123)
library(dendextend)
library(limma)
library(dplyr)
library(pheatmap)
library(gtools)       #Optional. Just to sort path of files read
library(doParallel)
CPU = 2
registerDoParallel(cores=CPU)
getDoParWorkers()


setwd(r"(C:\Users\mpariha1\Desktop\upload_to_google_drive\OneDrive - Johns Hopkins\New_Project\newResults)")

#Reading the list of sample-matched and confounders+qsv regressed rse-objects
rse_sva_cleaned_big      = readRDS("rse_sva_cleaned_big_newproject_noAgeSqr_noPCAremoved_qSVAremoved.rds")

#Reading the metafile with estimated power for each rse-object
pwr.df  = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\New_Project\\newResults\\sft_noAgeSqr_noPCAremoved_qSVAremoved\\soft-thresholding-pwrs_fulltable__newproject_noAgeSqr_PCAremoved_qSVAremoved.rds")


n.samples = unlist(map(rse_sva_cleaned_big,~dim(.)[2])) #Number of samples per set
n.genes   = unlist(map(rse_sva_cleaned_big,~dim(.)[1])) #Number of genes   per set

get_networks = function(rse,mod, net.type, assay.name, beta.pwr,pwr.type, nname){
  if(mod == "pearson"){
    li_more = list(nThreads = CPU, corType = "pearson")
    #Not specifying corFn defaults to pearson
  }

  WGCNA::allowWGCNAThreads(nThreads = CPU)
  WGCNAnThreads()
  
  pwr = unlist(beta.pwr[beta.pwr$tissue_grp == nname & beta.pwr$correlation == mod & beta.pwr$network == net.type, paste0(pwr.type,".pwr")])
  
  input_matrix = t(assays(rse)[[assay.name]]) #Genes in columns, subjects in rows
  
  li.args = list(datExpr = input_matrix, randomSeed = 123, verbose = 5, maxBlockSize = 9000,
                 TOMType = "signed", saveTOMs = FALSE, saveTOMFileBase = (nname), power = pwr,
                 minModuleSize = 40,
                 pamStage = TRUE, pamRespectsDendro = TRUE,  deepSplit = 4,                 
                 reassignThreshold = 1e-06,                                                 
                 mergeCutHeight    = 0.15,
                 numericLabels     = TRUE)
  
    
  print(nname)
  print(pwr)
  cor_default    = cor
  cor            = WGCNA::cor      #Overwrite default cor method as WGCNA::cor instead of stats::cor
  net            = do.call("blockwiseModules"                 , c(li.args    , li_more, networkType = net.type))
  
  li.args.IMC    = list(datExpr =  input_matrix, colors = net$colors, power = pwr, getWholeNetworkConnectivity = TRUE)
  IMC            = do.call("intramodularConnectivity.fromExpr", c(li.args.IMC, networkType = net.type))
  IMC$modules    = net$colors
  rownames(IMC)  = names(net$colors)
  net$IMC        = IMC
  
  net$KME        = WGCNA::signedKME(input_matrix, net$MEs)
  
  cor            = cor_default     #Reassign default cor method
  
  print(paste0(nname,"...sft.",mod,".",net.type,"|pwr:",pwr,"|pwr.type: ",pwr.type))
  saveRDS(net, paste0(nname,"...sft.",mod,".",net.type," & pwr-",pwr," & pwr.type-",pwr.type, ".rds"))
  #gc()
  return(net)
}

setwd(r"(C:\Users\mpariha1\Desktop\upload_to_google_drive\OneDrive - Johns Hopkins\New_Project\newResults\sft_noAgeSqr_noPCAremoved_qSVAremoved\net)")
#Run "WGCNA::blockwiseModules" to get network
#pearson correlation on the ranknorm assay, signed hybrid network, "by_all" power
assay.name = "ranknorm"
pwr.type   = "by_all"
mod        = "pearson"
net.pearson.shybrid   =  imap(rse_sva_cleaned_big,~ get_networks(.x, nname = .y, mod = mod , assay.name = assay.name, beta.pwr = pwr.df, pwr.type = pwr.type, net.type = "signed hybrid"))

dir      = getwd()
ls       = gtools::mixedsort(list.files(dir, pattern = ".rds"))
obj.name = substring(text = ls,first = 1,last = nchar(ls)-4)
fp       = file.path(dir,ls)

li     = list()
for (f in seq_along(fp)){
  li[[obj.name[f]]]     = readRDS(fp[f])
}

###Save blockwise module output for functional analysis
gene_module_list = imap(li,~ {
  all_modules = data.frame(module_numbers = .$colors, module_colors = WGCNA::labels2colors(.$colors), gene = names(.$colors), gene.clean = strsplit2(names(.$colors),"\\.")[,1], stringsAsFactors = F)
  module_li   = split(all_modules$gene.clean,all_modules$module_colors)
  return(module_li)
})

###Save Intramodular connectivity as list
IMC_list = imap(li,"IMC")

###Save KME as list
KME_list = imap(li,"KME")


############Get metainfo from list li
metadf = tibble(li.name = names(li))
metadf = metadf %>% tidyr::separate(li.name,into = c("full.name", "pwr.used", "match.used"),sep = " & ", remove = F)
metadf$match.used = gsub("pwr.type","",metadf$match.used)

metadf = merge(metadf,pwr.df, all.x = T, by = "full.name")  #Changes the order of x


metadf = metadf %>% mutate(
  genes         = n.genes[tissue_grp],
  grey          = map_int(li[match(li.name,names(li))],~table(.$colors)["0"]),
  N.mm          = map_int(li[match(li.name,names(li))],~length(table(.$colors))),
  modules       = map(li[match(li.name,names(li))],~table(.$colors)),           #This is a list-column
  modules_flat  = map_chr(modules,~paste0(.,collapse = "__")),                  #Flattened to character
  #modules       = NULL,
  pct.grey      = grey*100/genes,
  MEs           = map(li[match(li.name,names(li))], "MEs")
)

gene_module_metadf = metadf %>% mutate(test = map(gene_module_list[match(.$li.name,names(gene_module_list))],~.),
                                       IMC  = IMC_list,
                                       KME  = KME_list
                                       )

##Save the file with list of network metainfo (number of genes, number of grey, number of modules, power used etc)
saveRDS(gene_module_metadf, r"(C:\Users\mpariha1\Desktop\upload_to_google_drive\OneDrive - Johns Hopkins\New_Project\newResults\sft_noAgeSqr_noPCAremoved_qSVAremoved\net\metadf_noAgeSqr_noPCAremoved_qSVAremoved_deepSplit4_byAll.rds)")
