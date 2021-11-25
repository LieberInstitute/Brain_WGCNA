library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(SummarizedExperiment)
library(purrr)

set.seed(123)
library(dendextend)
library(limma)
library(tidyr)
library(dplyr)
library(gtools)       #Optional. Just to sort path of files read
library(doParallel)
library(future)
library(ggplot2)
CPU = 2
registerDoParallel(cores=CPU)
getDoParWorkers()
library(furrr)
options(mc.cores = 2, future.globals.maxSize= 3565158400)
library(slider)


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Sliding_Window_Networks/")

#Reading the metafile with estimated power for each rse-object
pwr.df = readRDS("results/soft-thresholding-pwrs_fulltable_Neuclean_sliding_window.rds")
pwr.df$outlier.window = ifelse((pwr.df$by_all.conn > unique(pwr.df$by_all.matched.conn)*3) |
                                            (pwr.df$by_all.conn < unique(pwr.df$by_all.matched.conn)/3), 1, 0)

#Reading the list of non-parsed rse-objects (confounders regressed)
rse_sva_cleaned_big         = readRDS("results/rse_sva_cleaned_big_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")


rse_sva_cleaned_big_sorted  = map(rse_sva_cleaned_big, ~ .x[, order(.x$Age)])

samples.window = 40

ordered.windows = imap(rse_sva_cleaned_big_sorted, ~{
  rse       = .x
  rse.name  = .y
  
  sanity.check = cbind.data.frame(col = colnames(rse), RNum = rse$RNum, BRNUM = rse$BrNum, Age = rse$Age, stringsAsFactors = F)
  
  ss_samples = slider::slide(colnames(rse)           ,.complete = T,.after = samples.window-1,~.) %>% .[!map_lgl(.,is.null)]
  ss_index   = slider::slide(seq_along(colnames(rse)),.complete = T,.after = samples.window-1,~.) %>% .[!map_lgl(.,is.null)]
  
  names(ss_samples) = map2_chr(ss_samples, ss_index, ~ paste(rse.name, .y[1],.y[samples.window],.x[1],.x[samples.window], sep = "__"))
  names(ss_index) = map2_chr(ss_samples, ss_index, ~ paste(rse.name, .y[1],.y[samples.window],.x[1],.x[samples.window], sep = "__"))
  
  return(lst(ss_samples, ss_index))
  
}) %>% transpose

n.samples = samples.window



get_networks = function(rse,mod, net.type, assay.name, beta.pwr,pwr.type, nname){
  
  if(mod == "pearson"){
    li_more = list(nThreads = CPU, corType = "pearson")
    #Not specifying corFn defaults to pearson
  }
  
  WGCNA::allowWGCNAThreads(nThreads = CPU)
  WGCNAnThreads()
  
  pwr = unlist(beta.pwr[beta.pwr$short.name == nname & beta.pwr$correlation == mod & beta.pwr$network == net.type, paste0(pwr.type,".pwr")])
  
  input_matrix = t(assays(rse)[[assay.name]]) #Genes in columns, subjects in rows
  
  li.args = list(datExpr = input_matrix, randomSeed = 123, verbose = 5, maxBlockSize = 13000,   
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
  
  li.args.IMC    = list(datExpr =  input_matrix, colors = net$colors, power = pwr,getWholeNetworkConnectivity = TRUE)	
  IMC            = do.call("intramodularConnectivity.fromExpr", c(li.args.IMC, networkType = net.type))	
  IMC$modules    = net$colors	
  rownames(IMC)  = names(net$colors)	
  
  net$IMC        = IMC
  net$KME        = WGCNA::signedKME(input_matrix, net$MEs)
  cor            = cor_default     #Reassign default cor method
   
  net$samples     = colnames(rse)
  net$age         = rse$Age
  net$sample.name = nname
  net$exp.mat     = assay.name
  net$net.type    = net.type
  net$pwr.used    = pwr 	
  net$match.used  = pwr.type
  
  print(paste0(nname,"...sft.",mod,".",net.type,"|pwr:",pwr,"|pwr.type: ",pwr.type))
  saveRDS(net, paste0(nname,"...sft.",mod,".",net.type," & pwr-",pwr," & pwr.type-",pwr.type, ".rds"))
  #gc()
  return(NULL)
}


setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Sliding_Window_Networks\\results\\net")
#Run "WGCNA::blockwiseModules" to get network
#pearson correlation on the ranknorm assay, signed hybrid network, "by_all" power


#future::plan(multiprocess)

pwr.type   = "by_all"
assay.name = "ranknorm"
mod        = "pearson"

future::plan(transparent)	
setup.net = function(rse_li, rse_windows_li, pwr.df){	
  walk2(rse_li, rse_windows_li,~{	
    rse         = ..1	
    li_window   = ..2	

      tmp = imap(li_window, ~{	
      window.name = ..2	
      sub.rse = rse[,.x]	
      print(paste0(window.name,": median age-", median(sub.rse$Age)))	
      
      net.pearson.shybrid   =  get_networks(sub.rse, nname = window.name, mod = mod , assay.name = assay.name, beta.pwr = pwr.df, pwr.type = pwr.type, net.type = "signed hybrid")	
      
      return(paste0(window.name,": median age-", median(sub.rse$Age)))	
    })	
    future::plan(transparent)	
  })	
}	

setup.net(rse_sva_cleaned_big_sorted, ordered.windows$ss_samples, pwr.df = pwr.df)	


#################
###########################
#####################################

dir = getwd()
ls  = list.files(dir, pattern = ".rds")  %>% gtools::mixedsort(.)
ls = ls[!grepl("metadf|big_list",ls)]

metadf = data.frame(file.name = ls, stringsAsFactors = F) %>%
  separate(col=file.name,sep = "\\.\\.\\.", into = c("short.name"), extra = "drop", remove = F, convert = T)

metadf = merge(metadf, pwr.df, by = "short.name", all.x = T) %>% dplyr::select(file.name, full.name, short.name, everything()) 
metadf = metadf[mixedorder(metadf$short.name),]

metadf = metadf %>% filter(outlier.window == 0)

li  = imap(metadf$file.name,~{print(.x);print(.y);readRDS(.x)}) %>% set_names(metadf$file.name)

metadf$KTotal  = map(li,~{.x$IMC[,"kTotal", drop = F]}) %>% set_names(metadf$short.name)
metadf$KWithin = map(li,~{.x$IMC[,"kWithin", drop = F]}) %>% set_names(metadf$short.name)

metadf$test    = map(li,"colors") %>% set_names(metadf$short.name)

big_list = map(metadf$test, ~{
  xx = data.frame(modules = .x, genes = names(.x), stringsAsFactors = F) %>% 
    mutate(genes = strsplit2(genes,"\\.")[,1], modules = WGCNA::labels2colors(modules))
  
  zz = split(xx$genes, xx$modules)
})


#Add back sample IDs for each window
samples.info.to.add.back = purrr::flatten(ordered.windows$ss_samples)
metadf$sampleInfo = samples.info.to.add.back[match(metadf$short.name, names(samples.info.to.add.back))]

saveRDS(metadf,file   = "metadf (noquantile, Neuclean, sliding, NC).rds")
saveRDS(big_list,file = "big_list (noquantile, Neuclean, sliding, NC).rds")

####
library(biomaRt)
library(tidyverse)

set.seed(123)

geneMap_fun = function(genes, ensembl = "grch37.ensembl.org", attributes){
  require(biomaRt)
  ensembl_archive = useMart("ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl", host=ensembl)
  Filters = "ensembl_gene_id"
  
  df = getBM(attributes = attributes,Filters,values = genes,mart = ensembl_archive)
  df = df[df$chromosome_name %in% c(1:22,"X","Y","M"),]
  df$chromosome_name = factor(df$chromosome_name,levels = c(1:22,"X","Y","M"))
  df = df[order(df$chromosome_name,df$start_position,df$end_position),]
  return(df)
}

all.genes = unique(unlist(big_list))

gene_grch38 = geneMap_fun(genes = all.genes, ensembl = "jul2016.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gc_content","entrezgene","version"))
gene_grch38$strand = ifelse(gene_grch38$strand > 0, "+","-")

all.genes_grch38 = unique(gene_grch38$ensembl_gene_id)
big_list_grch38 = map(big_list,~ {
  map(.x, ~{intersect(.x,all.genes_grch38)})
})


saveRDS(big_list_grch38,file = "big_list (noquantile, Neuclean, sliding, NC)[grch38].rds")