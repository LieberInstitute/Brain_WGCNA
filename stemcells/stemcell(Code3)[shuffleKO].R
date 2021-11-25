###############################################################################################
#                             Script to create shoffle/KO networks                            #
#                Script to Estimate beta for WGCNA via connectivity match                     #
#                                  iPSC data[shuffle/KO]                                      #
###############################################################################################


#Get sft
get_sft_parameters = function(rse, mod, net.type, assay.name, nname){
  input_matrix = t(assays(rse)[[assay.name]]) #Genes in columns, subjects in rows
  
  #browser()
  if(mod == "pearson"){
    li_more = list(
      corOpt = list(nThreads = CPU)
      #Not specifying corFn defaults to pearson
    )
  }
  
  WGCNA::allowWGCNAThreads(nThreads = CPU)
  WGCNAnThreads()
  
  li.args = list(data = input_matrix, powerVector = 1:35, RsquaredCut = 0.8,verbose = 5, blockSize = 3000, moreNetworkConcepts = TRUE)
  cor = WGCNA::cor
  ret     = do.call(WGCNA::pickSoftThreshold, c(li.args, li_more, networkType = net.type))
  cor = stats::cor
  saveRDS(ret, paste0(nname,"...sft.",mod,".",net.type,"_protectDx_noConsensus(redone).rds"))
  gc()
  print(paste0(nname,"...sft.",mod,".",net.type))
  return(ret)
}


library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(SummarizedExperiment)
library(purrr)
#options(stringsAsFactors = FALSE)
set.seed(123)
library(dendextend)
library(limma)
library(tidyr)
library(dplyr)
library(gtools)       #Optional. Just to sort path of files read
library(doParallel)
CPU = 3
registerDoParallel(cores=CPU)
getDoParWorkers()
library(furrr)
options(mc.cores = 3, future.globals.maxSize= 3565158400)
library(magrittr)


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/")
#setwd("C:/Users/Madhur/OneDrive - Johns Hopkins")
rse_sva_cleaned_full      = readRDS("redone/results/rse_sva_cleaned_protectDx(redone).rds")


n.samples = unlist(map(rse_sva_cleaned_full,~dim(.)[2])) #Number of samples per set

our_bins = c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")

ne = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\PGC.gene.lists_grch38[125 PGC ensemblIDs].RData", envir = ne)
li_consensus = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/SNPextension/consensus_genes/consensus_gene_list.rds")
consensus_noPGC_genes = setdiff(c(li_consensus$SCZ_prenatal_positive, li_consensus$SCZ_postnatal_positive), unlist(ne$PGC3.all.biotypes[our_bins]))

##Load set of genes  to remove from Leo
nf = new.env()
load("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/results/permuted.KO.genes.Leo.RData", envir = nf)
names(nf$permuted_consensus_genes) = paste0("sim",1:100)

seed = 11072021 #Old seed: 10142021
shuffle.expression = function(rse, genes) {
  exp = assays(rse)$ranknorm
  genes_sub = rownames(rse)[strsplit2(rownames(rse),"\\.")[,1] %in% genes]
  
  for (i in genes_sub){
    set.seed(seed)
    exp[i,] = sample(exp[i,])
    seed<<- seed+1
  }
  assays(rse)$ranknorm = exp
  return(rse)
}
rse_sva_cleaned_new = list()
rse_sva_cleaned_new[["sim0"]] = shuffle.expression(rse_sva_cleaned_full$rse_new , consensus_noPGC_genes)
for (sim in names(nf$permuted_consensus_genes)){
  rse_sva_cleaned_new[[sim]] = shuffle.expression(rse_sva_cleaned_full$rse_new , nf$permuted_consensus_genes[[sim]])
}


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/redone/results/KO sft estimation shuffeled/")
library(furrr)

assay.name          = "ranknorm"
mod                 = "pearson"
plan(list(tweak(multisession, workers = 3), tweak(multisession, workers = 3)))
sft.pearson.shybrid =  future_imap(rse_sva_cleaned_new,~ get_sft_parameters(.x, nname = .y, mod = mod , assay.name = assay.name, net.type = "signed hybrid"))
plan(sequential)

# 
# ##Reference for following code block: https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# ## Plot sft parameters ----
# sft = sft.pearson.shybrid$rse_new$fitIndices
# par(mfrow = c(1,2));
# cex1 = 0.9;
# 
# plot(sft[, 1], -sign(sft[, 3]) * sft[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", main = paste("Scale independence"), type = "n" , ylim = c(-1.5,1.5))
# text(sft[, 1], -sign(sft[, 3]) * sft[, 2], labels = 1:35, cex = cex1, col = "red", ylim = c(-1.5,1.5))
# abline(h = 0.80, col = "red")
# 
# plot(sft[, 1], sft[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
# text(sft[, 1], sft[, 5], labels = 1:35, cex = cex1, col = "red")


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/redone/results/KO sft estimation shuffeled/")

li  = list()
dir = getwd()
ls  = gtools::mixedsort(list.files(dir, pattern = "sft.pearson.signed hybrid_protectDx_noConsensus\\(redone\\).rds"))
fp  = file.path(dir,ls)

for (l in ls){
  li[[substring(text = l,first = 1,last = nchar(l)-4)]] = readRDS(file.path(dir,l))
}


li_fitIndices                = map(li,"fitIndices")           ##List of sft dataframe 


metadf = tibble(full.name = names(li))
metadf = metadf %>% mutate(
  rse_name       = strsplit2(full.name,"\\.\\.\\.")[,1]                         , 
  correlation    = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,2]    ,
  network        = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,3]    ,
  samples        = n.samples                                        ,
  powerEstimate  = unlist(map(li,"powerEstimate"))[.$full.name]                 ,
  nl = map(li_fitIndices, ~{
    scalefreeIndex = -sign(.x$slope) * .x$SFT.R.sq                                   
    
    if(any(scalefreeIndex > 0.8)){
      allowed.pwrs   = .x$Power[scalefreeIndex > 0.8]                   
      allowed.conn   = .x$median.k.[scalefreeIndex > 0.8]
      
    } else {                              #No pwr makes network scalefree
      allowed.pwrs   = .x$Power[which.min(abs(scalefreeIndex - 0.8))]                   
      allowed.conn   = .x$median.k.[which.min(abs(scalefreeIndex - 0.8))]
    }
    
    return(dplyr::lst(allowed.pwrs, allowed.conn))
  })) %>% unnest_wider(nl) %>% arrange(correlation, network)


mm = metadf
mm = mm %>% mutate(
  by_individual = pmap(list(allowed.pwrs, allowed.conn) ,~{
    pwr = min(.x)
    conn = round(.y[which(.x == pwr)],3)
    return(dplyr::lst(pwr,conn))
  }),
  notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}) %>% unnest_wider(by_individual,names_sep = ".")

mm = mm %>% group_by(correlation,network) %>% 
  mutate(by_all = { 
    #min.conn = min(map_dbl(allowed.conn,max))
    min.conn = 1.123148
    index = map(allowed.conn,~ which.min(abs(. - min.conn)))
    
    pmap(list(allowed.pwrs,allowed.conn,index, min.conn),~{
      pwr = ..1[..3]
      conn = round(..2[..3],3)
      return(dplyr::lst(pwr,conn, matched.conn = ..4))
    })},
    notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
  )  %>% unnest_wider(by_all,names_sep = ".") %>% ungroup()


mm$network = gsub("_protectDx_noConsensus\\(redone\\)","",mm$network)

dd = map_dfc(li,~ {.$fitIndices$SFT.R.sq * -sign(.$fitIndices$slope)})
dd = dd[ ,grep("pearson.signed hybrid",colnames(dd))]
matplot(1:35,dd , pch = letters[1:15],type = "b", lty = 1, col = "red", main = "Pearson,signed hybrid"); abline(h=0.8)


saveRDS(mm, "soft-thresholding-pwrs_fulltable_protectDx_KOsims(redone).rds")
saveRDS(li, "soft-thresholding-list_prtoectDx_KOsims(redone).rds")
saveRDS(rse_sva_cleaned_new, "rse_sva_cleaned_protectDx_KOsims(redone).rds")
