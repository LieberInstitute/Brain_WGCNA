###############################################################################################
#                 Script to estimate beta for WGCNA based on connecitivity match              #
#                      "Regional-coexpression" and "Age-period" studies                       #
#                                   (NC samples only)                                         #
###############################################################################################



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
CPU = 2
registerDoParallel(cores=CPU)
getDoParWorkers()
library(furrr)
options(mc.cores = 3, future.globals.maxSize= 3565158400)

setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/")

#Reading the list of age-parsed and non-parsed rse-objects (confounders regressed)
rse_sva_cleaned_big      = readRDS("results/rse_sva_cleaned_big_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")
rse_sva_cleaned_age_grps = readRDS("results/rse_sva_cleaned_age_grps_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")

rse_sva_cleaned_full = c(rse_sva_cleaned_big, unlist(rse_sva_cleaned_age_grps))
rse_sva_cleaned_full = rse_sva_cleaned_full[!grepl("dentate.gr__-1_100",names(rse_sva_cleaned_full))]

n.samples = unlist(map(rse_sva_cleaned_full,~dim(.)[2])) #Number of samples per set

#Get sft
get_sft_parameters = function(rse, mod, net.type, assay.name, nname){
        input_matrix = t(assays(rse)[[assay.name]]) #Genes in columns, subjects in rows
        
        
        if(mod == "pearson"){
                li_more = list(
                        corOpt = list(nThreads = CPU)
                        #Not specifying corFn defaults to pearson
                )
        }
        
        WGCNA::allowWGCNAThreads(nThreads = CPU)
        WGCNAnThreads()
        
        li.args = list(data = input_matrix, powerVector = 1:35, RsquaredCut = 0.8,verbose = 5, blockSize = 13000, moreNetworkConcepts = TRUE)
        
        ret     = do.call(WGCNA::pickSoftThreshold, c(li.args, li_more, networkType = net.type))
        saveRDS(ret, paste0(nname,"...sft.",mod,".",net.type,".rds"))
        gc()
        print(paste0(nname,"...sft.",mod,".",net.type))
        return(ret)
}


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/results/noquantile_CellPropclean_GEclean_newcuts_noPCA--neuCells/")

#Run "WGCNA::pickSoftThreshold" to get estimate of beta for each rse-object 
#pearson correlation on the ranknorm assay, signed hybrid network
assay.name = "ranknorm"
mod        = "pearson"

future::plan(multiprocess)
sft.pearson.shybrid   =  furrr::future_imap(rse_sva_cleaned_full,~ get_sft_parameters(.x, nname = .y, mod = mod , assay.name = assay.name, net.type = "signed hybrid"))
future::plan(transparent)


setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/results/noquantile_CellPropclean_GEclean_newcuts_noPCA--neuCells/")
dir      = getwd()
ls       = gtools::mixedsort(list.files(dir, pattern = ".rds"))
obj.name = substring(text = ls,first = 1,last = nchar(ls)-4)
fp       = file.path(dir,ls)
fp       = fp[grepl("hippo|dentate|caudate|dlpfc",fp)]
li = list()
for (f in seq_along(fp)){
        li[[obj.name[f]]] = readRDS(fp[f])
}
li = li[!grepl("dentate.gr__-1_100",names(li))]

li_fitIndices                = map(li,"fitIndices")           ##List of sft dataframe 


metadf = tibble(full.name = names(li))
metadf = metadf %>% mutate(
        tissue_grp     = strsplit2(full.name,"\\.\\.\\.")[,1]                         ,
        tissue         = strsplit2(tissue_grp,"\\.")[,1]                              ,
        age_grp        = strsplit2(tissue_grp,"__")[,2]                               ,
        correlation    = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,2]    ,
        network        = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,3]    ,
        samples        = n.samples[tissue_grp]                                        ,
        powerEstimate  = unlist(map(li,"powerEstimate"))[.$full.name]                 ,
        class          = ifelse(grepl("gr",full.name),"parsed","non-parsed")          ,
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
        })) %>% unnest_wider(nl) %>% arrange(correlation, network, tissue)


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
                min.conn = min(map_dbl(allowed.conn,max))
                index = map(allowed.conn,~ which.min(abs(. - min.conn)))
                
                pmap(list(allowed.pwrs,allowed.conn,index, min.conn),~{
                        pwr = ..1[..3]
                        conn = round(..2[..3],3)
                        return(dplyr::lst(pwr,conn, matched.conn = ..4))
                })},
               notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
               )  %>% unnest_wider(by_all,names_sep = ".") %>% ungroup()


mm_parsed = mm %>% group_by(correlation,network,tissue) %>% filter(class == "parsed") %>%
        mutate(by_tissue.pwr = ifelse(length(purrr::reduce(allowed.pwrs,intersect))==0, NA, min(purrr::reduce(allowed.pwrs,intersect))),
               by_tissue.conn  = pmap_dbl(list(allowed.pwrs, allowed.conn, by_tissue.pwr) ,~{
                       round(..2[which(..1 == ..3)],3)}),
               notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
        ) %>% ungroup()

mm_nonparsed = mm %>% group_by(correlation,network) %>% filter(class == "non-parsed") %>%
        mutate(by_tissue.pwr = min(purrr::reduce(allowed.pwrs,intersect)),
               by_tissue.conn  = pmap_dbl(list(allowed.pwrs, allowed.conn, by_tissue.pwr) ,~{
                       round(..2[which(..1 == ..3)],3)}),
               notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
        ) %>% ungroup()

mm = bind_rows(mm_parsed,mm_nonparsed)  %>% arrange(correlation, network, tissue)              

##Save the file with list of soft-thresholding powers (by_individual, by_tissue, by_all) and other relevant info
saveRDS(mm, "soft-thresholding-pwrs_fulltable_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")
