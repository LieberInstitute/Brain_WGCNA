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

setwd(r"(C:\Users\mpariha1\Desktop\upload_to_google_drive\OneDrive - Johns Hopkins\New_Project\newResults)")

#Reading the list of sample-matched and confounders regressed rse-objects
rse_sva_cleaned_big      = readRDS("rse_sva_cleaned_big_newproject_noAgeSqr_noPCAremoved_redone.rds")

n.samples = map(rse_sva_cleaned_big,~dim(.)[2]) #Number of samples per set


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
        
        li.args = list(data = input_matrix, powerVector = 1:25, RsquaredCut = 0.8,verbose = 5, blockSize = 9000, moreNetworkConcepts = TRUE)
        
        ret     = do.call(WGCNA::pickSoftThreshold, c(li.args, li_more, networkType = net.type))
        saveRDS(ret, paste0(nname,"...sft.",mod,".",net.type,".rds"))
        gc()
        print(paste0(nname,"...sft.",mod,".",net.type))
        return(ret)
}


setwd(r"(C:\Users\mpariha1\Desktop\upload_to_google_drive\OneDrive - Johns Hopkins\New_Project\newResults\sft_noAgeSqr_noPCAremoved_redone)")

#Run "WGCNA::pickSoftThreshold" to get estimate of beta for each rse-object 
#pearson correlation on the ranknorm assay, signed hybrid network
assay.name = "ranknorm"
mod        = "pearson"
sft.pearson.shybrid   =  imap(rse_sva_cleaned_big,~ get_sft_parameters(.x, nname = .y, mod = mod , assay.name = assay.name, net.type = "signed hybrid"))


setwd(r"(C:\Users\mpariha1\Desktop\upload_to_google_drive\OneDrive - Johns Hopkins\New_Project\newResults\sft_noAgeSqr_noPCAremoved_redone)")
dir      = getwd()
ls       = gtools::mixedsort(list.files(dir, pattern = ".rds"))
obj.name = substring(text = ls,first = 1,last = nchar(ls)-4)
fp       = file.path(dir,ls)

li = list()
for (f in seq_along(fp)){
        li[[obj.name[f]]] = readRDS(fp[f])
}
li = li[!grepl("dentate.gr__-1_100",names(li))]
li_fitIndices               = map(li,"fitIndices")

metadf = tibble(full.name = names(li))
metadf = metadf %>% mutate(
        tissue_grp     = strsplit2(full.name,"\\.\\.\\.")[,1]                         ,
        tissue         = strsplit2(tissue_grp,"\\.")[,1]                              ,
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


mm = metadf %>% group_by(correlation,network) %>% 
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


# dd = map_dfc(li,~ {.$fitIndices$SFT.R.sq * -sign(.$fitIndices$slope)})
# dd = dd[ ,grep("pearson.signed hybrid",colnames(dd))]
# matplot(1:25,dd , pch = letters[1:15],type = "b", lty = 1, col = interaction(metadf$tissue), main = "Pearson,signed hybrid"); abline(h=0.8)

##Save the file with list of soft-thresholding powers and other relevant info
saveRDS(mm, "soft-thresholding-pwrs_fulltable__newproject_noAgeSqr_noPCAremoved_redone.rds")
