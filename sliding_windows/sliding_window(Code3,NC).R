###############################################################################################
#            Script to Arrange samples in sliding age windows and estimate beta for WGCNA     #
#                                      via connectivity match                                 # 
#                                                                                             #
#                                     "Sliding window" study                                  #
#                                      (NC samples only)                                      #
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
library(future)
library(ggplot2)
CPU = 2
registerDoParallel(cores=CPU)
getDoParWorkers()
library(furrr)
options(mc.cores = 2, future.globals.maxSize= 3565158400)
library(slider)


setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\results")

#Reading the list of non-parsed rse-objects (confounders regressed)
rse_sva_cleaned_big         = readRDS("rse_sva_cleaned_big_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")

#Sort each object by Age
rse_sva_cleaned_big_sorted  = map(rse_sva_cleaned_big, ~ .x[, order(.x$Age)])

#Specify number of samples per window
samples.window = 40

#Find index and name of samples for each window
ordered.windows = imap(rse_sva_cleaned_big_sorted, ~{
        rse       = .x
        rse.name  = .y
        
        sanity.check = cbind.data.frame(col = colnames(rse), RNum = rse$RNum, BRNUM = rse$BrNum, Age = rse$Age, stringsAsFactors = F)
        
        ss_samples = slider::slide(colnames(rse)           ,.complete = T,.after = samples.window-1,~.) %>% .[!map_lgl(.,is.null)]
        ss_index   = slider::slide(seq_along(colnames(rse)),.complete = T,.after = samples.window-1,~.) %>% .[!map_lgl(.,is.null)]
        
        names(ss_samples) = map2_chr(ss_samples, ss_index, ~ paste(rse.name, .y[1],.y[samples.window],.x[1],.x[samples.window], sep = "__"))
        names(ss_index)   = map2_chr(ss_samples, ss_index, ~ paste(rse.name, .y[1],.y[samples.window],.x[1],.x[samples.window], sep = "__"))
        
        return(lst(ss_samples, ss_index))
        
}) %>% transpose


#Sample age by windows (wide)
v_wide = map2_dfc(map(rse_sva_cleaned_big_sorted,~.x$Age), ordered.windows$ss_index,~{
        measure = ..1
        samples = ..2
        m = map_dfr(samples, ~ measure[.x])
}) 

#Sample age by windows (long)
v_long = v_wide %>% pivot_longer(cols = everything(),values_to = "measured.var", names_to = "tissue__set.full") %>% 
        separate(col=tissue__set.full, sep = "__|.rds", into = c("tissue","sample_index_start", "sample_index_last", "sample_ID_start", "sample_ID_last"), extra = "drop", remove = F) %>%
        mutate(tissue__set = paste(tissue, sample_index_start, sample_index_last, sep = "__"),
               test        = paste(tissue, sample_index_start, sep = "__"),
               test        = factor(test, levels = mixedsort(unique(test))))

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
        
        li.args = list(data = input_matrix, powerVector = 1:25, RsquaredCut = 0.8,verbose = 5, blockSize = 13000, moreNetworkConcepts = TRUE)	
        
        cor_default     = cor	
        cor             = WGCNA::cor      #Overwrite default cor method as WGCNA::cor instead of stats::cor #Alternative is to use session package	
        
        ret             = do.call(WGCNA::pickSoftThreshold, c(li.args, li_more, networkType = net.type))	
        ret$samples     = colnames(rse)	
        ret$age         = rse$Age
        ret$sample.name = nname
        ret$exp.mat     = assay.name
        ret$net.type    = net.type 	
        cor             = cor_default	
        
        saveRDS(ret, paste0(nname,"...sft.",mod,".",net.type,".rds"))	
        #gc()	
        print(paste0(nname,"...sft.",mod,".",net.type))	
        
        return(NULL)	
}

setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins/Sliding_Window_Networks/sft.sliding.window")

#Run "WGCNA::pickSoftThreshold" to get estimate of beta for each rse-object 
#pearson correlation on the ranknorm assay, signed hybrid network
assay.name = "ranknorm"
mod        = "pearson"

setup.sft = function(rse_li, rse_windows_li){
  walk2(rse_li, rse_windows_li,~{
                rse         = ..1
                li_window   = ..2
                
                tmp = imap(li_window, ~{
                        window.name = ..2
                        sub.rse = rse[,.x]
                        print(paste0(window.name,": median age-", median(sub.rse$Age)))
                        sft.pearson.shybrid   =  get_sft_parameters(sub.rse, nname = window.name, mod = mod , assay.name = assay.name, net.type = "signed hybrid")
                        return(paste0(window.name,": median age-", median(sub.rse$Age)))
                })
                future::plan(transparent)
        })
}

setup.sft(rse_sva_cleaned_big_sorted, ordered.windows$ss_samples)
###################################################################
###################################################################

n.samples = samples.window #Number of samples per set

setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Sliding_Window_Networks\\sft.sliding.window\\")

dir      = getwd()
ls       = list.files(dir, pattern = ".rds")  %>% gtools::mixedsort(.)
li.full  = map(ls,readRDS) %>% set_names(ls)
li.full  = li.full[!grepl("random", names(li.full))]

outliers_li = split(li.full, (strsplit2(names(li.full),"__")[,1])) %>% map(~{

        mn1 = map_df(.x,~ return(round(-sign(.x$fitIndices$slope) * .x$fitIndices$SFT.R.sq,2)))
        dd1 = dist(t(mn1))
        cc1 = rowSums(as.matrix(dd1))  %>% scale %>% abs
        dimnames(cc1)[[1]][cc1>3]

        mn2 = map_df(.x,~ .x$fitIndices$median.k.)
        dd2 = dist(t(mn2))
        cc2 = rowSums(as.matrix(dd2))  %>% scale %>% abs
        dimnames(cc2)[[1]][cc2>3]
        
        return(unique(c(dimnames(cc1)[[1]][cc1>3],dimnames(cc2)[[1]][cc2>3])))

}) %>% unlist %>% unname 

li            = li.full[!(names(li.full) %in% outliers_li)]
li = li.full

li_fitIndices = map(li,"fitIndices")           ##List of sft dataframe 


metadf = data.frame(full.name = names(li), stringsAsFactors = F) %>%
        separate(col=full.name,sep = "__|\\.\\.\\.|.rds", into = c("tissue","sample_index_start", "sample_index_last", "sample_ID_start", "sample_ID_last"), extra = "drop", remove = F, convert = T)


metadf = metadf %>% mutate(
        short.name         = strsplit2(full.name,"\\.\\.\\.")[,1],
        median.age         = colMedians(as.matrix(v_wide))[match(short.name, colnames(v_wide))],
        mean.age           = colMeans(as.matrix(v_wide))  [match(short.name, colnames(v_wide))],
        correlation        = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,2]    ,
        network            = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,3]    ,
        samples            = n.samples                                                    , #map_dbl(li,~ length(.x$samples))
        powerEstimate      = unlist(map(li,"powerEstimate"))[.$full.name]                 ,
        class              = "parsed"                                                     ,
        scalefreeIndex     = map(li_fitIndices, ~{return(round(-sign(.x$slope) * .x$SFT.R.sq,2))}) ,
        #min.scalefreeIndex = map(li_fitIndices, ~{return(min(round(-sign(.x$slope) * .x$SFT.R.sq,2)))}) ,
        nl = map(li_fitIndices, ~{
                scalefreeIndex = -sign(.x$slope) * .x$SFT.R.sq                                   
                
                if(any(scalefreeIndex > 0.8)){
                        allowed.pwrs   = .x$Power[scalefreeIndex > 0.8]                   
                        allowed.conn   = round(.x$median.k.[scalefreeIndex > 0.8], 3)
                        
                } else {                              #No pwr makes network scalefree
                        allowed.conn   = .x$median.k.[which.min(abs(scalefreeIndex - 0.8))]
                        #allowed.pwrs   = .x$Power[which.min(abs(scalefreeIndex - 0.8))]                   
                        allowed.pwrs   = .x$Power[which(.x$median.k. == allowed.conn)]                   
                        
                }
                
                return(dplyr::lst(allowed.pwrs, allowed.conn))
        })) %>% unnest_wider(nl) %>% arrange(correlation, network, tissue)



tissue = "dlpfc"
pp = metadf[metadf$tissue == tissue,]
par(mfrow = c(2,1), mar = c(4,3,2,2))
plot(pp$median.age, sapply(pp$allowed.conn, max), type = "o", col = "blue", ylab = "Max. scale free connectivity", xlab = "Median Age", cex.lab = 0.7, main = tissue);grid()
plot(pp$median.age, sapply(pp$allowed.pwrs, min), type = "o", col = "red" , ylab = "Max. scale free connectivity", xlab = "Median Age", cex.lab = 0.7, main = tissue);grid()


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


dd = map_dfc(li,~ {.$fitIndices$SFT.R.sq * -sign(.$fitIndices$slope)})
dd = dd[ ,grep("pearson.signed hybrid",colnames(dd))]
matplot(1:25,dd , pch = letters[1:15],type = "b", lty = 1, col = interaction(metadf$tissue), main = "Pearson,signed hybrid"); abline(h=0.8)


saveRDS(mm, "soft-thresholding-pwrs_fulltable_Neuclean_sliding_window.rds")
saveRDS(li, "soft-thresholding-list_Neuclean_sliding_window.rds")
