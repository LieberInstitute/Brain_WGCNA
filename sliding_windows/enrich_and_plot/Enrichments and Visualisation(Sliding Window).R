#Function to perform hypergeometric enrichment on genesets
Enrich = function(ll.target, ll.components){
  ###ll.components.size argument is not needed
  
  all.genes = unique(unlist(ll.components)) 
  all.genes.except.greys = unique(unlist(ll.components[!(names(ll.components) %in% "grey")]))   #All genes except Grey genes (already exlcuded)
  
  
  ## create list of dataframes for hypergeometric test
  
  enrich.df.list = sapply(ll.target, ll.components = ll.components, function(target, ll.components) {
    
    hit.genes     = sapply(ll.components, intersect, target)
    hit.count     = sapply(hit.genes    , length)
    hit.pop       = length(intersect(all.genes.except.greys, target))
    fail.pop      = length(all.genes.except.greys) - hit.pop
    drawn         = sapply(ll.components,length)
    df            = data.frame(hit.count = hit.count ,
                               hit.pop    = hit.pop   ,
                               fail.pop   = fail.pop  ,
                               drawn      = drawn     )
    
    #   hit.pop.grey  = length(intersect(all.genes , target))
    #   fail.pop.grey = length(all.genes) - hit.pop.grey
    
    
    df$hit.pop [rownames(df) == "grey"] =  length(intersect(all.genes , target)) 
    df$fail.pop[rownames(df) == "grey"] =  length(all.genes) - df$hit.pop[rownames(df) == "grey"]      
    
    
    #Compute hypergeomertic test
    hg = apply(df,1, function(r){
      stats::phyper(q = r["hit.count"]-1  , 
                    m = r["hit.pop"]         ,
                    n = r["fail.pop"]        ,
                    k = r["drawn"]           ,
                    lower.tail = FALSE)
    })
    
    if("grey" %in% names(hg)){
      hg["grey"] = apply(df[rownames(df) == "grey",,drop=F], 1 , function(r){
        stats::phyper(q = r["hit.count"]  , 
                      m = r["hit.pop"]         ,
                      n = r["fail.pop"]        ,
                      k = r["drawn"]           ,
                      lower.tail = TRUE)
      })
    }
    
    return(list(df = df, hit.genes = hit.genes, hg = hg))
  }, simplify = F, USE.NAMES = T) 
  
  
  enrich.matrix      = lapply(enrich.df.list,function(x) x$df)
  hit.genes          = lapply(enrich.df.list,function(x) x$hit.genes)
  enrich.matrix.pval = data.frame(sapply(enrich.df.list,function(x) x$hg), check.names = F)
  
  
  enrich.matrix = sapply(1:length(enrich.matrix),hgg = enrich.matrix.pval, function(i, hgg){
    df = enrich.matrix[[i]]
    #hgg = enrich.matirx.pval
    df$hyper.geo.pval                          = df$hyper.geo.fdr = df$hyper.geo.bonf = hgg[,names(enrich.matrix)[i]]
    
    df$hyper.geo.fdr [rownames(df) != "grey"]  = p.adjust(df$hyper.geo.pval[rownames(df) != "grey"], method = "fdr")
    df$hyper.geo.bonf[rownames(df) != "grey"]  = p.adjust(df$hyper.geo.pval[rownames(df) != "grey"], method = "bonferroni") #x =(df$hyper.geo.pval)*nrow(df); x = ifelse(x>1,1,x)
    df$hyper.geo.log10.pval  = -log10(df$hyper.geo.pval)
    return(df)
  }, simplify = F, USE.NAMES = T)
  
  names(enrich.matrix) = names(enrich.df.list)
  
  
  enrich.matrix.fdr  = as.data.frame(sapply(enrich.matrix,"[[", "hyper.geo.fdr"))
  row.names(enrich.matrix.fdr) = row.names(enrich.matrix[[1]])
  
  enrich.matrix.bonf = as.data.frame(sapply(enrich.matrix,"[[", "hyper.geo.bonf"))
  row.names(enrich.matrix.bonf) = row.names(enrich.matrix[[1]])
  
  
  return(dplyr::lst(enrich.matrix, enrich.matrix.pval, enrich.matrix.bonf, enrich.matrix.fdr, ll.components, hit.genes))
  
  
}


###Begin ----
options(java.parameters = "-Xmx8000m")  #Required for xlsx package
#library(WGCNA)
library(SummarizedExperiment)
library(purrr)

set.seed(123)

library(limma)
library(dplyr)
library(pheatmap)
library(gtools)       #Optional. Just to sort path of files read
library(org.Hs.eg.db)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(directlabels)
library(patchwork)
library(tidyr)
library(cowplot)
library(ggtext)
library(SummarizedExperiment)
library(limma)
library(future)
library(purrr)
library(furrr)
library(magrittr)
library(biomaRt)
library(tidyverse)
library(qs)

our_bins          = c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")


##Read updated network names
new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\old-to-new Network names all.csv") %>% tibble::deframe()

#Read list of reference files for different pathologies [AD/ADHD/ALS/ASD/BIP/CD/MDD/MS/OCD/PD/PTSD/RA/SA/SCZ/SCZ.neg/UC][all.biotypes/protein.coding]
all.patho = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\all.patho.lists_grch38[PGC125].rds")
all.patho = all.patho[!grepl("75",names(all.patho))]
all.patho = unlist(unname(all.patho), recursive = F)

li_final = c(all.patho[c("PGC3.all.biotypes")])


##Read gene list: Sliding window (NC)
gm_AB_NC = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (sliding_window_NC)[grch38].rds")

##Read gene list: Sliding window (SchizoNew)
gm_AB_SC = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (sliding_window_SCZ)[grch38].rds")

## Read sliding window metadata (for median age for networks);
metadf = qread("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\metadf(NC+SchizoNew).qs")
metadf$new.ID = gsub("SchizoNew","SC", metadf$new.ID)

# ####################################
# #######################################
# ##All network modules  enrichment (hypergeo) for PGC3 list

#Sliding window (NC) enrichment
res_NC = imap(li_final,~{
  ref.list = .x
  print(.y)
  
  imap(gm_AB_NC,~{
    print(.y)
    tmp = Enrich(ll.target = ref.list, ll.components = .x)
  })
}) %>% transpose()

#Sliding window (Schizo) enrichment
res_SC = imap(li_final,~{
  ref.list = .x
  print(.y)
  
  imap(gm_AB_SC,~{
    print(.y)
    tmp = Enrich(ll.target = ref.list, ll.components = .x)
  })
}) %>% transpose()



#########################
#########################
###Save results on disk
#save(res_NC, res_SC, file = "C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/all_results(sliding window).RData")

##Load saved results from disk
load("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/all_results(sliding window).RData")
names(res_NC) = paste0("NC__",names(res_NC))
names(res_SC) = paste0("SC__",names(res_SC))

############
all.files = c(res_NC,res_SC) %>% imap_dfr(.id = "network",~{
  map_dfr(.id = "enrichment",.x,~{
    load_li = list(
      matrix.pval   = .x$enrich.matrix.pval,
      matrix.bonf   = .x$enrich.matrix.bonf,
      modules       = .x$ll.components,
      hits          = .x$hit.genes
    )
    
    load_li[["hits.count.df"]] = load_li$hits %>% purrr::transpose() %>% 
      map_depth(2,~length(.)) %>% map_dfr(.id = "bins", ~.) %>% data.frame(row.names = "bins", check.names = F)
    
    load_li[["hits.genes.df"]] = load_li$hits %>% imap_dfc(.id = "modules",~ {
      dd = tibble::enframe(.x)[,"value"] %>% set_names(.y)
      return(dd)
    }) %>% as.data.frame()
    
    rownames(load_li[["hits.genes.df"]]) = rownames(load_li[["hits.count.df"]])
    
    df.internal = cbind(modules = rownames(load_li$matrix.bonf), module.length = lengths(load_li$modules), bonf. = load_li$matrix.bonf, pvals. = load_li$matrix.pval, hits. = load_li$hits.count.df, hitgenes. = load_li$hits.genes.df, stringsAsFactors=F)
    
    
    fold.factor.no.grey   = colSums(df.internal[df.internal$modules != "grey", grep("hits..", colnames(df.internal))])/sum(df.internal[df.internal$modules != "grey","module.length"])
    fold.factor.all.genes = colSums(df.internal[                             , grep("hits..", colnames(df.internal))])/sum(df.internal[                             ,"module.length"])
    
    if (any(fold.factor.all.genes == 0)) {warning("check Fold Factor");}; fold.factor.all.genes[fold.factor.all.genes == 0] = 1;  #if(!fold.factor.no.grey)   fold.factor.no.grey   =1
    if (any(fold.factor.no.grey   == 0)) {warning("check Fold Factor")}; fold.factor.no.grey  [fold.factor.no.grey   == 0] = 1;  #if(!fold.factor.all.genes) fold.factor.all.genes =1
    Fold.df.no.grey   = df.internal[, grep("hits..", colnames(df.internal))] %>% sweep(1,df.internal$module.length,"/") %>% sweep(2,fold.factor.no.grey  ,"/")
    Fold.df.all.genes = df.internal[, grep("hits..", colnames(df.internal))] %>% sweep(1,df.internal$module.length,"/") %>% sweep(2,fold.factor.all.genes,"/")
    Fold.df           = Fold.df.no.grey
    Fold.df[rownames(Fold.df.no.grey) == "grey",] = Fold.df.all.genes[rownames(Fold.df.all.genes) == "grey",]
    colnames(Fold.df) = gsub("hits","Fold.hits",colnames(Fold.df))
    
    df.internal = cbind(df.internal, Fold.df)
    #df.internal = cbind(df.internal, Fold = (df.internal[,grep("hits..", colnames(df.internal))]/df.internal$module.length)/fold.factor)
    #browser()
    df.tmp = df.internal[,grep("pvals..", colnames(df.internal))]
    df.tmp[df.tmp == 1] = 1 - 1e-15
    df.tmp = apply(df.tmp,2,qnorm)
    colnames(df.tmp) = gsub("pvals","zvals", colnames(df.tmp))
    df.internal    = cbind(df.internal, df.tmp)
    df.tmp = NULL
    
    return(df.internal)
    
  })
}) %>% dplyr::select(c(network,enrichment, modules, module.length), matches(paste0("bonf..",our_bins,"$")), matches(paste0("pvals..",our_bins,"$")), matches(paste0("hitgenes..",our_bins,"$")), matches(paste0("hits..",our_bins,"$")), matches(paste0("Fold.hits..",our_bins)), contains(paste0("zvals",our_bins,"$")))

all.files_long = all.files %>% dplyr::select(-contains("hitgenes..")) %>% 
  pivot_longer(cols = -c(network,enrichment,modules,module.length),names_to = "test", values_to = "val") %>%
  drop_na(val) %>%
  separate(test   , into = c("test","bins")                                                        , sep = "\\.\\.", remove = F)   %>%
  separate(network, into = c("Dx","tissue","start_index","end_index","start_sample","end_sample")  , sep = "__"    , remove = F)   

all.files_long$median.age = metadf$median.age[match(all.files_long$network, paste0(strsplit2(metadf$new.ID,"\\.")[,1],"__",metadf$short.name))]
all.files_long$mean.age   = metadf$mean.age  [match(all.files_long$network, paste0(strsplit2(metadf$new.ID,"\\.")[,1],"__",metadf$short.name))]

###Save results on disk
#save(all.files_long, file = "C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/all.files_long(sliding window).RData")
rm(res_NC,res_SC,all.files)
gc()

##Caudate NC first network is not plotted as its median age is far from second network 
dat = all.files_long %>% 
  filter(test %in% c("Fold.hits")                          &  
           bins %in% c("PGC","kbp_200","kbp_200_negative") & 
           (!grepl("NC__caudate__1__",network))) %>%
  group_by(network,bins)                        %>% 
  mutate(var = var(val), max = max(val), bins = factor(bins,levels = c("PGC","kbp_200","kbp_200_negative"))) %>% 
    ungroup()                                   %>%
    filter((!modules %in% "grey"))

##Fold change (smooth plot)
p.smooth =  dat %>% ggplot(aes(x=median.age, y=val, group = network, color = tissue)) +
   geom_smooth(data = dat %>% filter( tissue %in% "dentate"), method = "loess",aes(group = interaction(Dx,tissue), color = tissue, fill = tissue, linetype = Dx), span= 0.65 , alpha = 0.2, size = 1.1) + 
   geom_smooth(data = dat %>% filter(!tissue %in% "dentate"), method = "loess",aes(group = interaction(Dx,tissue), color = tissue, fill = tissue, linetype = Dx), span= 0.5  , alpha = 0.2, size = 1.1) + 
  facet_grid(bins ~Dx) + 
  theme_bw(base_size = 19)  + theme(strip.text = element_text(size = 15, face = "bold")) + ylab(label = "Fold")
p.smooth

##Max Fold change (smooth plot)
p.max = dat %>% ggplot(aes(x=median.age, y=max, group = network, color = tissue)) +
  geom_smooth(data = dat %>% filter( tissue %in% "dentate"), method = "loess",aes(group = interaction(Dx,tissue), fill = tissue, color = tissue, linetype = Dx), span = 0.65, alpha = 0.2, size = 1.1) +
  geom_smooth(data = dat %>% filter(!tissue %in% "dentate"), method = "loess",aes(group = interaction(Dx,tissue), fill = tissue, color = tissue, linetype = Dx), span = 0.5 , alpha = 0.2, size = 1.1) + 
  facet_grid(bins~Dx) + 
  theme_bw(base_size = 19)  + theme(strip.text = element_text(size = 15, face = "bold")) + ylab(label = "Max")
p.max


##Variance of Fold change (smooth plot)
p.var = dat %>% ggplot(aes(x = median.age)) + 
  geom_smooth(data = dat %>% filter( tissue %in% "dentate"), method = "loess",aes(color = tissue, y = var, group = interaction(Dx,tissue), linetype = Dx), span = 0.65) +
  geom_smooth(data = dat %>% filter(!tissue %in% "dentate"), method = "loess",aes(color = tissue, y = var, group = interaction(Dx,tissue), linetype = Dx), span = 0.5 ) + 
  facet_grid(bins~Dx) + 
  theme_bw(base_size = 19)  + 
  theme(strip.text = element_text(size = 15, face = "bold"))
p.var
