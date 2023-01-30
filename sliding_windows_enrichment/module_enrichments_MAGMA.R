library(dplyr)
library(purrr)
library(SummarizedExperiment); library(qs);library(ggplot2); library(gridExtra);library(pryr); library(tibble); library(furrr); library(tidyr); library(readr)
library(dplyr); library(data.table); library(grid); library(limma); library(pheatmap); library(RColorBrewer); library(broom); library(MASS)
rm(list = ls())
load("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/Cell_specificity.RData")
gi <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/madhur data/wide_form_test(v3.6)_forPaper.rds")
gi_sw <- read_csv("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/madhur data/wide_form_test_slidingwindow_NC_SchizoNew(v1.3)_Leo_MP.csv")
names(gi_sw)[1] = "ensemblID"
gi_sw_networks = as.list(unique(gsub(".kTotal|.kme|.kWithin|.kout|.modules","",grep("__",names(gi_sw),value = T)))) %>% setNames(.,.) 


sampleInfo_sliding_window <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/Across age/sampleInfo_sliding_window.rds")
sw_ages = map_dbl(sampleInfo_sliding_window[7 : length(sampleInfo_sliding_window)], ~ {
  # print(which(..1 == 0))
  median(sampleInfo_sliding_window$Age[which(..1 == 1)])
})


genes.ref = gi[c("ensemblID","gencodeID", "Symbol")]

get_mod_genes = function(gi, net){
  mod_list = as.list(unique(gi[[paste0(net, ".modules")]])) %>% setNames(.,.)
  mod_genes = map(mod_list, ~ {
    gi$ensemblID[gi[[paste0(net, ".modules")]] == ..1 & !is.na(gi[[paste0(net, ".modules")]])]
  })
  return(mod_genes)
}
# mod_genes = get_mod_genes(gi_sw,net)
# gi$SCZ.PGC3.kbp35.10.ZSTAT

get_mod_enrichments = function(mod_genes, with_grey = T){
  mod_genes_nogrey = mod_genes[!names(mod_genes) == "grey" & !is.na(names(mod_genes))]
  mod_genes = mod_genes[!is.na(names(mod_genes))]
  if (with_grey == T){
    mod_genes_df = mod_genes
  } else {mod_genes_df = mod_genes_nogrey}
  
  mod_enrichments= imap_dfr(mod_genes_df, ~ {
    genes = ..1
    mod_name = ..2
    
    SCZ_MAGMA = gi[which(gi$ensemblID %in% genes ), "SCZ.PGC3.kbp35.10.ZSTAT"]
    SCZ_MAGMA = SCZ_MAGMA[which(!is.na(SCZ_MAGMA))]
    c(
    # map_dfr(as.list(DRONC_types) %>% setNames(.,.), ~ {
      # cell_type = ..1
    #   DRONC_score = DRONC_human[which(DRONC_human$ensemblID %in% genes ), cell_type]
    #   DRONC_score = DRONC_score[which(!is.na(DRONC_score))]
    #   DRONC_score_median =   median(DRONC_score)
    #   DRONC_score_max =   max(DRONC_score)
    #   DRONC_score_mean =   mean(DRONC_score)
    # }),
    SCZ_MAGMA =geneSetTest(index = which(gi$ensemblID %in% genes ), statistics = gi[, "SCZ.PGC3.kbp35.10.ZSTAT"],
                  alternative = "up", type= "t", ranks.only = T),
    SCZ_MAGMA_median =   median(SCZ_MAGMA),
    SCZ_MAGMA_max =   max(SCZ_MAGMA),  
    SCZ_MAGMA_mean =   mean(SCZ_MAGMA),  
    mod_size=length(genes), 
    mod_name=mod_name)
  }) %>% column_to_rownames( var = "mod_name")
  
  return(mod_enrichments)
}

pb <- progress::progress_bar$new(total = length(gi_sw_networks), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
gi_sw_modenr = imap(gi_sw_networks, ~ {
  pb$tick()
  # Sys.sleep(.001)
  mod_genes = get_mod_genes(gi_sw,..2)
  mod_enr = get_mod_enrichments(mod_genes, with_grey=T)
})

saveRDS(gi_sw_modenr,"results/sw_MAGMA.rds")

save(sw_modenrHMAGMA_corr_list,sw_modenrHMAGMA_corr_wide,sw_modenrHMAGMA_corr_long
     ,file = "results/sw_modenrHMAGMA_DRONCmedians_corr.RData")

