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

DRONC_human  = as.data.frame(Cell_specificity[["human"]][["DRONC_human"]])
DRONC_human  = tibble::rownames_to_column(DRONC_human, "Symbol")
DRONC_human  = merge(DRONC_human, genes.ref)    
DRONC_human <- DRONC_human[, 
                           c("ensemblID", "gencodeID", "Symbol" ,  "exCA",  "exDG", "exPFC",  "GABA", "NSC", "OPC", "ODC", "ASC", "MG" , "END" )]
DRONC_types   = colnames(DRONC_human[,4:(length(DRONC_human))])
# gi$CN_Juvenile.modules

GO_PGC3_genelists = readRDS("../data/GO_PGC3_genelists")
GO_PGC3_allgenes = unique(unlist(GO_PGC3_genelists))
# GO_terms = grep("central_nervous_system", names(GO_PGC3_genelists), value = T)
GO_terms = names(GO_PGC3_genelists)[which(lapply(GO_PGC3_genelists, length) >400)]


gi_sw_HMAGMA = readRDS(gi_sw_HMAGMA, file = "../data/gi_sw_HMAGMA.rds")

net = "NC.caudate__129"

get_mod_genes = function(gi, net){
  mod_list = as.list(unique(gi[[paste0(net, ".modules")]])) %>% setNames(.,.)
  mod_genes = map(mod_list, ~ {
    gi$ensemblID[gi[[paste0(net, ".modules")]] == ..1 & !is.na(gi[[paste0(net, ".modules")]])]
  })
  return(mod_genes)
}
mod_genes = get_mod_genes(gi_sw,net)

get_mod_enr_GO= function(mod_genes, net_name,  with_grey = T){
  mod_genes_nogrey = mod_genes[!names(mod_genes) == "grey" & !is.na(names(mod_genes))]
  mod_genes = mod_genes[!is.na(names(mod_genes))]
  if (with_grey == T){
    mod_genes_df = mod_genes
  } else {mod_genes_df = mod_genes_nogrey}
  mod_enrichments= imap_dfr(mod_genes_df, ~ {
    genes = ..1
    mod_name = ..2
    c(
      map_dfr(as.list(GO_terms) %>% setNames(.,.), ~ {
        GO_term = ..1
        fold_change = length(which(genes %in% GO_PGC3_genelists[[GO_term[1]]]))/length(genes)
      }),
      # SCZ_HMAGMA = geneSetTest(index = which(gi$ensemblID %in% genes ), statistics = gi[, "SCZ.PGC3.Adult_brain.ZSTAT"],
      #                         alternative = "up", type= "t", ranks.only = T),
      
      SCZ_HMAGMA_median = median(gi[which(gi$ensemblID %in% genes ), "SCZ.PGC3.Adult_brain.ZSTAT"]), 
      mod_size= length(genes), 
      mod_name= mod_name)
  }) %>% column_to_rownames( var = "mod_name")  %>% bind_cols(gi_sw_HMAGMA[[net_name]])
  
  return(mod_enrichments)
}

pb <- progress::progress_bar$new(total = length(gi_sw_networks), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
gi_sw_modenr = imap(gi_sw_networks, ~ {
  pb$tick()
  # Sys.sleep(.001)
  mod_genes = get_mod_genes(gi_sw,..2)
  mod_enr = get_mod_enr_GO(mod_genes, ..2,with_grey=T)
})

saveRDS(gi_sw_modenr, file = "results/sw_modenr_over400")

goterm = "GO:0045202_synapse"
goterm = "GO:0007155_cell_adhesion"
goterm = "GO:0043005_neuron_projection"  
goterm = "GO:0008324_cation_transmembrane_transporter_activity" 
goterm = "GO:0048468_cell_development" 

GO_terms_list = as.list(GO_terms) %>% setNames(.,.)
pb <- progress::progress_bar$new(total = length(GO_terms_list), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
sw_modenrHMAGMA_corr_list = map(GO_terms_list, ~ {
  goterm = ..1
  pb$tick()
  # Sys.sleep(.001)
  gi_sw_modenrHMAGMA_corr = imap_dbl(gi_sw_modenr, ~ {
    mod_enr_nogrey = ..1[!rownames(..1) == "grey",]
    broom::tidy(rlm(-log10(mod_enr_nogrey[["SCZ_HMAGMA"]]) ~ mod_enr_nogrey[["mod_size"]] + mod_enr_nogrey[[goterm]] )) %>% .[grep("goterm",.$term),] %>% pull(statistic)
    # cell_type = lm(-log10(mod_enr_nogrey$`GO:0007417_central_nervous_system_development`) ~ mod_enr_nogrey$mod_size)$residuals
    # SCZ_HMAGMA = lm(-log10(mod_enr_nogrey$SCZ_HMAGMA) ~ mod_enr_nogrey$mod_size)$residuals
    # cor.test(cell_type,SCZ_HMAGMA)
  })
})

sw_modenrHMAGMA_corr_list = readRDS( file = "results/sw_modenrHMAGMA_corr_list.rds")
sw_modenrHMAGMA_corr_list[["median_age"]] = sw_ages

sw_modenrHMAGMA_corr_wide = as.data.frame(sw_modenrHMAGMA_corr_list) %>% rownames_to_column( var = "network") 
sw_modenrHMAGMA_corr_long = data.table::melt(setDT(sw_modenrHMAGMA_corr_wide), 
                                 id.vars = c("network","median_age"), 
                                 variable.name = "GO_term")

sw_modenrHMAGMA_corr = save(sw_modenrHMAGMA_corr_list,sw_modenrHMAGMA_corr_wide,sw_modenrHMAGMA_corr_long
                            ,file = "results/sw_modenrHMAGMA_corr.RData")

gi_sw_modenrHMAGMA_corr = imap_dbl(gi_sw_modenr, ~ {
  mod_enr_nogrey = ..1[!rownames(..1) == "grey",]
  broom::tidy(rlm(-log10(mod_enr_nogrey[["SCZ_HMAGMA"]]) ~ mod_enr_nogrey[["mod_size"]] + mod_enr_nogrey[[goterm]] )) %>% .[grep("goterm",.$term),] %>% pull(statistic)
  # cell_type = lm(-log10(mod_enr_nogrey$`GO:0007417_central_nervous_system_development`) ~ mod_enr_nogrey$mod_size)$residuals
  # SCZ_HMAGMA = lm(-log10(mod_enr_nogrey$SCZ_HMAGMA) ~ mod_enr_nogrey$mod_size)$residuals
  # cor.test(cell_type,SCZ_HMAGMA)
})
tissue = "NC.hippo"
indices = grep(tissue, names(gi_sw_networks)) 
plot(sw_ages[indices],gi_sw_modenrHMAGMA_corr[indices])

mod_enr_nogrey[mod_enr_nogrey == 0] = min(mod_enr_nogrey[!mod_enr_nogrey == 0])

mod_enr_nogrey = gi_sw_modenr[[50]]
cell_type = mod_enr_nogrey$`GO:0007417_central_nervous_system_development`
SCZ_HMAGMA = lm(-log10(mod_enr_nogrey$SCZ_HMAGMA) ~ mod_enr_nogrey$mod_size)$residuals
cor.test(cell_type,SCZ_HMAGMA, method = "kendall")
plot(cell_type, SCZ_HMAGMA)
