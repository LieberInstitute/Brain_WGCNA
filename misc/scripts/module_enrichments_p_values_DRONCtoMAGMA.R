library(MASS);library(sfsmisc);library(dplyr);library(purrr)
library(SummarizedExperiment); library(qs);library(ggplot2); library(gridExtra);
library(pryr); library(tibble); library(furrr); library(tidyr); library(readr)
library(dplyr); library(data.table); library(grid); library(limma); 
library(pheatmap); library(RColorBrewer); library(broom)

rm(list = ls())

sw_genesettest_celltype_enr <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/wgcna age/results/sw_genesettest_celltype_wMAGMA.rds")

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
DRONC_types   = colnames(DRONC_human[,4:(length(DRONC_human))])
DRONC_types_list = as.list(DRONC_types) %>% setNames(.,.)

pb <- progress::progress_bar$new(total = length(DRONC_types_list), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
sw_modenrMAGMA_corr_list = map(DRONC_types_list, ~ {
  cell_type = ..1
  pb$tick()
  # Sys.sleep(.001)
  gi_sw_modenrHMAGMA_corr = imap(sw_genesettest_celltype_enr, ~ {
    mod_enr_nogrey = ..1[!rownames(..1) == "grey",]
    log10celltype = -log10(mod_enr_nogrey[[cell_type]])
    # print(log10celltype)
    log10celltype[log10celltype==Inf] = max(log10celltype[!log10celltype==Inf])
    rsl <- rlm(-log10(mod_enr_nogrey[["SCZ_MAGMA"]]) ~ mod_enr_nogrey[["mod_size"]] + log10celltype )
    # Z_stat = broom::tidy(rsl) %>% .[grep("log10celltype",.$term),] %>% pull(statistic)
    pvalue = f.robftest(rsl, var = "log10celltype")
    list(rsl = rsl,
      pvalue = pvalue)
  })
})

sw_modenrMAGMA_rlm_stats = map(sw_modenrMAGMA_corr_list, ~ {
  celltype_rlms = ..1
  out = imap_dfr(celltype_rlms, ~ {
    network_rlm = ..1
    c(
      Z_stat = broom::tidy(..1$rsl) %>% .[grep("log10celltype",.$term),] %>% pull(statistic),
      f.robftest_Fstat = network_rlm[["pvalue"]][["statistic"]][["F"]],
      f.robftest_pvalue = network_rlm[["pvalue"]][["p.value"]],
      f.robftest_log10pvalue = -log10(network_rlm[["pvalue"]][["p.value"]]),
      network = ..2
    )
  }) %>% column_to_rownames("network")
  out$f.robftest_pvalue_fdr = p.adjust(out$f.robftest_pvalue, method = "fdr" )
  out$f.robftest_log10pvalue_fdr = -log10(out$f.robftest_pvalue_fdr)
  out$f.robftest_pvalue_bonf = p.adjust(out$f.robftest_pvalue, method = "bonferroni" )
  out$f.robftest_log10pvalue_bonf = -log10(out$f.robftest_pvalue_bonf)
  out
})

saveRDS(sw_modenrMAGMA_rlm_stats,"results/sw_modenrMAGMA_rlm_stats.rds")


