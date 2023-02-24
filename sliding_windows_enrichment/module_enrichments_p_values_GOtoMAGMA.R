####################################################################################################
#           Calculate accross module association of MAGMA enrichment to GO term ratio              #
#                                                                                                  #
####################################################################################################

library(MASS);library(sfsmisc);library(dplyr);library(purrr)
library(SummarizedExperiment); library(qs);library(ggplot2); library(gridExtra);
library(pryr); library(tibble); library(furrr); library(tidyr); library(readr)
library(dplyr); library(data.table); library(grid); library(limma); 
library(pheatmap); library(RColorBrewer); library(broom); library(sfsmisc)


# gi_sw <- read_csv("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/madhur data/wide_form_test_slidingwindow_NC_SchizoNew(v1.3)_Leo_MP.csv")
# names(gi_sw)[1] = "ensemblID"
# gi_sw_networks = as.list(unique(gsub(".kTotal|.kme|.kWithin|.kout|.modules","",grep("__",names(gi_sw),value = T)))) %>% setNames(.,.) 
# 
# sampleInfo_sliding_window <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/Across age/sampleInfo_sliding_window.rds")
# sw_ages = map_dbl(sampleInfo_sliding_window[7 : length(sampleInfo_sliding_window)], ~ {
#   # print(which(..1 == 0))
#   median(sampleInfo_sliding_window$Age[which(..1 == 1)])
# })

rm(list = ls())

sw_modenr_over400_wMAGMA <- readRDS("results/sw_modenr_over400_wMAGMA.rds")

# load("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/Cell_specificity.RData")
gi <- readRDS("../madhur data/wide_form_test(v3.6)_forPaper.rds")

genes.ref = gi[c("ensemblID","gencodeID", "Symbol")]

GO_PGC3_genelists = readRDS("../data/GO_PGC3_genelists")
GO_PGC3_allgenes = unique(unlist(GO_PGC3_genelists))
# GO_terms = grep("central_nervous_system", names(GO_PGC3_genelists), value = T)
GO_terms = names(GO_PGC3_genelists)[which(lapply(GO_PGC3_genelists, length) >400)]
GO_terms = grep("synapse|somato|neuron|nervous|chromatin|development|DNA|RNA|_ion_",GO_terms, value = T)

SCZ_stat = "SCZ_MAGMA"
GO_terms_list = as.list(GO_terms) %>% setNames(.,.)

pb <- progress::progress_bar$new(total = length(GO_terms_list), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
sw_corr_list = map(GO_terms_list, ~ {
  goterm = ..1
  pb$tick()
  # Sys.sleep(.001)
  gi_sw_corr = imap(sw_modenr_over400_wMAGMA, ~ {
    mod_enr_nogrey = ..1[!rownames(..1) == "grey",]
    enr_stat = mod_enr_nogrey[[goterm]]
    rsl <- rlm(-log10(mod_enr_nogrey[[SCZ_stat]]) ~ mod_enr_nogrey[["mod_size"]] + enr_stat )
    pvalue = f.robftest(rsl, var = "enr_stat")
    list(rsl = rsl,
         pvalue = pvalue)
  })
})

pb <- progress::progress_bar$new(total = length(sw_corr_list), format = " [:bar] :current/:total (:percent) eta: :eta", force = TRUE)
sw_rlm_stats = map(sw_corr_list, ~ {
  GO_rlms = ..1
  pb$tick()
  out = imap_dfr(GO_rlms, ~ {
    network_rlm = ..1
    c(
      Z_stat = broom::tidy(..1$rsl) %>% .[grep("enr_stat",.$term),] %>% pull(statistic),
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

broom::tidy(sw_corr_list[["GO:0006325_chromatin_organization"]][["NC.caudate__1"]]$rsl) %>% .[grep("mod_enr_nogrey\\[\\[goterm\\]\\]",.$term),] %>% pull(statistic)

saveRDS(sw_rlm_stats,"results/sw_modenrMAGMA_GOterms_CNS_rlm_stats.rds")

