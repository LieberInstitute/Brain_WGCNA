####################################################################################################
#             Add MAGMA module enrichments to GO and DRONC celltype enrichment dataframes          #
#                                                                                                  #
####################################################################################################

library(dplyr)
library(purrr)
library(SummarizedExperiment); library(qs);library(ggplot2); library(gridExtra);library(pryr); library(tibble); library(furrr); library(tidyr); library(readr)
library(dplyr); library(data.table); library(grid); library(limma); library(pheatmap); library(RColorBrewer); library(broom); library(MASS)

sw_MAGMA <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/wgcna age/results/sw_MAGMA.rds")
sw_modenr_over400 <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/wgcna age/results/sw_modenr_over400.rds")
sw_genesettest_celltype_enr <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/wgcna age/results/sw_genesettest_celltype_enr.rds")

sw_modenr_over400_wMAGMA = map2(sw_modenr_over400,sw_MAGMA, ~ {
  if(all(rownames(..1)==rownames(..2))) {
    MAGMA_df = sapply(..2, as.numeric,by="row.names")
    cbind(..1,MAGMA_df)
  } else {break}
})

sw_genesettest_celltype_wMAGMA = map2(sw_genesettest_celltype_enr,sw_MAGMA, ~ {
  if(all(rownames(..1)==rownames(..2))) {
    MAGMA_df = sapply(..2, as.numeric,by="row.names")
    cbind(..1,MAGMA_df)
  } else {break}
})

saveRDS(sw_modenr_over400_wMAGMA,"results/sw_modenr_over400_wMAGMA.rds")
saveRDS(sw_genesettest_celltype_wMAGMA,"results/sw_genesettest_celltype_wMAGMA.rds")
