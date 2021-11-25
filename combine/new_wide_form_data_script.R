library(purrr)
#library(SummarizedExperiment)
#library(WGCNA)
library(limma)
library(tidyr)
library(dplyr)
library(gtools)
library(purrr)
library(dplyr)
library(magrittr)

load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Network Group\\In progress\\SCZ_gene_importance_network_v3.5.RData")

colnames(WGCNA) = gsub("silver","grey60",colnames(WGCNA))
colnames(WGCNA) = gsub("\\.| ","_",colnames(WGCNA))
all.colors = paste0(c("grey",WGCNA::standardColors(150)), collapse = "|")
colnames(WGCNA)[grep(all.colors, colnames(WGCNA), value = F)] =  gsub("\\|","|KME",grep(all.colors, colnames(WGCNA), value = T))
colnames(WGCNA) = gsub("\\|",".",colnames(WGCNA))


networks = c("CN", "CN_Juvenile", "CN_Adult", "CN_Older_Adult",
             "DLPFC", "DLPFC_Perinatal", "DLPFC_Juvenile", "DLPFC_Adult", "DLPFC_Older_Adult",
             "HP", "HP_Perinatal", "HP_Juvenile", "HP_Adult", "HP_Older_Adult",
             "DG", "DG_noQSVA", "HP_noQSVA", "DG_QSVA", "HP_QSVA",
             "Fromer2016_case", "Fromer2016_control",
             "Gandal2018", "Gandal2018PE", "Gandal2018PE_cs",
             "Li2018",
             "Pergola2017", "Pergola2019", "Pergola2020",
             "Radulescu2020", "Walker2019", "Werling2020")

coltypes = c("modules","kDiff","kOut","kTotal","kWithin","KME")

grid = expand.grid(networks = networks, coltypes = coltypes) %>% mutate(pattern = paste0(networks,".",coltypes))

WGCNA1 = map_dfc(grid$pattern, ~ {
  cols = grep(.x,colnames(WGCNA),value = T)
  print(cols)
  WGCNA[,cols,drop=F]
})


metadata = rbind(metadata_grch38,
                metadata_grch37[metadata_grch37$ensemblID %in% setdiff(metadata_grch37$ensemblID, metadata_grch38$ensemblID), ]) %>% set_rownames(.$ensemblID)
metadata$ensemblID = NULL

li = dplyr::lst(metadata,
          metadata.MAGMA,
          SDA,
          WGCNA1)


dd = li %>% map(~ as.data.frame(.x) %>% tibble::rownames_to_column("ensemblID")) %>%
  purrr::reduce(dplyr::full_join, by = "ensemblID") %>%
  set_rownames(.$ensemblID)


#dd1 = dd[,c(1:29,44,69,94,119,3359:4596)]
dd1 = dd1[,!grepl("kDiff", colnames(dd1))]
saveRDS(dd1,"wide_form_test(v3.6)_forPaper.rds")
