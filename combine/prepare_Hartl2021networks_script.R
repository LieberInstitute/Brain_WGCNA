library(tidyverse)
library(magrittr)

library(readxl)

####################
##Prepare Hartl2021 network
##https://doi.org/10.1038/s41593-021-00887-5

#Reading gene-module assignments from Supplementary Table 1
path = "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Hartl2021\\supp1.xlsx"

net_names = excel_sheets(path)

dat_all = map_dfc(net_names,~ {          # Read all sheets to list
  dd = as.data.frame(read_excel(path, sheet = .x))[,1:2] %>% set_colnames(c("gene",.x)) %>% tibble::column_to_rownames("gene")
  }) %>% tibble::rownames_to_column("gene")

out = list()
for (cl in colnames(dat_all)[3:ncol(dat_all)]){
  dat = dat_all %>% select(gene,.data[[cl]]) %>% group_by(.data[[cl]]) %>% summarise(ensembl = list(gene))
  out[[cl]] = dat$ensembl %>% set_names(dat[[cl]])
}

setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis")
saveRDS(out,"Hartl2021-gene-module-list.rds")
