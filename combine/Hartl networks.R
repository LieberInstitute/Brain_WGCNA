library(tidyverse)
library(magrittr)

#rr = read.table("clipboard-128000", sep = "\t", header = T)
library(readxl)

path = "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Hartl2021\\supp1.xlsx"

net_names = excel_sheets(path)

dat_all = map_dfc(net_names,~ {          # Read all sheets to list
  dd = as.data.frame(read_excel(path, sheet = .x))[,1:2] %>% set_colnames(c("gene",.x)) %>% tibble::column_to_rownames("gene")
  }) %>% tibble::rownames_to_column("gene")

consensus_genes = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/SNPextension/consensus_genes/consensus_gene_list.rds")[c("SCZ_prenatal_positive", "SCZ_postnatal_positive")] %>% unlist
dat_all = dat_all %>% mutate(inSCZConsensus = ifelse(gene %in% consensus_genes, TRUE, FALSE), .before = gene)
write.table(dat_all,"clipboard-64000")

out = list()
for (cl in colnames(dat_all)[3:ncol(dat_all)]){
  dat = dat_all %>% select(gene,.data[[cl]]) %>% group_by(.data[[cl]]) %>% summarise(ensembl = list(gene))
  out[[cl]] = dat$ensembl %>% set_names(dat[[cl]])
}

setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis")
saveRDS(out,"Hartl2021-gene-module-list.rds")
