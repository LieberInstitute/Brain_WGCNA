####################################################################################################
#                Script to Get a gene_module_list for all networks from the wide_data_file         #
#                                to be used for enrichment analysis                                #
#           "Regional-coexpression", "Age-period" and "Cell-population enrichment" studies         #
#                                                                                                  #
####################################################################################################

library(dplyr)
library(tidyr)
library(MASS)
library(psych)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(limma)
library(future)
library(purrr)
library(furrr)
library(magrittr)
library(biomaRt)
library(tidyverse)

set.seed(123)

geneMap_fun = function(genes, ensembl = "grch37.ensembl.org", attributes){
  require(biomaRt)
  ensembl_archive = useMart("ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl", host=ensembl)
  Filters = "ensembl_gene_id"
  
  df = getBM(attributes = attributes,Filters,values = genes,mart = ensembl_archive)
  df = df[df$chromosome_name %in% c(1:22,"X","Y","M"),]
  df$chromosome_name = factor(df$chromosome_name,levels = c(1:22,"X","Y","M"))
  df = df[order(df$chromosome_name,df$start_position,df$end_position),]
  return(df)
}


df = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Network Group\\wide_form_test(v3.5).rds")
df = df[,grep("modules", colnames(df))] #%>% tibble::rownames_to_column("ENSEMBL")


big_list = df %>% map(~ {
  xx = data.frame(ENSEMBL = rownames(df), modules = ., stringsAsFactors = F) %>% 
    group_by(modules) %>% summarise(vals = list(ENSEMBL)) %>% filter(!is.na(modules))
  return(xx$vals %>% set_names(xx$modules))
}) %>% `names<-`(gsub("\\|modules","",names(.)))
  
#all.genes = unique(c(unlist(big_list),unlist(big_list_NC),unlist(big_list_SCZ)))
all.genes = unique(unlist(big_list))

gene_grch38 = geneMap_fun(genes = all.genes, ensembl = "jul2016.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gc_content","entrezgene","version"))
gene_grch38$strand = ifelse(gene_grch38$strand > 0, "+","-")

all.genes_grch38 = unique(gene_grch38$ensembl_gene_id)
big_list_grch38 = map(big_list,~ {
  map(.x, ~{intersect(.x,all.genes_grch38)})
})

saveRDS(big_list, "gene-module list (wide_form_test) (all networks)[grch38].rds")
