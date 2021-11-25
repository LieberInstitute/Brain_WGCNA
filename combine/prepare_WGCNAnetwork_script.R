library(purrr)
library(SummarizedExperiment)
library(WGCNA)
library(limma)
library(tidyr)
library(dplyr)
library(gtools)
library(purrr)
library(dplyr)
library(magrittr)


WGCNAcolors = c("grey",WGCNA::standardColors(200)) %>% set_names(0:200)
WGCNAfn = function(.) WGCNAcolors[match(.,names(WGCNAcolors))]

neuclean = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\results\\noquantile_CellPropclean_GEclean_newcuts_noPCA--neuCells (redone)\\net\\metadf_noquantile_GEclean_newcuts_noPCA_deepSplit4_Neuclean_byAll.rds")
neuclean = neuclean[(neuclean$correlation == "pearson" & neuclean$network == "signed hybrid"),]

sva <- readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins/New_Project/newResults/sft_noAgeSqr_noPCAremoved_qSVAremoved/net/metadf_noAgeSqr_noPCAremoved_qSVAremoved_deepSplit4_byAll.rds")
sva = sva[(sva$correlation == "pearson" & sva$network == "signed hybrid"),]
sva$full.name = paste0("QSVAremoved.",sva$full.name)

nosva <- readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins/New_Project/newResults/sft_noAgeSqr_noPCAremoved_redone/net/metadf_noAgeSqr_noPCAremoved_redone_deepSplit4_byAll.rds")
nosva = nosva[(nosva$correlation == "pearson" & nosva$network == "signed hybrid"),]
nosva$full.name = paste0("noQSVAremoved.",nosva$full.name)

metadf_hippo_dentate = rbind(nosva, sva)
cc = intersect(colnames(neuclean),colnames(metadf_hippo_dentate))

metadf = rbind(neuclean[,cc], metadf_hippo_dentate[,cc])


modules_colors = metadf$test %>% set_names(strsplit2(metadf$full.name,"\\.\\.\\.")[,1])

modules_long = map_dfr(modules_colors,.id = "tissue_age",~{
   imap_dfr(.x,~ {
      data.frame(gene = .x,module = .y, stringsAsFactors = F)}
   )
})

modules_long$gene = strsplit2(modules_long$gene,"\\.")[,1]
modules_wide = tidyr::pivot_wider(modules_long,id_cols = gene, names_from = c(tissue_age), values_from = c(module), names_glue = "{tissue_age}|{.value}") %>% set_colnames(gsub("module","modules",colnames(.)))

KME = metadf$KME %>% set_names(strsplit2(metadf$full.name,"\\.\\.\\.")[,1])
KME = map(KME, ~{
   rownames(.x) = strsplit2(rownames(.x),"\\.")[,1]
   colnames(.x) = unname(WGCNAfn(substring(colnames(.x),4)))
   return(.x)
})


KME_long = map_dfr(KME,.id = "tissue_age",~ {.x %>% tibble::rownames_to_column("gene")})
KME_wide = pivot_wider(KME_long,id_cols = gene, names_from = c(tissue_age), values_from = -c(tissue_age, gene), names_glue = "{tissue_age}|{.value}")
KME_wide = KME_wide[,colSums(is.na(KME_wide)) < nrow(KME_wide)]

IMC = metadf$IMC %>% set_names(strsplit2(metadf$full.name,"\\.\\.\\.")[,1])
IMC = map(IMC, ~{
   rownames(.x) = strsplit2(rownames(.x),"\\.")[,1]
   return(.x)
})

IMC_long = map_dfr(IMC,.id = "tissue_age",~ {.x %>% tibble::rownames_to_column("gene")})
IMC_wide = tidyr::pivot_wider(IMC_long,id_cols = gene, names_from = c(tissue_age), values_from = -c(tissue_age, gene), names_glue = "{tissue_age}|{.value}")
IMC_wide = IMC_wide[,colSums(is.na(IMC_wide)) < nrow(IMC_wide)]
IMC_wide = IMC_wide[,!grepl("modules", colnames(IMC_wide))]

final_merged = merge(modules_wide, merge(IMC_wide,KME_wide, by = "gene", all = T),by = "gene", all = T) %>% tibble::column_to_rownames("gene") 
colnames(final_merged) = gsub("\\.\\.\\.", "\\.", colnames(final_merged))
colnames(final_merged) = gsub("-",".", colnames(final_merged))
colnames(final_merged) = gsub("noQSVAremoved.dentate", "dentate.noQSVAremoved",colnames(final_merged))
colnames(final_merged) = gsub("noQSVAremoved.hippo"  , "hippo.noQSVAremoved",colnames(final_merged))
colnames(final_merged) = gsub("QSVAremoved.dentate"  , "dentate.QSVAremoved",colnames(final_merged))
colnames(final_merged) = gsub("QSVAremoved.hippo"    , "hippo.QSVAremoved",colnames(final_merged))

saveRDS(final_merged, file = "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\WGCNAnetworks.rds")


# new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\old-to-new Network names all.csv") %>% tibble::deframe()
# colnames(final_merged) = paste0(new.names[strsplit2(colnames(final_merged),"\\|")[,1]],"|",strsplit2(colnames(final_merged),"\\|")[,2])
