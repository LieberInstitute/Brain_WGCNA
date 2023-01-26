library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(tidyverse)

cc <- readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/SNPextension/consensus_genes/consensus_gene_list_new.withHartl.rds")
cc = c(cc$SCZ_prenatal_positive,cc$SCZ_postnatal_positive.new)

our_bins = c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")

ne = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\PGC.gene.lists_grch38[125 PGC ensemblIDs].RData", envir = ne)
pp = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\PGC.gene.lists_grch38[120 PGC ensemblIDs].RData", envir = pp)
PGC = ne$PGC3.all.biotypes[our_bins];PGC$PGC = pp$PGC3.all.biotypes$PGC

##Function to get gene annotation from Biomart
geneMap_fun = function(genes, ensembl, attributes){
  require(biomaRt)
  ensembl_archive = useMart("ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl", host=ensembl)
  Filters = "ensembl_gene_id"
 # browser()
  df = getBM(attributes = attributes,Filters,values = genes,mart = ensembl_archive)
  df = df[df$chromosome_name %in% c(1:22,"X","Y","M"),]
  df$chromosome_name = factor(df$chromosome_name,levels = c(1:22,"X","Y","M"))
  df = df[order(df$chromosome_name,df$start_position,df$end_position),]
  return(df)
}


gene_grch38 = geneMap_fun(genes = cc, ensembl = "dec2017.archive.ensembl.org",#ensembl = "jul2016.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gene_gc_content","entrezgene","version"))

gene_grch38 = gene_grch38[!duplicated(gene_grch38$external_gene_name),]

#setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\SNPextension\\consensus_genes")
#saveRDS(gene_grch38, "consensus_gene_list_full_annotation.withHartl.rds")

bb = gene_grch38$external_gene_name
bb.in.SCZ = sort(gene_grch38$external_gene_name[gene_grch38$ensembl_gene_id %in% unlist(PGC)])

library(simplifyEnrichment)
#Load correct pre and post consensus universe
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\SNPextension\\consensus_genes\\universe.pre.RData")
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\SNPextension\\consensus_genes\\universe.post.RData")
uu0 = c(universe.pre, universe.post)

gene_grch38 = geneMap_fun(genes = uu0, ensembl = "dec2017.archive.ensembl.org",#ensembl = "jul2016.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gene_gc_content","entrezgene","version"))

gene_grch38 = gene_grch38[!duplicated(gene_grch38$external_gene_name),]

uu = gene_grch38$external_gene_name[match(uu0,gene_grch38$ensembl_gene_id)]


setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Figures")
pdf("consensus_gene_GOplots.pdf", width = 16, height = 12)
for (oo in c("BP","MF","CC")){
  go1 = enrichGO(bb, OrgDb = "org.Hs.eg.db",keyType = "SYMBOL", ont = oo, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = 1000, universe = uu)
  ##fold = ifelse(go1@gene %in% bb.in.SCZ, 1,-1) %>% set_names(go1@gene)
  go1@result = go1@result %>% filter(p.adjust <0.05)
  go1@result$geneID  = strsplit(go1@result$geneID,"/") %>% map(~ {.x[.x %in% bb.in.SCZ] = paste0("(",.x[.x %in% bb.in.SCZ],")"); paste0(.x,collapse = "/")})
  #saveRDS(go1@result,"consensus_genes_enrichGO.rds")
  
  go2 = enrichGO(bb[!bb %in% bb.in.SCZ], OrgDb = "org.Hs.eg.db",keyType = "SYMBOL", ont = oo, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = 1000, universe = uu)
  ##fold = ifelse(go1@gene %in% bb.in.SCZ, 1,-1) %>% set_names(go1@gene)
  #go2@result = go2@result %>% filter(p.adjust <0.05)
  #go2@result$geneID  = strsplit(go2@result$geneID,"/") %>% map(~ {.x[.x %in% bb.in.SCZ] = paste0("(",.x[.x %in% bb.in.SCZ],")"); paste0(.x,collapse = "/")})
  
  #print(barplot(go1, title = paste0("GO:",oo),showCategory = 10))
  #print(cnetplot(go1, showCategory = 6, title = paste0("GO:",oo), colorEdge = T, circular = T, cex_category = 1.2, cex_gene = 1.2, cex_label_category = 0.9))
  #print(cnetplot(x = go1, foldChange = bb %in% bb.in.SCZ, showCategory = 15, title = paste0("GO:",oo), colorEdge = F, circular = F, layout = "gem", cex_category = 1.2, cex_gene = 1.2, cex_label_category = 0.9, ))
  
  svg("consensus_gene_GOplots7.svg", width = 13, height = 10)
  #print(cnetplot(x = go1, showCategory = 15, title = paste0("GO:",oo), colorEdge = F, circular = F, layout = "fr", cex_category = 1, cex_gene = 0.7, cex_label_category = 0.8, cex_gene_category = 0.9, shadowtext = "gene", color_gene = "red", seed = 1000, ont = oo) + theme(plot.margin = margin(3,6,3,6)))
  cc = cnetplot(x = go1, showCategory = 15, title = paste0("GO:",oo), colorEdge = F, circular = F, layout = "dh", cex_category = 1, cex_gene = 1, cex_label_category = 1, cex_gene_category = 1, shadowtext = "gene", color_gene = "red", seed = 1000, ont = oo)
  print(cc + labs(size = "hits (genes)") + theme(legend.position = "right", legend.title=element_text(size=12), legend.text=element_text(size=12),plot.margin = margin(3,6,3,6)))
  dev.off()
 # print(simplifyEnrichment::simplifyGO(go1@result$ID))
}
dev.off()
