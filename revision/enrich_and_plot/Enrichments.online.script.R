####################################################################################################
#                    Script to run various enrichments on gene_module_list:                        #
#                   SCZ, DEGs, DMGs, TWAS, GO, MAGMA, Cell Specificity etc.                        #
#                                                                                                  #
#                  Script to make sankey plots for DLPFC, HP and CN SCZ risk genes.                #
#                                                                                                  #
#                                                                                                  #
#                          Prepare and format data for other visualisations                        #
#           "Regional-coexpression", "Age-period" and "Cell-population enrichment" studies         #
#                                                                                                  #
####################################################################################################


##Function to get gene annotation from Biomart
geneMap_fun = function(genes, ensembl, attributes){
  require(biomaRt)
  ensembl_archive = useMart("ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl", host=ensembl)
  Filters = "ensembl_gene_id"   #"hgnc_symbol"
  
  df = getBM(attributes = attributes,Filters,values = genes,mart = ensembl_archive)
  df = df[df$chromosome_name %in% c(1:22,"X","Y","M"),]
  df$chromosome_name = factor(df$chromosome_name,levels = c(1:22,"X","Y","M"))
  df = df[order(df$chromosome_name,df$start_position,df$end_position),]
  return(df)
}

#Function to perform hypergeometric enrichment on genesets (over-representation for non-grey modules, under-representation for grey modules)
Enrich = function(ll.target, ll.components){
  
  all.genes = unique(unlist(ll.components)) 
  all.genes.except.greys = unique(unlist(ll.components[!(names(ll.components) %in% "grey")]))   #All genes except Grey genes (already exlcuded)
  
  
  ## create list of dataframes for hypergeometric test
  enrich.df.list = sapply(ll.target, ll.components = ll.components, function(target, ll.components) {
    
    hit.genes     = sapply(ll.components, intersect, target)
    hit.count     = sapply(hit.genes    , length)
    hit.pop       = length(intersect(all.genes.except.greys, target))
    fail.pop      = length(all.genes.except.greys) - hit.pop
    drawn         = sapply(ll.components,length)
    df            = data.frame(hit.count = hit.count ,
                               hit.pop    = hit.pop   ,
                               fail.pop   = fail.pop  ,
                               drawn      = drawn     )
    
    #   hit.pop.grey  = length(intersect(all.genes , target))
    #   fail.pop.grey = length(all.genes) - hit.pop.grey
    
    
    df$hit.pop [rownames(df) == "grey"] =  length(intersect(all.genes , target)) 
    df$fail.pop[rownames(df) == "grey"] =  length(all.genes) - df$hit.pop[rownames(df) == "grey"]      
    
    
    #Compute hypergeomertic test
    hg = apply(df,1, function(r){
      stats::phyper(q = r["hit.count"]-1  , 
                    m = r["hit.pop"]         ,
                    n = r["fail.pop"]        ,
                    k = r["drawn"]           ,
                    lower.tail = FALSE)
    })
    
    if("grey" %in% names(hg)){
      hg["grey"] = apply(df[rownames(df) == "grey",,drop=F], 1 , function(r){
        stats::phyper(q = r["hit.count"]  , 
                      m = r["hit.pop"]         ,
                      n = r["fail.pop"]        ,
                      k = r["drawn"]           ,
                      lower.tail = TRUE)
      })
    }
    
    return(list(df = df, hit.genes = hit.genes, hg = hg))
  }, simplify = F, USE.NAMES = T) 
  
  
  enrich.matrix      = lapply(enrich.df.list,function(x) x$df)
  hit.genes          = lapply(enrich.df.list,function(x) x$hit.genes)
  enrich.matrix.pval = data.frame(sapply(enrich.df.list,function(x) x$hg), check.names = F)
  
  
  enrich.matrix = sapply(1:length(enrich.matrix),hgg = enrich.matrix.pval, function(i, hgg){
    df = enrich.matrix[[i]]
    #hgg = enrich.matirx.pval
    df$hyper.geo.pval                          = df$hyper.geo.fdr = df$hyper.geo.bonf = hgg[,names(enrich.matrix)[i]]
    
    df$hyper.geo.fdr [rownames(df) != "grey"]  = p.adjust(df$hyper.geo.pval[rownames(df) != "grey"], method = "fdr")
    df$hyper.geo.bonf[rownames(df) != "grey"]  = p.adjust(df$hyper.geo.pval[rownames(df) != "grey"], method = "bonferroni") #x =(df$hyper.geo.pval)*nrow(df); x = ifelse(x>1,1,x)
    df$hyper.geo.log10.pval  = -log10(df$hyper.geo.pval)
    return(df)
  }, simplify = F, USE.NAMES = T)
  
  names(enrich.matrix) = names(enrich.df.list)
  
  
  enrich.matrix.fdr  = as.data.frame(sapply(enrich.matrix,"[[", "hyper.geo.fdr"))
  row.names(enrich.matrix.fdr) = row.names(enrich.matrix[[1]])
  
  enrich.matrix.bonf = as.data.frame(sapply(enrich.matrix,"[[", "hyper.geo.bonf"))
  row.names(enrich.matrix.bonf) = row.names(enrich.matrix[[1]])
  
  
  return(dplyr::lst(enrich.matrix, enrich.matrix.pval, enrich.matrix.bonf, enrich.matrix.fdr, ll.components, hit.genes))
  
  
}


#Function to perform hypergeometric TWAS enrichment on genesets
Enrich.TWAS = function(ll.target, ll.components, target.name, ll.universe, gene_map = gene_grch38){
  all.modules.names = names(ll.components)
  
  ## create list of dataframes for hypergeometric test
  enrich.df.list = sapply(names(ll.target), ll.components = ll.components, function(target.name, ll.components, universe = ll.universe) {
    
    #Use different universe based on the tissue type
    if (target.name %in% c("Gandal","Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Jaffe_HIPPO_pgc2.clozuk", "Jaffe_DLPFC_pgc2", "Jaffe_HIPPO_pgc2", "Hall", "Paquola_CAUDATE")){
      universe.name = dplyr::case_when(
        target.name == "Gandal"  ~ "Gandal" , #14750
        target.name == "Gusev"   ~ "Gusev" ,  
        #target.name == "Huckins" ~ "Huckins",#10930
        target.name == "Jaffe_DLPFC_pgc2.clozuk" ~ "Jaffe_DLPFC_pgc2.clozuk",
        target.name == "Jaffe_HIPPO_pgc2.clozuk" ~ "Jaffe_HIPPO_pgc2.clozuk",
        target.name == "Jaffe_DLPFC_pgc2"        ~ "Jaffe_DLPFC_pgc2", 
        target.name == "Jaffe_HIPPO_pgc2"        ~ "Jaffe_HIPPO_pgc2",  
        target.name == "Hall"                    ~ "Hall" ,              
        target.name == "Paquola_CAUDATE"         ~ "Paquola_CAUDATE" )
      
      target = ll.target[[target.name]]
      universe = ll.universe[[universe.name]]
      ll.components = sapply(ll.components, intersect, universe)
      ll.components = ll.components[lengths(ll.components) > 10]
      hit.genes     = sapply(ll.components, intersect, target)
      hit.count     = sapply(hit.genes    , length)
      hit.pop       = length(intersect(universe, target))
      fail.pop      = length(universe) - hit.pop
      drawn         = sapply(ll.components, function(x) {length(intersect(x, universe))})
      df            = data.frame(hit.count = hit.count ,
                                 hit.pop    = hit.pop   ,
                                 fail.pop   = fail.pop  ,
                                 drawn      = drawn     )
      
      df$hit.pop [rownames(df) == "grey"] =  length(intersect(universe , target)) 
      df$fail.pop[rownames(df) == "grey"] =  length(universe) - df$hit.pop[rownames(df) == "grey"]      
      
    }
    
    
    if (target.name %in% c("Huckins")){
      target = ll.target[[target.name]]
      g = gene_map$ensemblID[gene_map$gene_type == "protein_coding"]
      ll.components = sapply(ll.components, function(x) {
        return(intersect(x, g))})
      universe = 10930
      
      
      
      ll.components = ll.components[(lengths(ll.components) < universe) & (lengths(ll.components) > 10)]
      
      hit.genes     = sapply(ll.components, intersect, target)
      hit.count     = sapply(hit.genes    , length)
      hit.pop       = length(target)
      fail.pop      = universe - hit.pop
      drawn         = sapply(ll.components, function(x) {min(universe,length(x))})#sapply(ll.components, length)
      df            = data.frame(hit.count = hit.count ,
                                 hit.pop    = hit.pop   ,
                                 fail.pop   = fail.pop  ,
                                 drawn      = drawn     )
      
      
      df$hit.pop [rownames(df) == "grey"] =  length(target)
      df$fail.pop[rownames(df) == "grey"] =  universe - df$hit.pop[rownames(df) == "grey"]  
      if (any(sapply(ll.components,length) > universe)) {
        df$hit.pop[rownames(df) %in% names(ll.components)[sapply(ll.components,length) > universe]] = NA
        df$fail.pop[rownames(df) %in% names(ll.components)[sapply(ll.components, length) > universe]] = NA
      }
      
    }
    
    removed.modules = all.modules.names[!(all.modules.names %in% rownames(df))]
    if (length(removed.modules) > 0){
      df[removed.modules,] = NA
      df = df[match(all.modules.names, rownames(df)),]
      
      hit.genes[removed.modules] = NA
      hit.genes = hit.genes[match(all.modules.names, names(hit.genes))]
      
    }
    
    
    #Compute hypergeomertic test
    hg = apply(df,1, function(r){
      stats::phyper(q = r["hit.count"]-1  , 
                    m = r["hit.pop"]         ,
                    n = r["fail.pop"]        ,
                    k = r["drawn"]           ,
                    lower.tail = FALSE)
    })
    
    if("grey" %in% names(hg)){
      hg["grey"] = apply(df[rownames(df) == "grey",,drop=F], 1 , function(r){
        stats::phyper(q = r["hit.count"]  , 
                      m = r["hit.pop"]         ,
                      n = r["fail.pop"]        ,
                      k = r["drawn"]           ,
                      lower.tail = TRUE)
      })
    }
    return(list(df = df, hit.genes = hit.genes, hg = hg))
  }, simplify = F, USE.NAMES = T) 
  
  xx = sapply(enrich.df.list,function(x) x$hg, simplify = F)
  
  enrich.matrix      = lapply(enrich.df.list,function(x) x$df)
  hit.genes          = lapply(enrich.df.list,function(x) x$hit.genes)
  enrich.matrix.pval = data.frame(sapply(enrich.df.list,function(x) x$hg), check.names = F)
  
  enrich.matrix = sapply(1:length(enrich.matrix),hgg = enrich.matrix.pval, function(i, hgg){
    df = enrich.matrix[[i]]
    df$hyper.geo.pval                          = df$hyper.geo.fdr = df$hyper.geo.bonf = hgg[,names(enrich.matrix)[i]]
    df$hyper.geo.fdr [rownames(df) != "grey"]  = p.adjust(df$hyper.geo.pval[rownames(df) != "grey"], method = "fdr")
    df$hyper.geo.bonf[rownames(df) != "grey"]  = p.adjust(df$hyper.geo.pval[rownames(df) != "grey"], method = "bonferroni") #x =(df$hyper.geo.pval)*nrow(df); x = ifelse(x>1,1,x)
    df$hyper.geo.log10.pval  = -log10(df$hyper.geo.pval)
    
    return(df)
  }, simplify = F, USE.NAMES = T)
  
  names(enrich.matrix) = names(enrich.df.list)
  
  yy = map(enrich.matrix,~ .x$hyper.geo.fdr %>% set_names(rownames(.x)))
  enrich.matrix.fdr = map_dfr(.id = "list",yy,~.x) %>% tibble::column_to_rownames("list") %>% t
  
  yy = map(enrich.matrix,~ .x$hyper.geo.bonf %>% set_names(rownames(.x)))
  enrich.matrix.bonf = map_dfr(.id = "list",yy,~.x) %>% tibble::column_to_rownames("list") %>% t
  
  return(dplyr::lst(enrich.matrix, enrich.matrix.pval, enrich.matrix.bonf, enrich.matrix.fdr, ll.components, hit.genes))
  
}


#Function to get GO/KEGG/Reactome enrichment
#### GO enrichment ####
GO_enrich = function(Net,method.adjust = "bonf"){
  
  library(clusterProfiler, attach.required = T)
  library(org.Hs.eg.db)
  library(DOSE, attach.required = T)
  
  library(magrittr)
  library(purrr)
  library(limma)
  
  ll.modules = Net[!names(Net) %in% "grey"]
  universe   = as.character(unlist(ll.modules))
  
  ll.modules.ENTREZ = map(ll.modules, ~ {clusterProfiler::bitr(., fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID})
  universe.ENTREZ   = as.character(unlist(ll.modules.ENTREZ))
  
  
  ############ Main -------------
  result_all_modules_BP      = compareCluster(geneClusters = ll.modules, fun = "enrichGO_custom", universe = universe, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  result_all_modules_MF      = compareCluster(geneClusters = ll.modules, fun = "enrichGO_custom", universe = universe, ont = "MF", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  result_all_modules_CC      = compareCluster(geneClusters = ll.modules, fun = "enrichGO_custom", universe = universe, ont = "CC", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  result_all_modules_KEGG    = compareCluster(geneClusters = ll.modules.ENTREZ, fun = "enrichKEGG"   , universe = universe.ENTREZ, organism = "hsa"  , keyType = "kegg"     , pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  result_all_modules_Pathway = compareCluster(geneClusters = ll.modules.ENTREZ, fun = "enrichPathway", universe = universe.ENTREZ, organism = "human"                       , pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  
  final = dplyr::lst(result_all_modules_BP,result_all_modules_MF,result_all_modules_CC, result_all_modules_KEGG, result_all_modules_Pathway)
}


###Begin Main Code Here----
options(java.parameters = "-Xmx8000m")  #Required for xlsx package
#library(WGCNA)
library(SummarizedExperiment)
library(purrr)

set.seed(123)

library(limma)
library(dplyr)
library(pheatmap)
library(gtools)       #Optional. Just to sort path of files read
library(org.Hs.eg.db)
library(magrittr)
library(ggplot2)
library(ggrepel)
#library(directlabels)
library(patchwork)
library(tidyr)
#library(cowplot)
library(ggtext)
library(SummarizedExperiment)
library(limma)
library(future)
library(purrr)
library(furrr)
library(magrittr)
library(biomaRt)
library(tidyverse)

#Define bins/lists for different enrichments
our_bins          = c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")
patho_bins        = c("AD","ADHD","ALS","ASD","BIP","CD","MDD","MS","OCD","PD","PTSD","RA","SA","UC")
our_bins_complete = c(our_bins,
                      "Apua_CAUDATE_sczd", "Clozapine_mouse", "DEGs_human_Sousa", "Fromer_nosva", "Haloperidol", "Jaffe_DLPFC_HIPPO_sczd", "Jaffe_DLPFC_sczd", "Jaffe_HIPPO_sczd", "Jaffe_sczd_2018", "Risperidone_aripripazole_blood", #DEGs bin
                      "DMGs_Hannon", "DMGs_Jaffe", "DMGs_Kinoshita", "DMGs_Montano", "DMGs_Numata", "DMGs_Wockner", #DMGs bin
                      "Finan", "IDG_Schizo_Tclin", "Santos", "Santos_Tclin", "Wang", #Druggable_genes bin
                      "Gandal","Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Jaffe_HIPPO_pgc2.clozuk", "Jaffe_DLPFC_pgc2", "Jaffe_HIPPO_pgc2", "Hall","Paquola_CAUDATE", #TWAS bins
                      "LoF"   #LOF
)
our_networks = c("caudate","caudate__.1.25","caudate__25.50","caudate__50.100",
                 #"dentate",
                 "dentate.noQSVAremoved","dentate.QSVAremoved",
                 "hippo.noQSVAremoved","hippo.QSVAremoved",
                 "dlpfc","dlpfc__.1.6","dlpfc__6.25","dlpfc__25.50","dlpfc__50.100",
                 "hippo","hippo__.1.6","hippo__6.25","hippo__25.50","hippo__50.100",
                 "rse.pwr4.protectDx","rse.pwr14.protectDx", "stemcell.protectDx.pwr4","stemcell.protectDx.pwr14",
                 "rse.pwr4","rse.pwr14","stemcell.pwr4","stemcell.pwr14",
                 "stemcell.noConsensus.pwr4","stemcell.noConsensus.pwr14",
                 "stemcell.protectDx.noConsensus.pwr4","stemcell.protectDx.noConsensus.pwr14") 

published_networks = sort(c("Pergola2017","Pergola2019","Pergola2020",
                            "Fromer2016_case", "Fromer2016_control", "Gandal2018a", "Gandal2018b", "Gandal2018b_cs",#"Gandal2018", "Gandal201PE", "Gandal2018PE_cs"
                            "Radulescu2020","Walker2019", "Werling2020"))


##Read updated network names for consistent naming
new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\old-to-new Network names all.csv") %>% tibble::deframe()
new.names[grepl("Gandal2018$"     ,names(new.names))] = "Gandal2018a"
new.names[grepl("Gandal2018PE$"   ,names(new.names))] = "Gandal2018b"
new.names[grepl("Gandal2018PE_cs$",names(new.names))] = "Gandal2018b_cs"
new.names = c(new.names,c("Gandal2018a"="Gandal2018a","Gandal2018b"="Gandal2018b","Gandal2018b_cs"="Gandal2018b_cs"))
new.names = c(new.names,c("Gandal2018a"="Gandal2018a","Gandal2018b"="Gandal2018b","Gandal2018b_cs"="Gandal2018b_cs"))


#Read list of reference files for different pathologies [AD/ADHD/ALS/ASD/BIP/CD/MDD/MS/OCD/PD/PTSD/RA/SA/SCZ/SCZ.neg/UC][all.biotypes/protein.coding]
all.patho = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\all.patho.lists_grch38[PGC125].rds")
all.patho = all.patho[!grepl("75",names(all.patho))]
all.patho = unlist(unname(all.patho), recursive = F)

#Get the updated SCZ PGC120 lists and replace appropriate bins in the all.patho
pp = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\PGC.gene.lists_grch38[120 PGC ensemblIDs].RData", envir = pp)
all.patho$PGC3.all.biotypes$PGC        = pp$PGC3.all.biotypes$PGC       #Replace first/PGC AB bin only
all.patho$PGC3.negative.all.biotypes   = pp$PGC3.negative.all.biotypes  #Replace all negative AB bins 
all.patho$PGC3.protein.coding$PGC      = pp$PGC3.protein.coding$PGC     #Replace first/PGC PC bin only
all.patho$PGC3.negative.protein.coding = pp$PGC3.negative.protein.coding#Replace all negative PC bins

#Read further reference lists for DEGs/DMGs/Druggable_genes/TWAS/LOF etc
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\DEGs.RData"                       , envir = DEGS          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\DMGs.RData"                       , envir = DMGS          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\Druggable_genes.RData"            , envir = DRUG          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\TWAS.RData"                       , envir = TWAS          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\TWAS.universe.RData"              , envir = TWAS.universe <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\pleiotropic_Lof_Denovo_CNVs.RData", envir = LOF           <- new.env())
li_final = c(all.patho, as.list(DEGS),as.list(DMGS),as.list(DRUG),as.list(LOF))

##Read Hartl2021 gene list
Hartl2021 =readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Hartl2021\\Hartl2021-gene-module-list.rds")
#gm_AB = c(gm_AB, Hartl2021)
Hartl2021 = Hartl2021 %>% map(~{
  names(.x) = gsub("WHOLE_BRAIN\\.","BW.",names(.x))
  .x
})
module.color0 = WGCNA::standardColors(length(unique(unlist(map_depth(Hartl2021,1,names))))) %>% set_names(unique(unlist(map_depth(Hartl2021,1,names))))
module.color0["grey"] = "grey" #A named vector of Hartl2021 modules and psuedo-colors names picked from WGCNA colors for plotting. 
  
##Read gene list with all networks nets: parsed/non-parsed/published/sample-matched/metanets/replication/prenatal
gm_AB <- readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Paper new analysis/Revision/metanets/gm_AB_final_revision_metanets.rds")
gm_AB = gm_AB %>% map(~{
  names(.x) = gsub("WHOLE_BRAIN\\.","BW.",names(.x))
  .x
})
gm_AB = c(gm_AB,Hartl2021)
gm_AB = gm_AB[unique(c(names(gm_AB),names(Hartl2021)))]

names(gm_AB)[grep("Gandal2018$"     ,names(gm_AB))] = "Gandal2018a"
names(gm_AB)[grep("Gandal2018PE$"   ,names(gm_AB))] = "Gandal2018b"
names(gm_AB)[grep("Gandal2018PE_cs$",names(gm_AB))] = "Gandal2018b_cs"
names(gm_AB)[grep("META Perinatal"  ,names(gm_AB))] = "Replication Perinatal"
names(gm_AB)[grep("META Juvenile"   ,names(gm_AB))] = "Replication Juvenile"

#Keep ensemblIDs converted to Gene Symbol for future reference
allgenes = unique(unlist(gm_AB))
mapping  = clusterProfiler::bitr(allgenes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

#Get biotype for genes using biomart
gene_grch38 = geneMap_fun(genes = allgenes, ensembl = "dec2017.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gene_gc_content","entrezgene","version"))
gene_grch38 = gene_grch38[!duplicated(gene_grch38$ensembl_gene_id),]
gene_grch38$ensembl    = gene_grch38$ensemblID    = gene_grch38$ensembl_gene_id
gene_grch38$hgnc       = gene_grch38$Symbol       = gene_grch38$external_gene_name
gene_grch38$gene_type = gene_grch38$gene_biotype

#Get lists of protein coding only genes
gm_PC = map(gm_AB,~ {
  tmp = map(.x, ~{intersect(.x,gene_grch38$ensembl_gene_id[gene_grch38$gene_biotype == "protein_coding"])})
  tmp = tmp[lengths(tmp)>20]  ##Discarding modules with less than 21 genes when restricting to only protein.coding genes
})



# ####################################
# #######################################

###All network modules enrichment (hypergeo) for SCZ PGC3 (All biotype lists)
res_AB = imap(li_final[!grepl("protein.coding",names(li_final))],~{
  ref.list = .x
  ref.name = .y
  print(ref.name)
  
  imap(gm_AB,~{
    tmp = Enrich(ll.target = ref.list, ll.components = .x)
  })
}) %>% transpose()

###All network modules restricted to PC genes only, enrichment (hypergeo) for SCZ PGC3 (protein coding lists)
res_PC = imap(li_final[grepl("protein.coding",names(li_final))],~{
  ref.list = .x
  ref.name = .y
  print(ref.name)
  
  imap(gm_PC,~{
    tmp = Enrich(ll.target = ref.list, ll.components = .x)
  })
}) %>% transpose()


## TWAS enrichment with TWAS specific universe
res_AB_TWAS = imap(gm_AB,~{
  tmp = Enrich.TWAS(ll.target = as.list(TWAS$TWAS), ll.components = .x, ll.universe = as.list(TWAS.universe$TWAS.universe), target.name = ref.list.name)
})

res = pmap(list(res_AB, res_PC, res_AB_TWAS),~ c(..1,..2,TWAS = list(..3)))


##GO enrichment
#plan(multisession, workers = 7)
res_AB_GO = future_imap(gm_AB,~ {
  
  ##Use a modified clusterProfiler::enrichGO function to speed up enrichments (enrichGO_custom) by reading GO_DATA_ENSEMBL_XX files from environment instead of disk
  load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\GO_DATA_ENSEMBL.RData",envir = .GlobalEnv)
  source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/Reference files/enrichGO_custom_script.R")
  
  GO_enrich(Net = .x)
})


##Cell specificity enrichment
##
oe = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\Cell_specificity.RData", envir = oe)
source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/Reference files/Cell_specificity_function_script.R")
res_AB_CellSpecificity = imap(gm_AB,~ {
  Cell_Specificity_function(Net = .x, network.name = .y, genemap = gene_grch38, specificity = oe$Cell_specificity$human$DRONC_human)
})

##MAGMA enrichment
source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/MAGMA_enrich_preprocess_script.R")
source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/MAGMA_enrich_script.R")
res_AB_MAGMA = iwalk(gm_AB,~ {
  #MAGMA_enrich_preprocess(genemap = gene_grch38,pathology = "SCZ")
  
  MAGMA_enrich(Net = .x,
               pathology        = "SCZ",
               Net_name         = .y,
               genemap          = gene_grch38,
               target.directory = paste0(getwd(),"/SCZ/"),
               output.directory = paste0(getwd(),"/SCZ/"),
               magma.directory  = getwd())
})


# #Read MAGMA results from disk
# ll = list.files(path = "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Enrichments",pattern = '.*_MAGMA.RData', recursive = T,full.names = T)
# names(ll) = gsub("_"," ",gsub("\\.RData|_MAGMA","",basename(ll)))
# 
# nm = map_dfr(.id = "net", ll,~{
#   ee = new.env()
#   load(.x,ee)
#   tmp1 = ee$PGC_MAGMA_out
#   tmp2 = map_dfr(.id = "biotype",tmp1,~{
#     bin = .x
#     map_dfr(.id = "bins",bin,~.x$PGC3)})
# })
# nm1 = nm %>% pivot_wider(id_cols = c(net,VARIABLE), values_from = "P_bonf", names_from = c(biotype,bins))


#########################
#########################
###Save results on disk
#save(res,res_AB_GO, res_AB_CellSpecificity, file = "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\all_results(metanets).RData")

##Load saved results from disk
load("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Paper new analysis/Revision/metanets/all_results(metanets).RData")

############
#Compile results from different enrichments into a df with network_name, module_name, enrichment and a list-col of results
all.files = res %>% imap_dfr(.id = "network",~{
  map_dfr(.id = "enrichment",.x,~{
    load_li = list(
      matrix.pval   = .x$enrich.matrix.pval,
      matrix.bonf   = .x$enrich.matrix.bonf,
      modules       = .x$ll.components,
      hits          = .x$hit.genes
    )
    
    load_li[["hits.count.df"]] = load_li$hits %>% purrr::transpose() %>% 
      map_depth(2,~length(.)) %>% map_dfr(.id = "bins", ~.) %>% data.frame(row.names = "bins", check.names = F)
    
    load_li[["hits.genes.df"]] = load_li$hits %>% imap_dfc(.id = "modules",~ {
      dd = tibble::enframe(.x)[,"value"] %>% set_names(.y)
      return(dd)
    }) %>% as.data.frame()
    
    rownames(load_li[["hits.genes.df"]]) = rownames(load_li[["hits.count.df"]])
    
    df.internal = cbind(modules = rownames(load_li$matrix.bonf), module.length = lengths(load_li$modules), bonf. = load_li$matrix.bonf, pvals. = load_li$matrix.pval, hits. = load_li$hits.count.df, hitgenes. = load_li$hits.genes.df, stringsAsFactors=F)
    
    
    fold.factor.no.grey   = colSums(df.internal[df.internal$modules != "grey", grep("hits..", colnames(df.internal))])/sum(df.internal[df.internal$modules != "grey","module.length"])
    fold.factor.all.genes = colSums(df.internal[                             , grep("hits..", colnames(df.internal))])/sum(df.internal[                             ,"module.length"])
    
    if (any(fold.factor.all.genes == 0)) {warning("check Fold Factor");}; fold.factor.all.genes[fold.factor.all.genes == 0] = 1;  #if(!fold.factor.no.grey)   fold.factor.no.grey   =1
    if (any(fold.factor.no.grey   == 0)) {warning("check Fold Factor")}; fold.factor.no.grey  [fold.factor.no.grey   == 0] = 1;  #if(!fold.factor.all.genes) fold.factor.all.genes =1
    Fold.df.no.grey   = df.internal[, grep("hits..", colnames(df.internal))] %>% sweep(1,df.internal$module.length,"/") %>% sweep(2,fold.factor.no.grey  ,"/")
    Fold.df.all.genes = df.internal[, grep("hits..", colnames(df.internal))] %>% sweep(1,df.internal$module.length,"/") %>% sweep(2,fold.factor.all.genes,"/")
    Fold.df           = Fold.df.no.grey
    Fold.df[rownames(Fold.df.no.grey) == "grey",] = Fold.df.all.genes[rownames(Fold.df.all.genes) == "grey",]
    colnames(Fold.df) = gsub("hits","Fold.hits",colnames(Fold.df))
    
    df.internal = cbind(df.internal, Fold.df)
    #df.internal = cbind(df.internal, Fold = (df.internal[,grep("hits..", colnames(df.internal))]/df.internal$module.length)/fold.factor)
    #browser()
    df.tmp = df.internal[,grep("pvals..", colnames(df.internal))]
    df.tmp[df.tmp == 1] = 1 - 1e-15
    df.tmp = apply(df.tmp,2,qnorm)
    colnames(df.tmp) = gsub("pvals","zvals", colnames(df.tmp))
    df.internal    = cbind(df.internal, df.tmp)
    df.tmp = NULL
    
    return(df.internal)
    
  })
}) %>% dplyr::select(c(network,enrichment, modules, module.length), matches(paste0("bonf..",our_bins_complete,"$")), matches(paste0("pvals..",our_bins_complete,"$")), matches(paste0("hitgenes..",our_bins_complete,"$")), matches(paste0("hits..",our_bins_complete,"$")), matches(paste0("Fold.hits..",our_bins_complete)), contains(paste0("zvals",our_bins_complete,"$")))

new.names = strsplit2(unique(all.files$network),"\\.\\.\\.")[,1] %>% set_names(.,.)

#Turn df with list-col above into a long-form df suitabel for plotting
all.files_long = all.files %>% dplyr::select(-contains("hitgenes..")) %>% 
  pivot_longer(cols = -c(network,enrichment,modules,module.length),names_to = "test", values_to = "val") %>%
  drop_na(val) %>%
  separate(enrichment, into = c("enrichment","biotype"), sep = "\\.(?=(all.biotypes|protein.coding))", remove = F, fill = "warn") %>%
  separate(test      , into = c("test","bins")         , sep = "\\.\\."                              , remove = F)   %>%
  mutate(
    network            = case_when(
      grepl("Gandal2018$"     ,network) ~ "Gandal2018a",
      grepl("Gandal2018PE$"   ,network) ~ "Gandal2018b",
      grepl("Gandal2018PE_cs$",network) ~ "Gandal2018b_cs",
      TRUE                              ~ network),
    network.type = case_when(grepl("caudate|dlpfc|hippo|dentate",network)~"our",#ifelse(grepl("caudate|dlpfc|hippo|dentate",network),"our","published"),
                             network %in% published_networks~"published",
                             network %in% names(Hartl2021) ~ "Hartl"),
    new.network        = ifelse(network.type %in% "Hartl",paste0("Hartl2021_",network),new.names[network]),
    new.network.module = paste0(new.network,">",modules),
    biotype            = replace_na(biotype,"all.biotypes"),
    module.color       = ifelse(network.type %in% "Hartl",module.color0[modules],modules)
  )
all.files_long$network.type[all.files_long$network %in% names(Hartl2021)] = "Hartl"
all.files_long$network.type[all.files_long$network %in% c("Hartl2021_BRNCTX","BRNCTX")] = "published"


all.files_long$new.network = all.files_long$network
all.files_long$new.network.module = paste0(all.files_long$network,">",all.files_long$modules)

#Save long-form file with modulewise enrichment results for SCZ/Pathologies/DEGs/DMGs/Druggable genes/LOF/TWAS
saveRDS(all.files_long, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\all.files_long_metanets.rds")

###Prepare enrichment results file for visualisation script. For several non-SCZ enrichments, enrichments in multiple bins are combined by tissue types. 
cols_from_p_all = c("shortname",              "filename",               "list.type",              "biotype",                "tissue_age",             "enrichment",
                    "class",                  "net.type",               "new_tissue_age",         "original_tissue_age",    "modules",                "module.length",
                    "bins",                   "bonf.vals",              "pvals.vals",             "hits.count",             "hitgenes",               "Fold.vals",
                    "Zscore.vals",            "sig.modules",            "tissue_age_sig.modules", "temp_tissue_age",        "module.color")

#Reformatting to template df by adding NA columns to be used with legacy visualisation script
all.files_long_metanets_test = all.files_long %>% pivot_wider(id_cols = c(biotype, enrichment, new.network, new.network.module, network,modules, bins, module.length), names_from = test, values_from = val ) %>%
  mutate(shortname = NA, filename = NA, list.type = enrichment, biotype = biotype, tissue_age = new.network, enrichment = enrichment, class = NA, net.type = NA, new_tissue_age = new.network, original_tissue_age = new.network, modules = modules, module.length = module.length, bins = bins, bonf.vals = bonf, pvals.vals = NA, hits.count = hits, hitgenes = NA, Fold.vals = NA, Zscore.vals = NA, sig.modules = NA, tissue_age_sig.modules = NA, temp_tissue_age = NA, module.color = modules) %>%
  filter(!is.na(bonf)) %>% dplyr::select(all_of(cols_from_p_all))

all.files_long_metanets_test = all.files_long_metanets_test %>% mutate(
  list.type = case_when(
    list.type == "AD"              ~ "ad"              ,
    list.type == "ADHD"            ~ "adhd"            ,
    list.type == "ALS"             ~ "als"             ,
    list.type == "ASD"             ~ "asd"             ,
    list.type == "BIP"             ~ "bip"             ,
    list.type == "CD"              ~ "cd"              ,
    list.type == "DEGs"            ~ "DEGs"            ,
    list.type == "DMGs"            ~ "DMGs"            ,
    list.type == "Druggable_genes" ~ "Druggable_genes" ,
    list.type == "LOF"             ~ "LOF"             ,
    list.type == "MDD"             ~ "mdd"             ,
    list.type == "MS"              ~ "ms"              ,
    list.type == "OCD"             ~ "ocd"             ,
    list.type == "PD"              ~ "pd"              ,
    list.type == "PTSD"            ~ "ptsd"            ,
    list.type == "RA"              ~ "ra"              ,
    list.type == "SA"              ~ "sa"              ,
    list.type == "PGC3"            ~ "PGC3"            ,
    list.type == "PGC3.negative"   ~ "PGC3.negative"   ,
    list.type == "TWAS"            ~ "TWAS"            ,
    list.type == "UC"              ~ "uc"              ,
    TRUE                           ~ list.type),
  enrichment = case_when(
    enrichment == "AD"              ~ "ad"              ,
    enrichment == "ADHD"            ~ "adhd"            ,
    enrichment == "ALS"             ~ "als"             ,
    enrichment == "ASD"             ~ "asd"             ,
    enrichment == "BIP"             ~ "bip"             ,
    enrichment == "CD"              ~ "cd"              ,
    enrichment == "DEGs"            ~ "DEGs"            ,
    enrichment == "DMGs"            ~ "DMGs"            ,
    enrichment == "Druggable_genes" ~ "Druggable_genes" ,
    enrichment == "LOF"             ~ "LOF"             ,
    enrichment == "MDD"             ~ "mdd"             ,
    enrichment == "MS"              ~ "ms"              ,
    enrichment == "OCD"             ~ "ocd"             ,
    enrichment == "PD"              ~ "pd"              ,
    enrichment == "PTSD"            ~ "ptsd"            ,
    enrichment == "RA"              ~ "ra"              ,
    enrichment == "SA"              ~ "sa"              ,
    enrichment == "PGC3"            ~ "SCZ"             ,
    enrichment == "PGC3.negative"   ~ "SCZ"             ,
    enrichment == "TWAS"            ~ "TWAS"            ,
    enrichment == "UC"              ~ "uc"              ,
    TRUE                            ~ enrichment)
)

#For some enrichments, we group some bins
temp = all.files_long_metanets_test %>% group_by(list.type, enrichment, biotype, new_tissue_age)
gk = group_keys(temp)
temp1 = pmap_dfr(list(group_split(temp), gk$list.type, gk$enrichment, gk$biotype, gk$new_tissue_age),~ {
  dat        = ..1
  enrichment = ..3
  biotype    = ..4
  tissue     = ..5

  ID = paste(enrichment,biotype,tissue)
  print(ID)

  if (enrichment %in% c("SCZ", "ad", "adhd", "als", "asd", "bip", "cd", "mdd", "ms", "ocd", "pd", "ptsd", "ra", "sa", "uc")){
    out = NULL
  }

  if(enrichment %in% "DEGs") {
    out = dat %>% filter(
      (grepl("CN"                ,new_tissue_age) & grepl("Apua_CAUDATE_sczd",bins)) |
        (grepl("HP"              ,new_tissue_age) & grepl("Jaffe_HIPPO_sczd" ,bins)) |
        (grepl("DG"              ,new_tissue_age) & grepl("Jaffe_HIPPO_sczd" ,bins)) |
        (grepl("DLPFC"           ,new_tissue_age) & grepl("Jaffe_DLPFC_sczd" ,bins)) |
        (!grepl("CN|HP|DG|DLPFC" ,new_tissue_age) & grepl("Jaffe_DLPFC_sczd" ,bins))
    )
    out$bins = "DEGs SCZ" #DEGs SCZ uses enrichment bin based upon network tissue type
  }

  if (enrichment %in% c("DMGs")){
    vec = c("DMGs_Hannon", "DMGs_Jaffe", "DMGs_Kinoshita", "DMGs_Montano", "DMGs_Numata", "DMGs_Wockner")
    out = dat %>% filter(bins %in% vec)
    out = out %>% group_by(list.type, enrichment, biotype, new_tissue_age, modules) %>%
      mutate(bonf.vals = metap::sumlog(bonf.vals)$p, bins = "DMGs general") %>% dplyr::slice_head(n=1) %>% ungroup #DMGs general uses sumlog of enrichment from multiple bins

  }

  if (enrichment %in% c("pleiotropic_Lof_Denovo_CNVs")){
    vec = c("LoF")
    #vec = c("LoF")
    out = NULL #out = dat %>% filter(bins %in% vec)
    #dat$bins = factor(dat$bins)
  }

  if (enrichment %in% c("Druggable_genes")){
    vec = c("Finan", "IDG_Schizo_Tclin", "Santos", "Santos_Tclin", "Wang")
    out = dat %>% filter(bins %in% vec)

    out = out %>% group_by(list.type, enrichment, biotype, new_tissue_age, modules) %>%
      mutate(bonf.vals = metap::sumlog(bonf.vals)$p, bins = "Drug general") %>% dplyr::slice_head(n=1) %>% ungroup  #Drug general uses sumlog of enrichment from multiple bins
  }


  if (enrichment %in% c("TWAS")){
    vec = c("Gandal","Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Jaffe_HIPPO_pgc2.clozuk", "Jaffe_DLPFC_pgc2", "Jaffe_HIPPO_pgc2", "Hall","Paquola_CAUDATE")
    out = dat %>% filter(bins %in% vec)

    out = out %>% filter(
      (grepl("caudate"                       ,new_tissue_age) & grepl("Paquola_CAUDATE"    ,bins))                                  |
        (grepl("hippo"                       ,new_tissue_age) & grepl("Jaffe_HIPPO_pgc2.clozuk" ,bins))                             |
        (grepl("dentate"                     ,new_tissue_age) & grepl("Jaffe_HIPPO_pgc2.clozuk" ,bins))                             |
        (grepl("dlpfc"                       ,new_tissue_age) & (bins %in% c("Gandal","Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Hall"))) |
        (!grepl("caudate|hippo|dentate|dlpfc",new_tissue_age) & (bins %in% c("Gandal","Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Hall")))
    )


    out = out %>% group_by(list.type, enrichment, biotype, new_tissue_age, modules) %>% group_split %>% map_df(~{
      if (dim(.x)[1] > 1){
        #if (ID ==  "TWAS all.biotypes ALL") browser()
        res = .x %>% mutate(bonf.vals = ifelse(any(is.na(bonf.vals)), NA,metap::sumlog(bonf.vals)$p), bins = "TWAS general") %>% dplyr::slice_head(n=1)   #TWAS general uses sumlog of enrichment from multiple bins, bins select according to network tissue type
        return(res)
      } else {
        res = .x %>% mutate(bins = "TWAS general")
        return(res)
      }
    })
  }

  return(out)
})

all.files_long_metanets_test = rbind(all.files_long_metanets_test, temp1)
#Save long-form file with modulewise enrichment results for SCZ/Pathologies/DEGs/DMGs/Druggable genes/LOF/TWAS summarised enrichments suitable for visualization
saveRDS(all.files_long_metanets_test, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\all.files_long_metanets_test.rds")

#Read list of module previously identified in publications for SCZ
priortised_modules = list(
  "Fromer2016_case"   = c("yellow","red","blue","lightyellow","greenyellow","cyan","grey60")              ,              
  "Fromer2016_control"= c("magenta","blue","greenyellow","salmon","turquoise","bisque")                   ,
  "Radulescu2020"     = c("blue","brown","green")                                                         , 
  "Walker2019"        = c("red","blue")                                                                   ,                                                                   
  "Werling2020"       = c("yellow","red","brown","magenta","cyan")                                        ,                                        
  "Li2018"            = c("green","skyblue3","blue","violet","royalblue","black")                         ,                         
  "Gandal2018PE"      = c("lightyellow","brown","red","black","blue","darkred","saddlebrown","turquoise") , 
  "Gandal2018PE_cs"   = c("pink","lightcyan")                                                             ,                                                             
  "Gandal2018"        = c("turquoise","green","yellow","salmon","purple")                                 ,
  "Gandal2018b"       = c("lightyellow","brown","red","black","blue","darkred","saddlebrown","turquoise") , 
  "Gandal2018b_cs"    = c("pink","lightcyan")                                                             ,                                                             
  "Gandal2018a"       = c("turquoise","green","yellow","salmon","purple")                                 ,
  "Pergola2017"       = c(NA)                                                                             ,
  "Pergola2019"       = c("darkgreen")                                                                    ,
  "Pergola2020"       = c("darkorange")                                                                   ,
  "BRNACC"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNAMY"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNCBH"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,                                                                    
  "BRNCBL"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,                                                                    
  "BRNCDT"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNCTX"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,                                                                    
  "BRNCTXB24"         = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNCTXBA9"         = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNHIP"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNHYP"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BRNPUT"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,                                                                    
  "BRNSNA"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "ALL"               = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BGA"               = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "BROD"              = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "CEREBELLUM"        = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "CTX"               = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "NS.SCTX"           = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "STR"               = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "SUBCTX"            = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                           ,
  "WHOLE_BRAIN"       = c("BW.M1","BW.M4","WHOLE_BRAIN.M1", "WHOLE_BRAIN.M4","CEREBELLUM.M3","CTX.M3")                                 
)

#Generate a long-form file with only SCZ protein-coding enrichments (excludes grey modules) 
all.files_long_SCZ = all.files_long %>% filter(enrichment %in% "PGC3" & biotype %in% "protein.coding" & !(modules %in% "grey") & test %in% "bonf") %>%
  mutate(significant = ifelse(-log10(val) > -log10(0.05),T,F)) %>%
  rowwise() %>%
  mutate(priortised = ifelse(modules %in% priortised_modules[[as.character(network)]],T,F)) %>%
  ungroup()

#Save long-form file with modulewise enrichment results for SCZ
saveRDS(all.files_long_SCZ, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\all.files_long_SCZ_metanets.rds")


#Generate a long-form file with only negative SCZ protein-coding enrichments 
all.files_long_SCZ.negative.pvals = all.files_long %>% filter(enrichment %in% "PGC3.negative" & biotype %in% "protein.coding" & !(modules %in% "grey") & test %in% "pvals") %>%
  mutate(significant = ifelse(-log10(val) > -log10(0.05),T,F)) %>%
  rowwise() %>%
  mutate(priortised = ifelse(modules %in% priortised_modules[[as.character(network)]],T,F)) %>%
  ungroup()

#Make a list of SCZ risk modules based on 3 significantly enriched windows for SCZ.protein.coding
SCZ_risk_modules_df = all.files_long_SCZ %>% filter(!new.network %in%  c("DG","DG.QSVA","HP.noQSVA"))               %>% 
  mutate(final.val = -log10(val)#, 
  )                         %>%
  group_by(enrichment,biotype,new.network.module)                                                                   %>%
  summarise(sig.bins = list(as.character(bins)[significant]),
            n.sig.bins = lengths(sig.bins),
            network.type = network.type[1])                                                                         %>%
  ungroup()                                                                                                         %>%
  pivot_wider(id_cols = c(new.network.module, network.type), values_from = c(sig.bins, n.sig.bins), names_from = c(enrichment,biotype)) %>%
  separate(new.network.module, into = c("new.network","module"), sep = ">", remove = F)                             %>%
  mutate(isSCZriskmodules   = ifelse(n.sig.bins_PGC3_protein.coding >=3,1,0))  %>%
  filter(isSCZriskmodules == 1)

#Make a final list of SCZ risk modules based on 3 significantly enriched windows for SCZ.protein.coding. Exclude Prenatal networks and Hartl networks except BRNCTX
SCZ_risk_modules_df = SCZ_risk_modules_df[!SCZ_risk_modules_df$network.type %in% "Hartl" & !grepl("Prenatal",SCZ_risk_modules_df$new.network.module),]
##Add back Walker2019>brown and Li2018>yellow module for the consensus analysis only 
SCZ_risk_modules_df.added = rbind(SCZ_risk_modules_df, data.frame(new.network.module = c("Walker2019>brown","Li2018>yellow"), new.network = c("Walker2019","Li2018"), module = c("brown","yellow"), network.type = NA, sig.bins_PGC3_protein.coding = NA, n.sig.bins_PGC3_protein.coding = NA, isSCZriskmodules = 0))

##Save Risk module list
#saveRDS(SCZ_risk_modules_df , file =  "SCZ.risk.modules.df.rds")
##saveRDS(SCZ_risk_modules_df %>% pull(new.network.module), file =  "SCZ.risk.modules.rds")

##For negative risk modules, if a network has no significant negative module, find most significant module from (min pval) nn_200
SCZ_risk_modules_df.correct.negative = all.files_long_SCZ.negative.pvals %>% filter(!new.network %in%  c("DG","DG.QSVA","HP.noQSVA") & bins %in% "kbp_200" & modules != "grey")               %>% 
  mutate(final.val = -log10(val)#, 
         #sig.module = ifelse(final.val >-log10(0.05),1,0), 
         #bins = factor(bins, levels = our_bins_complete),
         #new.network = new.names[network]
  ) %>%  
  group_by(new.network)  %>% 
  arrange(desc(val),desc(module.length)) %>% #group_split()
  slice_min(val)

#Save list of raw input files needed for the consensus analysis
saveRDS(list(gene_module_list_AB = gm_AB, gene_module_list_PC = gm_PC, risk.module.positive.added = SCZ_risk_modules_df.added, risk.module.negative.corrected = SCZ_risk_modules_df.correct.negative),
        file = "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\new_consensus_raw_data_metanets (corrected).rds")


#################################
#################################
##### Code for sankey plots (hitgenes) -----
make_sankeys = function(all.files0 = all.files, SCZ_risk_modules_df0 = SCZ_risk_modules_df){
  
  #Prepare nodes df and edges df for sankey plots by subsetting
  sankey = all.files0 %>% dplyr::select(network, enrichment, modules, module.length, hitgenes..kbp_200) %>% 
    mutate(new.network        = network,
           new.network.module = paste0(new.network,">",modules),
           hitgenes..kbp_200  = hitgenes..kbp_200 %>% set_names(new.network.module))                   %>%
    filter(enrichment %in% "PGC3.all.biotypes" & grepl(" ",new.network))                               %>%
    split(.,.$new.network)
  
  sankey.hitgenes    = map(sankey,"hitgenes..kbp_200") %>% purrr::flatten(.)  #CAUTION: Purrr::flatten is changing in upcoming purrr versions
  
  sankey.pairs = list( #Define age-transitions networks for each tissue 
    DLPFC = list(
      t1 = crossing(source = sankey$`DLPFC Perinatal`$new.network.module, target = sankey$`DLPFC Juvenile`$new.network.module)  ,
      t2 = crossing(source = sankey$`DLPFC Juvenile`$new.network.module , target = sankey$`DLPFC Adult`$new.network.module)     ,
      t3 = crossing(source = sankey$`DLPFC Adult`$new.network.module    , target = sankey$`DLPFC Older Adult`$new.network.module)
    ),
    
    HP = list(
      t1 = crossing(source = sankey$`HP Perinatal`$new.network.module   , target = sankey$`HP Juvenile`$new.network.module)  ,
      t2 = crossing(source = sankey$`HP Juvenile`$new.network.module    , target = sankey$`HP Adult`$new.network.module)     ,
      t3 = crossing(source = sankey$`HP Adult`$new.network.module       , target = sankey$`HP Older Adult`$new.network.module)
    ),
    
    CN = list(
      t1 = crossing(source = sankey$`CN Juvenile`$new.network.module    , target = sankey$`CN Adult`$new.network.module)     ,
      t2 = crossing(source = sankey$`CN Adult`$new.network.module       , target = sankey$`CN Older Adult`$new.network.module)
    )
  ) %>% map_dfr(.id = "network",~ bind_rows(.x,.id = "transitions")) %>%
    rowwise()                                                        %>%
    mutate(hitgenes.intersect       = {list(intersect(sankey.hitgenes[[source]],sankey.hitgenes[[target]]))}) %>%
    ungroup()                                                                                                 %>%
    mutate(hitgenes.intersect.count = lengths(hitgenes.intersect))                                            %>%
    filter(hitgenes.intersect.count >0)
  
  sankey.nodes.list  = data.frame(node = unique(c(sankey.pairs$source,sankey.pairs$target))) %>% separate(node, into = c("new.network","module"),sep =">",remove = F)
  
  ##Plot sankeys for three tissues
  walk(c("DLPFC","HP","CN"),~{
    tissue = .x
 
    nodes.df = sankey.nodes.list %>% filter(grepl(tissue, node)) %>%
      mutate(module2 = ifelse(node %in% SCZ_risk_modules_df0$new.network.module, toupper(module),"")) %>%
      .[c(grep(paste0(tissue," Perinatal"),.$node),grep(paste0(tissue," Juvenile"),.$node),grep(paste0(tissue," Adult"),.$node),grep(paste0(tissue," Older Adult"),.$node)),]
    edges.df = sankey.pairs %>% filter(network %in% tissue) %>%
      mutate(source.index  = match(source, nodes.df$node)-1 ,
             target.index  = match(target, nodes.df$node)-1,
             source.module = strsplit2(source,">")[,2],
             target.module = strsplit2(target,">")[,2],)
    
    library(sankeyD3)
    
    mod.colors = paste0('"',unique(tolower(nodes.df$module)),'"', collapse = ',')
    my_color = paste0("d3.scaleOrdinal() .domain([", mod.colors, "]) .range([", mod.colors, "])")
    
    
    sankey1 = sankeyD3::sankeyNetwork(Links          = edges.df, 
                                      Source         = "source.index", 
                                      Target         = "target.index", 
                                      Value          = "hitgenes.intersect.count", 
                                      LinkGroup      = "source.module", 
                                      linkGradient   = F, 
                                      linkColor      = "lightgrey",
                                      
                                      Nodes          = nodes.df, 
                                      NodeGroup      = "module", 
                                      NodeID         = "module2", 
                                      showNodeValues = F, 
                                      #NodeValue = "totalgenes.count",
                                      nodeCornerRadius = 5,
                                      nodeLabelMargin  = 6,
                                      nodeStrokeWidth  = 1, 
                                      units            = "genes", 
                                      #NodePosX = "col.x.index",
                                      #nodePadding = 13, #scaleNodeBreadthsByString = T, nodeShadow = F,
                                      nodeWidth = 20, curvature = 0.5, linkOpacity = 0.4,
                                      #NodeFontSize = 15, 
                                      fontSize            = 21, 
                                      iterations          = 10, 
                                      colourScale         = my_color, 
                                      title               = paste0("SCZ risk genes (200kbp)"), 
                                      margin              = list(top = 20, bottom = 40,left = 75, right = 55), 
                                      numberFormat        = ",d",
                                      align               = "center", 
                                      dragX               = T, 
                                      dragY               = T, 
                                      zoom                = T, 
                                      highlightChildLinks = T, 
                                      doubleclickTogglesChildren = F, 
                                      linkType            = "path1", 
                                      xAxisDomain         = as.list(unique(nodes.df$new.network)),
                                      yOrderComparator = JS("function(a,b) {if (a.name == 'grey'||a.name == 'GREY') {return  -1;}; return a.y - b.y; }"),
                                      xScalingFactor = .95
    )
    
    #Increase font size of x-axis text
    sankey1 = htmlwidgets::onRender(sankey1, 'function(el){
                                alltext   = el.querySelector("g.x.axis");
                                alltext.setAttribute("font-size","19");
                                    }')
    #Update text for links for tooltips
    sankey1$x$links$test = paste0(edges.df$source.module," -> ",edges.df$hitgenes.intersect.count," genes -> ",edges.df$target.module)
    
    sankey1 <- htmlwidgets::onRender(
      sankey1,
      '
  function(el, x) {
  d3.selectAll(".link").select("*").text(function(d) { return d.test ; });
  }
  '
    )
    print(sankey1)
    
  })
  
  rm(sankey, sankey.hitgenes, sankey.pairs, sankey.nodes.list) 
}

make_sankeys(all.files, SCZ_risk_modules_df) #Makes sankey plots for DLPFC, CN, HIPPO hitgenes



