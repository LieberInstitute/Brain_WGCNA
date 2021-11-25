####################################################################################################
#                    Script to run various enrichments on gene_module_list:                        #
#                   SCZ, DEGs, DMGs, TWAS, GO, MAGMA, Cell Specificity etc.                        #
#                                                                                                  #
#                  Script to make sankey plots for DLPFC, HP and CN SCZ risk genes.                #
#                                                                                                  #
#                  Script to make binwise PGC3 enrichment scatterplots for networks.               #
#                                                                                                  #
#                          Prepare and format data for other visualisations                        #
#           "Regional-coexpression", "Age-period" and "Cell-population enrichment" studies         #
#                                                                                                  #
####################################################################################################
 

##Function to get gene annotation from Biomart
geneMap_fun = function(genes, ensembl, attributes){
  require(biomaRt)
  ensembl_archive = useMart("ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl", host=ensembl)
  Filters = "ensembl_gene_id"
  
  df = getBM(attributes = attributes,Filters,values = genes,mart = ensembl_archive)
  df = df[df$chromosome_name %in% c(1:22,"X","Y","M"),]
  df$chromosome_name = factor(df$chromosome_name,levels = c(1:22,"X","Y","M"))
  df = df[order(df$chromosome_name,df$start_position,df$end_position),]
  return(df)
}

#Function to perform hypergeometric enrichment on genesets
Enrich = function(ll.target, ll.components){
  ###ll.components.size argument is not needed
  
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
      #print(length(ll.components))
      #print(lengths(ll.components))
      # if (any(lengths(ll.components) < 10)) {
      #   print(names(ll.components)[lengths(ll.components)<10]);print(target.name);print(network_name)
      # }
      #View(ll.components)
      
      #if (any(ll.components.size > length(universe))) {print(target.name);print(names(ll.components)[ll.components.size > length(universe)])}
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
      
      # if (any(ll.components.size > length(universe))) {
      #   df$hit.pop[rownames(df) %in% names(ll.components)[ll.components.size > length(universe)]] = NA
      #   df$fail.pop[rownames(df) %in% names(ll.components)[ll.components.size > length(universe)]] = NA
      # }
      
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
      #df$drawn[rownames(df) == "grey"] = min(14750, df$drawn[rownames(df) == "grey"])
      #if (df$drawn[rownames(df) == "grey"] > 10390) {df$drawn[rownames(df) == "grey"] = NA}
      if (any(sapply(ll.components,length) > universe)) {
        #print(names(ll.components)[sapply(ll.components,length) > universe])
        
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
      
      #   hg[removed.modules] = NA
      #   hg = hg[match(all.modules.names, names(hg))]
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
    # if(sum(is.nan(hg)) > 0) print(target.name)
    #  hg[is.nan(hg)] = NA
    #View(hg,target.name)
    
    
    return(list(df = df, hit.genes = hit.genes, hg = hg))
  }, simplify = F, USE.NAMES = T) 
  
  xx = sapply(enrich.df.list,function(x) x$hg, simplify = F)
  
  enrich.matrix      = lapply(enrich.df.list,function(x) x$df)
  hit.genes          = lapply(enrich.df.list,function(x) x$hit.genes)
  #enrich.matrix.pval = data.frame(map_dfr(.id = "list",xx,~.x) %>% tibble::column_to_rownames("list") %>% t)
  enrich.matrix.pval = data.frame(sapply(enrich.df.list,function(x) x$hg), check.names = F)
  
  
  #View(enrich.matrix)
  #View(hit.genes)
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
  
  yy = map(enrich.matrix,~ .x$hyper.geo.fdr %>% set_names(rownames(.x)))
  enrich.matrix.fdr = map_dfr(.id = "list",yy,~.x) %>% tibble::column_to_rownames("list") %>% t
  #enrich.matrix.fdr  = as.data.frame(sapply(enrich.matrix,"[[", "hyper.geo.fdr"))
  #row.names(enrich.matrix.fdr) = row.names(enrich.matrix[[1]])
  
  yy = map(enrich.matrix,~ .x$hyper.geo.bonf %>% set_names(rownames(.x)))
  enrich.matrix.bonf = map_dfr(.id = "list",yy,~.x) %>% tibble::column_to_rownames("list") %>% t
  #enrich.matrix.bonf = as.data.frame(sapply(enrich.matrix,"[[", "hyper.geo.bonf"))
  #row.names(enrich.matrix.bonf) = row.names(enrich.matrix[[1]])
  

  return(dplyr::lst(enrich.matrix, enrich.matrix.pval, enrich.matrix.bonf, enrich.matrix.fdr, ll.components, hit.genes))
  
}


#Function to get GO enrichment
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
  
  result_all_modules_BP     = compareCluster(geneClusters = ll.modules, fun = "enrichGO_custom", universe = universe, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  result_all_modules_MF     = compareCluster(geneClusters = ll.modules, fun = "enrichGO_custom", universe = universe, ont = "MF", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  result_all_modules_CC     = compareCluster(geneClusters = ll.modules, fun = "enrichGO_custom", universe = universe, ont = "CC", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = method.adjust)
  
  final = dplyr::lst(result_all_modules_BP,result_all_modules_MF,result_all_modules_CC)
}

###Begin ----
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
library(directlabels)
library(patchwork)
library(tidyr)
library(cowplot)
library(ggtext)
library(SummarizedExperiment)
library(limma)
library(future)
library(purrr)
library(furrr)
library(magrittr)
library(biomaRt)
library(tidyverse)

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

##Read updated network names
new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\old-to-new Network names all.csv") %>% tibble::deframe()

##Read reference files
# #Read PGC3 reference list (all.biotypes/protein.coding/all.biotypees.negative/protein.coding.negative) for all bins
# load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\PGC.gene.lists_grch38[125 PGC ensemblIDs].RData", envir = PGC <- new.env())

#Read list of reference files for different pathologies [AD/ADHD/ALS/ASD/BIP/CD/MDD/MS/OCD/PD/PTSD/RA/SA/SCZ/SCZ.neg/UC][all.biotypes/protein.coding]
all.patho = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\all.patho.lists_grch38[PGC125].rds")
all.patho = all.patho[!grepl("75",names(all.patho))]
all.patho = unlist(unname(all.patho), recursive = F)

load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\DEGs.RData"                       , envir = DEGS          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\DMGs.RData"                       , envir = DMGS          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\Druggable_genes.RData"            , envir = DRUG          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\TWAS.RData"                       , envir = TWAS          <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\TWAS.universe.RData"              , envir = TWAS.universe <- new.env())
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\pleiotropic_Lof_Denovo_CNVs.RData", envir = LOF           <- new.env())
li_final = c(all.patho, as.list(DEGS),as.list(DMGS),as.list(DRUG),as.list(LOF))



##Read gene list: parsed/non-parsed/published/sample-matched networks
gm_AB = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (wide_form_test) (all networks)[grch38].rds")

###Read gene list: Sliding window (NC)
#gm_AB = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (sliding_window_NC)[grch38].rds")

###Read gene list: Sliding window (SchizoNew)
#gm_AB = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (sliding_window_SCZ)[grch38].rds")


allgenes = unique(unlist(gm_AB))
mapping  = clusterProfiler::bitr(allgenes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

gene_grch38 = geneMap_fun(genes = allgenes, ensembl = "jul2016.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gc_content","entrezgene","version"))
gene_grch38 = gene_grch38[!duplicated(gene_grch38$ensembl_gene_id),]
gene_grch38$ensembl    = gene_grch38$ensemblID    = gene_grch38$ensembl_gene_id
gene_grch38$hgnc       = gene_grch38$Symbol       = gene_grch38$external_gene_name
gene_grch38$gene_type = gene_grch38$gene_biotype

gm_PC = map(gm_AB,~ {
  tmp = map(.x, ~{intersect(.x,gene_grch38$ensembl_gene_id[gene_grch38$gene_biotype == "protein_coding"])})
  tmp = tmp[lengths(tmp)>20]  ##Discarding modules with less than 21 genes when restricting to only protein.coding genes
})



# ####################################
# #######################################
# ##All network modules  enrichment (hypergeo) for PGC3+consensus_list
res_AB = imap(li_final[!grepl("protein.coding",names(li_final))],~{
  ref.list = .x
  print(.y)
  
  imap(gm_AB,~{
    tmp = Enrich(ll.target = ref.list, ll.components = .x)
  })
}) %>% transpose()

res_PC = imap(li_final[grepl("protein.coding",names(li_final))],~{
  ref.list = .x
  print(.y)
  
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
##Use a modified clusterProfiler::enrichGO function to speed up enrichments (enrichGO_custom)
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\GO_DATA_ENSEMBL.RData",envir = .GlobalEnv)
source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/enrichGO_custom_script.R")


res_AB_GO = imap(gm_AB,~ {
  GO_enrich(Net = .x)
})


##Cell specificity enrichment
##
oe = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\Reference files\\Cell_specificity.RData", envir = oe)

source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/Cell_specificity_function_script.R")
res_AB_CellSpecificity = imap(gm_AB,~ {
  Cell_Specificity_function(Net = .x, network.name = .y, genemap = gene_grch38, specificity = oe$Cell_specificity$human$DRONC_human)
})

source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/MAGMA_enrich_preprocess_script.R")
source("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/MAGMA_enrich_script.R")
res_AB_MAGMA = iwalk(gm_AB,~ {
  #Cell_Specificity_function(Net = .x, network.name = .y, genemap = gene_grch38, specificity = oe$Cell_specificity$human$DRONC_human)
  #MAGMA_enrich_preprocess(genemap = gene_grch38,pathology = "SCZ")

  MAGMA_enrich(Net = .x,
               pathology        = "SCZ",
               Net_name         = .y,
               genemap          = gene_grch38,
               target.directory = paste0(getwd(),"/SCZ/"),
               output.directory = paste0(getwd(),"/SCZ/"),
               magma.directory  = getwd())
})

res_AB_HMAGMA = iwalk(gm_AB,~ {
  Cell_Specificity_function(Net = .x, network.name = .y, genemap = gene_grch38, specificity = oe$Cell_specificity$human$DRONC_human)
})


#########################
#########################
###Save results on disk
#save(res,res_AB_GO, res_AB_CellSpecificity, file = "C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/all_results.RData")

##Load saved results from disk
load("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Code for WGCNA paper/all_results.RData")


############
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

all.files_long = all.files %>% dplyr::select(-contains("hitgenes..")) %>% 
  pivot_longer(cols = -c(network,enrichment,modules,module.length),names_to = "test", values_to = "val") %>%
  drop_na(val) %>%
  separate(enrichment, into = c("enrichment","biotype"), sep = "\\.(?=(all.biotypes|protein.coding))", remove = F, fill = "warn") %>%
  separate(test      , into = c("test","bins")         , sep = "\\.\\."                              , remove = F)   %>%
  mutate(network.type = ifelse(grepl("caudate|dlpfc|hippo|dentate",network),"our","published"),
         new.network        = new.names[network],
         new.network.module = paste0(new.network,">",modules),
         biotype            = replace_na(biotype,"all.biotypes"))



my_new_theme = theme(axis.title        = element_text(face = "bold", size = 12),
                     axis.title.x      = element_text(margin = margin(t = 3,unit = "pt")),
                     axis.title.y      = element_text(margin = margin(r = 5,unit = "pt")),
                     axis.text.y         = element_text(size = 11),
                     axis.text.x       = element_blank(), #element_text(angle = -45, hjust = 0, vjust = 0.5),
                     legend.position   = "none",
                     plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
                     strip.text.y      = element_text(size = 8, angle = 0), #Text for facet_wrap labels
                     # panel.background  = element_blank()
                     panel.background  = element_rect(fill = "white", colour = NA), 
                     # panel.border      = element_rect(fill = NA, colour="grey50"), 
                     panel.grid.major  = element_line(colour = alpha("grey90",0.3), size = 0.2),
                     panel.grid.minor  = element_line(colour = alpha("white",0.3) , size = 0.5),
                     #panel.grid.minor  = element_line(colour = "grey98", size = 0.5),
                     #panel.spacing.x   = unit(2, "pt"),
                     #panel.spacing.y   = unit(2, "pt"),
                     plot.margin        = unit(c(0.5,0.5,0.5,0.5), "pt")
) 

pos  <- position_jitter(seed = 1, height = 0  , width = 0.45)

for (type in c("our","published")){
  for (btype in c("all.biotypes","protein.coding")){
    dat = all.files_long %>% filter(enrichment %in% "PGC3"  &
                                      biotype    %in% btype &
                                      test       %in% c("bonf","hits")  & 
                                      network.type == type    &
                                      !modules   %in% "grey") %>%
      pivot_wider(id_cols = -test, values_from = val, names_from = test)
    
    gg = dat                                                      %>% 
      mutate(final.bonf = -log10(bonf), 
             sig.module = ifelse(final.bonf >-log10(0.05),1,0), 
             bins = factor(bins, levels = our_bins))                    %>%
      ggplot(aes(x = factor(""), y = final.bonf)) +
      geom_jitter(aes(size = ifelse(sig.module, hits,NA), color = I(modules)), pch = "circle", alpha = 0.7, inherit.aes = T, position = pos)+
      geom_text(aes(label = ifelse(sig.module, modules,NA)), position = pos, size = rel(2.5), alpha = 0.9, vjust = 0.5, hjust = "inward") +
      facet_grid(new.network ~ bins, space = "free") +
      geom_hline(aes(yintercept = -log10(0.05)), color= "black", size=0.5,linetype=2) +
      #scale_size_manual(values = c(1,5))+
      my_new_theme+
      xlab(label = "")+
      ylab(label = "-log10(bonf)")+
      ggtitle(label = paste0("PGC ", btype))
    
    print(gg)
  }
}


SCZ_risk_modules_df = all.files_long %>% filter(test %in% "bonf" & !network %in%  c("dentate","dentate.QSVAremoved"))                 %>% 
  mutate(final.val = -log10(val), 
         sig.module = ifelse(final.val >-log10(0.05),1,0), 
         bins = factor(bins, levels = our_bins_complete), new.network = new.names[network])                         %>%
  group_by(enrichment,biotype,new.network.module)                                                                   %>%
  summarise(sig.bins = list(as.character(bins)[sig.module==1]),
            n.sig.bins = lengths(sig.bins))                                                                         %>%
  ungroup()                                                                                                         %>%
  pivot_wider(id_cols = c(enrichment,new.network.module), values_from = c(sig.bins, n.sig.bins), names_from = c(enrichment,biotype)) %>%
  separate(new.network.module, into = c("new.network","module"), sep = ">", remove = F)                             %>%
  mutate(across(contains("n.sig.bins_PGC3"),~replace_na(.x,0)),
         isSCZriskmodules   = ifelse(new.network %in% c("HP.QSVA","HP.noQSVA","DG.noQSVA"),  
                                     ifelse(n.sig.bins_PGC3_all.biotypes >1|n.sig.bins_PGC3_protein.coding >1,1,0),          ##For sample-matched hippo-dentate, use threshold of 2 significant windows instead of 3
                                     ifelse(n.sig.bins_PGC3_all.biotypes >2|n.sig.bins_PGC3_protein.coding >2,1,0)))  %>%
  group_by(new.network)                                                                                                 %>%
  mutate(is_nonGrey_SCZriskmodule_inNetwork = sum(isSCZriskmodules[module != "grey"]) > 0)                         %>%    ##Only keep networks with atleast 1 non-grey SCZ risk module 
  ungroup()                                                                                                         %>%  
  filter(isSCZriskmodules & is_nonGrey_SCZriskmodule_inNetwork)


saveRDS(SCZ_risk_modules_df , file =  "SCZ.risk.modules.df.rds")
saveRDS(SCZ_risk_modules_df %>% pull(new.network.module), file =  "SCZ.risk.modules.rds")


##### Code for sankey plots -----
sankey = all.files %>% dplyr::select(network, enrichment, modules, module.length, hitgenes..kbp_200) %>% 
  mutate(new.network        = new.names[network],
         new.network.module = paste0(new.network,">",modules),
         hitgenes..kbp_200  = hitgenes..kbp_200 %>% set_names(new.network.module))                   %>%
  filter(enrichment %in% "PGC3.all.biotypes" & grepl(" ",new.network))                               %>%
  split(.,.$new.network)

sankey.hitgenes    = map(sankey,"hitgenes..kbp_200") %>% flatten

sankey.pairs = list(
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
    mutate(module2 = ifelse(node %in% SCZ_risk_modules_df$new.network.module, toupper(module),"")) %>%
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
  
  
  sankey1 = htmlwidgets::onRender(sankey1, 'function(el){
                                alltext   = el.querySelector("g.x.axis");
                                alltext.setAttribute("font-size","19");
                                    }')
  
  print(sankey1)
})
  
rm(sankey, sankey.hitgenes, sankey.pairs, sankey.nodes.list) 

##################################################################################
####### Prepare files for SCZ Risk Modules enrichment visualization ##############

##Read list of SCZ risk modules:
SCZ_risk_modules_df = readRDS(file =  "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\SCZ.risk.modules.df.rds")

##load preprocessed results for (PGC3 PERM.ALL, PCG3 LOCI.WINDOWS,MAGMA & HMAGMA)
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\preprocessed_results(some).RData", other.enrichments <- new.env())

## Prepare PGC HYPER.ALL, PGC3 HYPER.PC, PGC3 HYPER.NEG.ALL, PGC3 HYPER.NEG.PC
dat.block1 = SCZ_risk_modules_df %>% tibble::column_to_rownames("new.network.module") %>% 
  dplyr::select(n.sig.bins_PGC3_all.biotypes,n.sig.bins_PGC3_protein.coding,n.sig.bins_PGC3.negative_all.biotypes,n.sig.bins_PGC3.negative_protein.coding) %>%
  set_colnames(c("PGC3 HYPER.ALL", "PGC3 HYPER.PC","PGC3 HYPER.NEG.ALL", "PGC3 HYPER.NEG.PC"))

## Prepare PGC.PERM.ALL and PGC3.LOCI.WINDOWS
dat.block2 = cbind.data.frame(
  PGC.PERM.ALL = other.enrichments$top.gwas.df %>% set_rownames(gsub("\n"," ",rownames(.)))      %>%
    .[SCZ_risk_modules_df$new.network.module,]                                                     %>%
    dplyr::select(sig.wins)                                                                        %>%
    rename(PGC.PERM.ALL = sig.wins)                                                                 ,
  
  PGC.LOCI.WINDOWS = other.enrichments$snp.based.df %>% set_rownames(gsub("\n"," ",rownames(.))) %>%
    .[SCZ_risk_modules_df$new.network.module,]                                                   %>%
    dplyr::select(sig.wins)                                                                      %>%
    rename(PGC.LOCI.WINDOWS = sig.wins)
)                                                                                                %>%
  set_colnames(c("PGC3 PERM.ALL","PGC3 LOCI.WINDOWS"))


##Prepare HMAGMA enrichment results:
HMAGMA.df = other.enrichments$HMAGMA.df %>% set_rownames(gsub("\n"," ",rownames(.)))             %>% 
  .[SCZ_risk_modules_df$new.network.module,]                                                     %>%
  dplyr::select(Adult_brain.ALL)                                                                 %>%
  tibble::rownames_to_column("new.network.module")                                               %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) 


##Prepare MAGMA enrichment results:
MAGMA.df = other.enrichments$MAGMA.df %>% set_rownames(gsub("\n"," ",rownames(.)))               %>% 
  .[SCZ_risk_modules_df$new.network.module,]                                                     %>%
  dplyr::select(`35.10kbp.ALL`,`35.10kbp.PC`)                                                    %>%
  tibble::rownames_to_column("new.network.module")                                               %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) 


##Prepare TWAS SCZ results:
TWAS.SCZ.df = all.files_long %>% filter(enrichment %in% c("TWAS") & test %in% "bonf" & new.network.module %in% SCZ_risk_modules_df$new.network.module) %>% 
  split(.,.$new.network.module) %>%
  imap_dfr(.id = "new.network.module",~{
    col = case_when(
      grepl("CN"    ,.y) ~ c("Paquola_CAUDATE"),
      grepl("HP"    ,.y) ~ c("Jaffe_HIPPO_pgc2.clozuk"),
      grepl("DG"    ,.y) ~ c("Jaffe_HIPPO_pgc2.clozuk"),
      grepl("DLPFC" ,.y) ~ c("Gandal", "Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Hall"),
      TRUE               ~ c("Gandal", "Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Hall")
    )
  
    tmp = .x %>% filter(bins %in% col) %>% group_by(new.network.module) %>%
      summarise("TWAS SCZ" = ifelse(lengths(val) > 1, -log10(metap::sumlog(val)$p), -log10(val))) %>%
      slice_head()
    
  }) %>% 
#  tibble::rownames_to_column("new.network.module")                                               %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) 


## Prepare DEGs.SCZ
DEGs.SCZ.df = all.files_long %>% filter(enrichment %in% c("DEGs") & test %in% "bonf" & new.network.module %in% SCZ_risk_modules_df$new.network.module) %>% 
  split(.,.$new.network.module) %>%
  imap_dfr(.id = "new.network.module",~{
    col = case_when(
      grepl("CN"    ,.y) ~ "Apua_CAUDATE_sczd",
      grepl("HP"    ,.y) ~ "Jaffe_HIPPO_sczd" ,
      grepl("DG"    ,.y) ~ "Jaffe_HIPPO_sczd" ,
      grepl("DLPFC" ,.y) ~ "Jaffe_DLPFC_sczd" ,
      TRUE               ~ "Jaffe_DLPFC_sczd"
    )
    
    tmp = .x %>% filter(bins %in% col) %>% group_by(new.network.module) %>%
      summarise("DEGS SCZ" = ifelse(lengths(val) > 1, -log10(metap::sumlog(val)$p), -log10(val))) %>%
      slice_head()
    
  }) %>% 
  #  tibble::rownames_to_column("new.network.module")                                               %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) 


## Prepare DMGs.SCZ
DMGs.SCZ.df = all.files_long %>% filter(enrichment %in% c("DMGs") & test %in% "bonf" & new.network.module %in% SCZ_risk_modules_df$new.network.module) %>% 
  group_by(new.network.module)                                                         %>%
  summarise("DMGS SCZ" = -log10(metap::sumlog(val)$p))                                  %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) %>%
  dplyr::select(new.network.module,`DMGS SCZ`)


## Prepare DEGs.HUMAN.APES
DEGs.HUMAN.APES.df = all.files_long %>% filter(enrichment %in% c("DEGs") & test %in% "bonf" & new.network.module %in% SCZ_risk_modules_df$new.network.module & bins %in% c("DEGs_human_Sousa")) %>% 
  mutate("DEGS HUMAN.APES" = -log10(val))                                                                   %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.)              %>%
  dplyr::select(new.network.module,`DEGS HUMAN.APES`)


## Prepare LOF
LOF.df = all.files_long %>% filter(enrichment %in% c("pleiotropic_Lof_Denovo_CNVs") & test %in% "bonf" & new.network.module %in% SCZ_risk_modules_df$new.network.module & bins %in% c("LoF")) %>% 
  mutate(LOF = -log10(val))                                                      %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) %>%
  dplyr::select(new.network.module,LOF)


## Prepare DRUG
DRUG.df = all.files_long %>% filter(enrichment %in% c("Druggable_genes") & test %in% "bonf" & new.network.module %in% SCZ_risk_modules_df$new.network.module) %>% 
  group_by(new.network.module)                                                         %>%
  summarise(DRUG = -log10(metap::sumlog(val)$p))                                  %>%
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) %>%
  dplyr::select(new.network.module,DRUG)


## Prepare Cell Specificity enrichment:
cell.specificity.df= res_AB_CellSpecificity %>% map_dfr(.id = "network",~ as.data.frame(.x) %>% dplyr::select(contains("bonf..")) %>% tibble::rownames_to_column("modules")) %>% 
  mutate(new.network = new.names[network],
         new.network.module = paste0(new.network,">",modules),
         across(contains("bonf.."),~ -log10(.)))                                       %>%
  filter(new.network.module %in% SCZ_risk_modules_df$new.network.module)               %>%
  dplyr::select(-network, -new.network, -modules)                                      %>%   
  left_join(data.frame(new.network.module = SCZ_risk_modules_df$new.network.module),.) %>%
  set_colnames(gsub("bonf..","",colnames(.)))


## Prepare Pathology enrichment
Patho.df = SCZ_risk_modules_df %>% set_colnames(gsub("n.sig.bins_|_all.biotypes","",colnames(.))) %>%
  dplyr::select(new.network.module, any_of(patho_bins))


## Prepare GO enrichment:
GO.df= res_AB_GO %>% map_dfr(.id = "network", ~ {.x %>% map_dfr(.id = "GO",~ .x@compareClusterResult)}) %>% 
  rename(module = Cluster)                       %>%
  mutate(GO = gsub("result_all_modules_","",GO),
         new.network = new.names[network],
         new.network.module = paste0(new.network,">",module),
         bonf.val = -log10(p.adjust))            %>%
  filter(p.adjust < 0.05 & new.network.module %in% SCZ_risk_modules_df$new.network.module & (!module %in% "grey"))  %>%   {
    tmp = .
    IDs = tmp %>% group_by(GO,new.network.module)    %>%
      slice_min(p.adjust,n = 3)                      %>%
      ungroup()                                      %>%
      pull(Description)                              %>%
      unique(.)
    
    tmp %>% unite("newID",GO,ID,Description,sep = ">",remove = F) %>%filter(Description %in% IDs)    %>% 
      pivot_wider(id_cols = newID, values_from = bonf.val, names_from = new.network.module)       %>%
      arrange(newID)                                                                              %>%
      tibble::column_to_rownames("newID")                                                         %>%
      t(.)                                                                                        %>%
      as.data.frame(.)
  }
  
  
## Prepare TF enrichment:
TF.df= other.enrichments$TF.df %>% pull(data)    %>% bind_rows(.id = "new.network.module") %>% 
  mutate(new.network.module = gsub("\n"," ",new.network.module),
         bonf.val = -log10(p.adjust))                                                      %>%
  filter(p.adjust < 0.05 & new.network.module %in% SCZ_risk_modules_df$new.network.module) %>% {
    tmp = .
    IDs = tmp %>% group_by(new.network.module)       %>%
      slice_min(p.adjust,n = 3)                      %>%
      ungroup()                                      %>%
      pull(Description)                              %>%
      unique(.)
      
    tmp %>% filter(Description %in% IDs)                                                          %>% 
      pivot_wider(id_cols = Description, values_from = bonf.val, names_from = new.network.module) %>%
      arrange(Description)                                                                        %>%
      tibble::column_to_rownames("Description")                                                   %>%
      t(.)                                                                                        %>%
      as.data.frame(.)
  }   
              

  
save(dat.block1, dat.block2, 
     HMAGMA.df, MAGMA.df, 
     TWAS.SCZ.df, DEGs.SCZ.df, DMGs.SCZ.df, DEGs.HUMAN.APES.df, LOF.df, DRUG.df,
     cell.specificity.df, 
     Patho.df, GO.df, TF.df,
     file = "all.figure.df.RData")
