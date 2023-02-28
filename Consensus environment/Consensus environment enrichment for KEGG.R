####################################################################################################
#                    Script to run various enrichments on gene_module_list:                        #
#                   SCZ, DEGs, DMGs, TWAS, GO, MAGMA, Cell Specificity etc.                        #
#                                                                                                  #
#                  Script to make sankey plots for DLPFC, HP and CN SCZ risk genes.                #
#                                                                                                  #
#                  Script to make binwise PGC3 enrichment scatterplots for networks.               #
#                                                                                                  #
#                          Prepare and format data for other visualisations                        #
#                                "Consensus Environment " plots                                    #
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


con_symbol = c(
  "CHRNA4"
  ,"DENND1A"
  ,"ADCY9"
  ,"ADCY5"
  ,"MKL1"
  ,"DLGAP2"
  ,"RTL10"
  ,"HECW1"
  ,"RAB3A"
  ,"ATP2A2"
  ,"UBE2O"
  ,"KCNQ3"
  ,"NOS1AP"
  ,"KIF3C"
  ,"DPP6"
  ,"SRRM4"
  ,"KSR2"
  ,"SLITRK1"
  ,"OPCML"
  ,"PI4KA"
  ,"KIAA1549L"
  ,"UNC79"
  ,"CLSTN3"
  ,"UNC80"
  ,"TENM2"
  ,"FAM135B"
  ,"PTPRN2"
  ,"RGAG4"
)
con_ensembl = c(
  "ENSG00000101204"
  ,"ENSG00000119522"
  ,"ENSG00000162104"
  ,"ENSG00000173175"
  ,"ENSG00000196588"
  ,"ENSG00000198010"
  ,"ENSG00000215012"
  ,"ENSG00000002746"
  ,"ENSG00000105649"
  ,"ENSG00000174437"
  ,"ENSG00000175931"
  ,"ENSG00000184156"
  ,"ENSG00000198929"
  ,"ENSG00000084731"
  ,"ENSG00000130226"
  ,"ENSG00000139767"
  ,"ENSG00000171435"
  ,"ENSG00000178235"
  ,"ENSG00000183715"
  ,"ENSG00000241973"
  ,"ENSG00000110427"
  ,"ENSG00000133958"
  ,"ENSG00000139182"
  ,"ENSG00000144406"
  ,"ENSG00000145934"
  ,"ENSG00000147724"
  ,"ENSG00000155093"
  ,"ENSG00000242732"
)
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





library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(tidyverse)
library(DOSE)


### Load gene modules (gm_AB) and GWAS bin lists (li_final)
load("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/enrichment.shorter.RData")

### Load GO environments
load("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/GO_DATA_ENSEMBL.RData",envir = .GlobalEnv)

module.color0 = WGCNA::standardColors(length(unique(unlist(map_depth(Hartl2021,1,names))))) %>% set_names(unique(unlist(map_depth(Hartl2021,1,names))))
module.color0["grey"] = "grey"
names(gm_AB)[grep("Gandal2018$"     ,names(gm_AB))] = "Gandal2018a"
names(gm_AB)[grep("Gandal2018PE$"   ,names(gm_AB))] = "Gandal2018b"
names(gm_AB)[grep("Gandal2018PE_cs$",names(gm_AB))] = "Gandal2018b_cs"

Hartl2021_gm <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/Hartl2021-gene-module-list.rds")
gm_AB$Hartl2021_BRNCTX =Hartl2021_gm$BRNCTX

### Convert gene module lists to gene symbols
allgenes = unique(unlist(gm_AB))
mapping  = clusterProfiler::bitr(allgenes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

gene_grch38 = geneMap_fun(genes = allgenes, ensembl = "dec2017.archive.ensembl.org",
                          attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype", "transcript_count","percentage_gene_gc_content","entrezgene","version"))
gene_grch38 = gene_grch38[!duplicated(gene_grch38$ensembl_gene_id),]
gene_grch38$ensembl    = gene_grch38$ensemblID    = gene_grch38$ensembl_gene_id
gene_grch38$hgnc       = gene_grch38$Symbol       = gene_grch38$external_gene_name
gene_grch38$gene_type = gene_grch38$gene_biotype

gm_PC = map(gm_AB,~ {
  tmp = map(.x, ~{intersect(.x,gene_grch38$ensembl_gene_id[gene_grch38$gene_biotype == "protein_coding"])})
  tmp = tmp[lengths(tmp)>20]  ##Discarding modules with less than 21 genes when restricting to only protein.coding genes
})

gm_AB_sym = map(gm_AB,~ {
  tmp = map(.x, ~{gene_grch38$external_gene_name[gene_grch38$ensembl_gene_id %in% .x]})
  tmp = tmp[lengths(tmp)>20]  ##Discarding modules with less than 21 genes when restricting to only protein.coding genes
})

gm_PC_sym = map(gm_PC,~ {
  tmp = map(.x, ~{gene_grch38$external_gene_name[gene_grch38$ensembl_gene_id %in% .x]})
  tmp = tmp[lengths(tmp)>20]  ##Discarding modules with less than 21 genes when restricting to only protein.coding genes
})

### Identify genes in modules with consensus genes
geneswithcon = function(gm, con) {
  map(gm, ~ {
    net = ..1 
    con_ind = map_lgl(net, ~ {
      if(sum(con %in% ..1)>0) {
        out = TRUE
      } else {out = FALSE}
    })
    # print(con_ind)
    not_con_ind = map_lgl(net, ~ {
      if(sum(con %in% ..1)>0) {
        out = FALSE
      } else {out = TRUE}
    })
    
    con_mod_genes = net[con_ind]
    con_mod_genes = unlist(con_mod_genes[!names(con_mod_genes) == "grey"])
    
    not_con_mod_genes = net[not_con_ind]
    not_con_mod_genes = unlist(not_con_mod_genes[!names(not_con_mod_genes) == "grey"])
    
    grey_genes = net[["grey"]]
    
    list(Genes_With_Consensus = con_mod_genes,
         Genes_Not_With_Consensus = not_con_mod_genes,
         grey = grey_genes)
  })
}

gm_PC_con = geneswithcon(gm_PC, con_ensembl)
gm_AB_con = geneswithcon(gm_AB, con_ensembl)

gm_PC_con = geneswithcon(gm_PC_sym, con_symbol)
gm_AB_con = geneswithcon(gm_AB_sym, con_symbol)


### Compute KEGG enrichment
kegg_sets_sym = readRDS("../data/Ontology sets/kegg_genelists_nonest.rds")
kegg_sets_sym_new = map(kegg_sets_sym , ~ {
  intersect(..1, unlist())
})
allkegggenes_sym = unique(unlist(kegg_sets_sym))
n_allkegggenes_sym =length(allkegggenes_sym)

gm_AB_con_kegg = imap_dfr(gm_AB_con, ~ {
  # print(..2)
  net_name = ..2
  net = ..1
  con_set = net$Genes_With_Consensus
  # n_con_set = length(con_set)
  universe = c(net$Genes_With_Consensus, net$Genes_Not_With_Consensus)
  n_universe = length(universe)
  print(n_universe)
  # n_universe = length(universe)
  
  kegg_sets_sym_new = map(kegg_sets_sym , ~ {
    intersect(universe,..1)
  })
  
  universe = intersect(universe,unlist(kegg_sets_sym_new))
  n_universe = length(universe)
  # print(n_universe)
  
  con_set = intersect(con_set,unlist(kegg_sets_sym_new))
  n_con_set = length(con_set)
  
  out = imap_dfr(kegg_sets_sym_new, ~ {
    title = ..2
    gs = ..1
    n_gs = length(gs)
    # print(gs)
    con_in_gs = intersect(con_set, gs)
    n_con_in_gs = length(con_in_gs)
    
    # print(n_con_in_gs)
    
    
    out = if(n_con_in_gs > -1) {
      1.0-phyper(n_con_in_gs-1, n_gs, n_universe-n_gs, n_con_set)
    } else {1}
    
    c(network = net_name,
      name = title,
      genes  = paste(con_in_gs, collapse = ","),
      n_universe = n_universe,
      n_genes_con = n_con_set,
      n_genes_ont = n_gs,
      n_genes_con_in_ont = n_con_in_gs,
      p = as.numeric(out)
    )
  }) %>% mutate_at(c('p'), as.numeric)
  
  out$p.adjust = p.adjust(out$p, method = "BH")
  out
})


### Format for plotting
new.names = read_csv("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/madhur data/old-to-new Network names all.csv") %>% tibble::deframe()
new.names[grepl("Gandal2018$"     ,names(new.names))] = "Gandal2018a"
new.names[grepl("Gandal2018PE$"   ,names(new.names))] = "Gandal2018b"
new.names[grepl("Gandal2018PE_cs$",names(new.names))] = "Gandal2018b_cs"
new.names = c(new.names,c("Gandal2018a"="Gandal2018a","Gandal2018b"="Gandal2018b","Gandal2018b_cs"="Gandal2018b_cs"))
new.names = c(new.names,c("Gandal2018a"="Gandal2018a","Gandal2018b"="Gandal2018b","Gandal2018b_cs"="Gandal2018b_cs"))
new.names = c(new.names,c("Hartl2021_BRNCTX"="Hartl2021_BRNCTX"))

KEGG.df.og= as.data.frame(gm_AB_con_kegg) %>% 
  mutate(new.network = new.names[network],
         # KEGG = gsub(".*00","",name),
         # new.network.module = paste0(new.network,">",module),
         bonf.val = -log10(p.adjust+.00000000000001)
         )            %>%
  filter(p.adjust < 0.05& 
           # !new.network %in% c("DLPFC Perinatal"  ,  "DLPFC Juvenile"   ,  "DLPFC Adult"      ,  "DLPFC Older Adult") &
           !grepl("CN|HP|DG|2018b",new.network) )  %>%   {
             tmp = .
             IDs = tmp %>% group_by(name,network)    %>%
               slice_min(p.adjust,n =5)                      %>%
               ungroup()                                      %>%
               # pull(Description)                              %>%
               unique(.)
             
             tmp %>% unite("newID",name,sep = "",remove = F)   %>%
               pivot_wider(id_cols = newID, values_from = bonf.val, names_from = new.network)       %>%
               arrange(newID)                                                                              %>%
               tibble::column_to_rownames("newID")                                                         %>%
               t(.)                                                                                        %>%
               as.data.frame(.)
           }

KEGG.df.og = rbind(KEGG.df.og, Pergola2017 = 0)
KEGG.df = KEGG.df.og[!rownames(KEGG.df.og) == "DLPFC",]%>% filter(if_any(everything(), ~ !is.na(.)))
KEGG.df = KEGG.df[match(c(
  "DLPFC Perinatal"   ,  "Werling2020"   ,     "Walker2019" ,     "Li2018"        , 
  "DLPFC Juvenile" ,    "DLPFC Adult"    ,    "DLPFC Older Adult" ,       "Hartl2021_BRNCTX",
  "Gandal2018a"   ,"Fromer2016_case" ,   "Fromer2016_control", 
  "Pergola2017" ,
  "Pergola2020" ,       "Pergola2019"     ,   "Radulescu2020"     
  
),rownames(KEGG.df) ),]
KEGG.df = KEGG.df[,colSums(is.na(KEGG.df))<nrow(KEGG.df)]
KEGG.df = as.data.frame(t(KEGG.df))
KEGG.df[is.na(KEGG.df)] <- 0
rownames(KEGG.df) = gsub(".*>","", rownames(KEGG.df))
KEGG.df$median =  apply(KEGG.df, 1, median, na.rm=T)
KEGG.df = KEGG.df[order(KEGG.df$median, decreasing = TRUE),]
KEGG.df = KEGG.df[1:(ncol(KEGG.df)-1)]
KEGG.df = t(KEGG.df)

zero_th = 7
KEGG.df = KEGG.df[, colSums(KEGG.df != 0) > zero_th]
newnames = colnames(KEGG.df) 
newnames = gsub("\\[.*", "", newnames)
newnames = gsub('.*[[:digit:]]+', "", newnames)
colnames(KEGG.df) = newnames


library(corrplot)
library(RColorBrewer)

x11();corrplot((as.matrix(KEGG.df[])), method = 'color'
               ,is.corr = FALSE
               # , p.mat = 10^-t(KEGG.df[1:4])
               # , order = 'alphabet'
               , tl.col = 'black',
               cl.cex = .7,
               tl.cex = 1.2,
               # cl.pos = 'n',
               # cl.ratio = 0,
               tl.srt = 90,  col=colorRampPalette(c("white","darkcyan"))(200))

