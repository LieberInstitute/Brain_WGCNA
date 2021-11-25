Enrich = function(ll.target, ll.components, ll.components.size, network_name, target_name){
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


options(java.parameters = "-Xmx8000m")  #Required for xlsx package
library(SummarizedExperiment)
library(purrr)
set.seed(123)
library(limma)
library(dplyr)
library(pheatmap)
library(gtools)       #Optional. Just to sort path of files read
library(furrr)
options(mc.cores = 3, future.globals.maxSize= 3565158400)

library(org.Hs.eg.db)
library(qs)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(directlabels)
library(patchwork)
library(tidyr)
library(cowplot)
library(ggtext)

our_bins = c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")
new_bins = c(our_bins,"MAGMA","consensus")
our_networks = c("caudate","caudate__.1.25","caudate__25.50","caudate__50.100",
                 "dentate",
                 "dentate.noQSVAremoved","dentate.QSVAremoved",
                 "hippo.noQSVAremoved","hippo.QSVAremoved",
                 "dlpfc","dlpfc__.1.6","dlpfc__6.25","dlpfc__25.50","dlpfc__50.100",
                 "hippo","hippo__.1.6","hippo__6.25","hippo__25.50","hippo__50.100",
                 "rse.pwr7.protectDx","stemcell.protectDx.pwr7",
                 "rse.pwr7","stemcell.pwr7",
                 "stemcell.noConsensus.pwr7",
                 "stemcell.protectDx.noConsensus.pwr7",
                 "simnoKO.pwr7","sim0.pwr7"
                 )            

ne = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\PGC.gene.lists_grch38[125 PGC ensemblIDs].RData", envir = ne)
li_consensus = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/SNPextension/consensus_genes/consensus_gene_list.rds")

li_new = list(positive_all = c(li_consensus$SCZ_prenatal_positive,li_consensus$SCZ_postnatal_positive),
              negative_all = c(li_consensus$SCZ_prenatal_negative,li_consensus$SCZ_postnatal_negative))
li_new$consensus_noPGC = li_new$positive_all[!li_new$positive_all %in% unlist(ne$PGC3.all.biotypes[our_bins])]
li_MAGMA = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/MAGMA_consensus_gene_list.rds")
li_MAGMA_noPGC = li_MAGMA[!li_MAGMA %in% unlist(ne$PGC3.all.biotypes[our_bins])]

li_final = c(ne$PGC3.all.biotypes, consensus = list(li_new$consensus_noPGC), MAGMA = list(li_MAGMA_noPGC))

##Mega gene list: All networks + stem cell protectDx + stem cell protectDx shuffled
gm <- c(readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (wide_form_test) (all networks)[grch38].rds"),
        readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/redone/results/gene-module list(stem cells)_protectDx(redone).rds"),
        readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/redone/results/KO sft estimation shuffeled/KO nets/gene-module list(stem cells)_protectDx_KOsims(redone).rds")
        )
gm                               = gm[!grepl("pwr15",names(gm))]
names(gm)                        = gsub("\\(redone\\)","",names(gm))
names(gm)[grep("rse",names(gm))] = gsub("rse","stemcells",names(gm)[grep("rse",names(gm))])
names(gm)[grep("sim",names(gm))] = gsub("...sft.pearson.signed hybrid & pwr-7 & pwr.type-by_individual_protectDx",".pwr7" ,names(gm)[grep("sim",names(gm))])
names(gm)[grep("sim",names(gm))] = gsub("...sft.pearson.signed hybrid & pwr-15 & pwr.type-by_all_protectDx"      ,".pwr15",names(gm)[grep("sim",names(gm))])
names(gm)[grep("sim",names(gm))] = gsub("sim","stemcells.protectDx.sim",names(gm)[grep("sim",names(gm))])


allgenes = unique(unlist(gm))
mapping = clusterProfiler::bitr(allgenes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

####################################
#######################################
##All network modules  enrichment (hypergeo) for PGC3+consensus_list
res = imap(gm,~{
  #if (.y == "dlpfc__25.50") browser()
  tmp = Enrich(ll.target = li_final, ll.components = .x, network_name = names(.x), target_name = names(li_new))
})


############
all.files_SCZ = imap_dfr(.id = "network",res,~{
  
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

all.files_SCZ = all.files_SCZ %>% mutate(network = case_when(
  grepl("stemcells.pwr7.protectDx"     ,network)  ~ "simnoKO.pwr7" ,
  grepl("stemcells.protectDx.sim.*pwr7",network)  ~ gsub("stemcells.protectDx.","",network),
  TRUE                                            ~ network
))

# saveRDS(all.files_SCZ,"C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Stem_Cell_Project_2021\\redone\\results\\all.files_SCZ(redone).rds")
# all.files_SCZ = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Stem_Cell_Project_2021\\redone\\results\\all.files_SCZ(redone).rds")


stem_cell_results  = all.files_SCZ[grep("sim",all.files_SCZ$network),]
stem_cell_results  = stem_cell_results[,!grepl("hitgenes",colnames(stem_cell_results))]


stem_cell_results_long = stem_cell_results %>% filter(grepl("sim",network)) %>% 
  pivot_longer(cols = c(contains("bonf.."),contains("pvals.."),contains("zvals.."),contains("Fold.hits.."),starts_with("hits..")),names_to = "type") %>% 
  separate(type, into = c("enrichment","bins"), sep = "\\.\\.", remove = F) %>% 
  mutate(bins = factor(bins,levels = our_bins),
         pwr = case_when(grepl("pwr7" ,network)~"pwr7",
                         grepl("pwr15",network)~"pwr15",
                         TRUE                  ~ "NA"),
         new_network = case_when(
           grepl("noKO" ,network) ~ "noKO",
           grepl("sim0.",network) ~ "sim0",
           TRUE                   ~ "simALL"
         )) %>%
  filter(bins %in% new_bins)


##########Plotting enrichment in PGC lists vs enrichment in consensus list
gg1 = list()
cur_network = "simnoKO.pwr7"  ##no KO network
#cur_network = "sim0.pwr7"    ##shuffling expression for only 23 nonGWAS consensus genes
#cur_network = paste0("sim",1:100,".pwr7")   ##shuffling expression for 23 other genes (*100 times)

all_cols = c("network", "modules", "module.length", paste0("bonf..",new_bins),paste0("pvals..",new_bins),paste0("hits..",new_bins),paste0("Fold.hits..",new_bins))

dat_all = all.files_SCZ[,all_cols] %>% filter(network %in% cur_network & modules != "grey") %>% mutate(across(c(contains("pvals.."),contains("bonf..")), ~ -log10(.x)))
dat     = all.files_SCZ[,all_cols] %>% filter(network %in% cur_network & modules != "grey"& hits..consensus >0)  %>% mutate(across(c(contains("pvals.."),contains("bonf..")), ~ -log10(.x)))

sc1 = range(dat[,grepl("pvals.."    ,colnames(dat)) & !grepl("MAGMA",colnames(dat))])


for (bb in our_bins){
  
  ##Prepare table for Fisher test, all modules
  tb = matrix(nrow = 2,data =
                c(
                  sum(dat_all[,"pvals..consensus"]  < 0.5 & dat_all[,paste0("pvals..",bb)]  < 0.5), #Both yes
                  sum(dat_all[,"pvals..consensus"]  < 0.5 & dat_all[,paste0("pvals..",bb)] >= 0.5), #only consensus yes
                  sum(dat_all[,"pvals..consensus"] >= 0.5 & dat_all[,paste0("pvals..",bb)]  < 0.5), #only GWAS yes
                  sum(dat_all[,"pvals..consensus"] >= 0.5 & dat_all[,paste0("pvals..",bb)] >= 0.5)  #Both no
                )
  ) %>% set_rownames(c("GWASenriched:YES","GWASenriched:NO")) %>% set_colnames(c("ConsensusGenes:YES","ConsensusGenes:NO"))

  ff = broom::tidy(fisher.test(tb, alternative = "g"))[,c("estimate","p.value")]                      ##all modules
  rr = broom::tidy(robustbase::lmrob(paste0("pvals..consensus ~pvals..", bb), data = dat))$p.value[2] ##only modules with consensus gene hit
  

  gg1[[bb]] = all.files_SCZ %>% filter(network %in% cur_network & modules != "grey") %>%
    ggplot(aes(x = -log10(pvals..consensus), y = -log10(.data[[paste0("pvals..",bb)]]))) +
    geom_point(data = all.files_SCZ %>% filter(network %in% cur_network & modules != "grey" & hits..consensus > 0), alpha = 0.8, size = 3, aes(fill = I(modules)), pch = 21, stroke = 0.5, show.legend = F)+
    geom_smooth(method = "lm", alpha = 0.2, color = adjustcolor("black",alpha.f = 0.3),fullrange = T)+
    xlim(extendrange(sc1,f = c(0.2,0.1))) + ylim(extendrange(sc1,f = c(0.1,-0.25)))+ #coord_cartesian(xlim = sc1, ylim = sc1, expand = F, clip = "off")
    xlab(label = NULL) + ylab(label = NULL)+
    theme_bw() +
    theme(
      strip.text.x        = element_text(size = rel(1), angle = 90, hjust = 0.5),
      axis.text.x         = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.y         = element_text(size = rel(1.4)),
      plot.title          = element_text(hjust = 0.5, size = rel(1.45))
    )
  if (rr < 0.05 & !is.na(rr)) {
    gg1[[bb]] = gg1[[bb]] + theme(plot.title = element_text(colour = "red"))
  }
  gg1[[bb]] = gg1[[bb]] + labs(title = glue::glue("{bb}:  {round(ff$estimate,2)}[{format.pval(ff$p.value,digits = 3)}],  RLM:  {format.pval(rr, digits = 3)}"))
}
gg1_final = cowplot::plot_grid(plotlist = gg1) + 
  theme(plot.margin = unit(c(0.8,0,0.7,0.7),"cm")) + 
  annotate(geom = "text", label ="-log10(p.vals) for consensus genes (hypergeo)", x = 0.5  , y = -0.01, size = rel(5), angle = 0 ) +
  annotate(geom = "text", label ="-log10(p.vals) for SCZ PGC3 (hypergeo)"       , x = -0.01, y = 0.5  , size = rel(5), angle = 90) +
  annotate(geom = "text", label =paste0(cur_network)     , x = 0.5  , y = 1.025 , size = rel(5.5)  , angle = 0 ) 
print(gg1_final)


