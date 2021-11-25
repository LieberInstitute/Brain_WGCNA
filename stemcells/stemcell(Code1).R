options(java.parameters = "-Xmx8000m")  #Required for xlsx package
library(SummarizedExperiment)           # check this link: https://f1000research.com/articles/6-1558/v1 , https://www.bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
library(recount)                        # check this paper https://www.nature.com/articles/nbt.3838
library(jaffelab)                       # vignette http://127.0.0.1:13522/library/recount/doc/recount-quickstart.html
library(limma)
library(dendextend)
library(purrr)
library(dplyr)

library(pheatmap)
library(RNOmni)
library(quantro)
library(BRETIGEA)
library(ggplot2)
library(gghalves)
library(tidyr)
library(magrittr)

set.seed(123)

setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/")

#Read list of sample IDs for the neuronal samples
nn = readRDS("redone/neuronal_samples.rds")

getRSEdata = function(){
  rse = list()
  load("data/rse_gene_humanRat_crossSpecies_annotatedAndMerged_allData_n304.Rdata")
  rse[["ercc"]] = rse_ercc
  rse[["gene"]] = rse_gene
  
  return(rse)
}


#Load all RSE data in a list
rse_raw_list = getRSEdata()
rse_raw_list_updated = list(gene = rse_raw_list$gene)


#Pre-process rse-data
rse_preprocess_step1 = rse_raw_list_updated %>% purrr::imap(~{
  rse_merged = .
  #rse_merged = merge_rse_metrics(.)
  assays(rse_merged)$RPKM     = recount::getRPKM(rse_merged, length_var = "Length", mapped_var = "numMapped")
  assays(rse_merged)$logRPKM  = log2(assays(rse_merged)$RPKM+1)
  return(rse_merged)
})


#Subset RSE objects for desired Samples: : All male, all European ancestry 
rse_preprocess_step2 = rse_preprocess_step1 %>% purrr::imap(~{
  rse_merged_small = SummarizedExperiment::subset(
    x      = .                                     ,
    subset = !grepl("ENSRNOG",gencodeID)           , #Keep only human genes (remove rat genes)          
    select = SAMPLE_ID %in% nn)                      #Keep only neuronal DIV samples
  return(rse_merged_small)
})


#Outlier sample removal (from each subgroup)
remove_outlier_samples = function(rse, name.rse , sdout = 3){
  
  exp                 = assays(rse,withDimnames = F)$logRPKM
  IAC                 = as.matrix(dist(t(exp)))    #Sample-sample pairwise euclidean distance/disimilarity (based on all genes)
  
  scaled.IAC          = scale(matrixStats::rowMeans2(IAC))
  #scaled.IAC          = scale(matrixStats::rowSums2(IAC)/(ncol(IAC)-1)) #No need to consider distace of a sample from itself
  names(scaled.IAC)   = rownames(IAC)
  samp.outlier        = abs(scaled.IAC) > sdout
  names(samp.outlier) = rownames(IAC)
  
  IACnew              = IAC[!samp.outlier,!samp.outlier] #After removing outliers
  
  dend = IAC %>% as.dist %>% hclust (method = "average") %>% as.dendrogram
  #labels(dend)[!(labels(dend) %in% rownames(IAC)[samp.outlier])] = ""
  cll = rep("black", length(labels(dend))); names(cll) = labels(dend)
  cll[labels(dend) %in% names(samp.outlier)[samp.outlier]] = "red"
  dend = dend %>% set("labels_col", cll)
  dend = dend %>% set("labels_cex", 0.6)
  
  dendnew = IACnew %>% as.dist %>% hclust (method = "average") %>% as.dendrogram
  dendnew = dendnew %>% set("labels_cex", 0.6)
  
  par(mfrow = c(2,2), mar = c(2.5,4,1.5,1),oma = c(1,1,0.5,0.5))
  y.lim = c(0,max(attributes(dend)$height,attributes(dendnew)$height))
  plot(dend    ,ylim = y.lim)
  plot(dendnew ,ylim = y.lim)
  plot(scaled.IAC[match(labels(dend),names(scaled.IAC))], ylim = c(-6,6), pch = 19, col = cll, ylab = "z-scored sample")
  abline(h=c(sdout,-sdout), col = "red")
  title(name.rse,outer = T,line = -1)
  
  return(samp.outlier)
}

# Remove outlier samples from age-parsed data using logRPKM data
rse_outlier_removed = rse_preprocess_step2 %>% purrr::imap( ~{
  found.outliers = remove_outlier_samples(.x,.y,sdout = 3)
  rse_merged = SummarizedExperiment::subset(
    x      = .x                                     ,
    select = !found.outliers                        ,
    subset = TRUE
  )
  return(list(rse = rse_merged, removed = names(found.outliers)[found.outliers]))    
})


rse_outlier_removed_big = map(rse_outlier_removed, "rse")
samples_removed         = map(rse_outlier_removed, "removed")

#Outlier gene removal
rse_preprocess_step3 = imap(rse_outlier_removed_big, ~{
  expMat            = assays(.)$RPKM
  median.gene.exprs = matrixStats::rowMedians(expMat)
  names(median.gene.exprs) = rownames(expMat)
  
  
  #Use 'identical(all.equal(x,0),TRUE)' for numerical comparison on floating point numbers
  gene_filter0 = as.vector (         rowSums(expMat ==0) <= ncol(expMat)*0.2) #Genes to keep (not more than 20% samples are zero)
  gene_filter1 = as.vector(          median.gene.exprs   >= .1)               #Genes to keep (median exp > 0.1)
  #gene_filter2 = as.vector(abs(scale(median.gene.exprs)) <  3)            #Genes to keep (z-scoreed_median_exp < 3)
  
  median.gene.exprs_filtered = median.gene.exprs[gene_filter0 & gene_filter1]
  gene_filter2 = as.vector(abs(scale(median.gene.exprs_filtered)) <  3)            #Genes to keep (z-scored_median_exp < 3)
  
  
  #genes_to_keep = unique(rownames(expMat)[gene_filter0 & gene_filter1 & gene_filter2])
  genes_to_keep  = names(median.gene.exprs_filtered)[gene_filter2]
  
  print(.y)
  print(rbind(table(gene_filter0), table(gene_filter1), table(gene_filter2)))
  print(paste0(length(genes_to_keep),"/",dim(expMat)[1]))
  
  rse_merged2 = SummarizedExperiment::subset(
    x      = .                                     ,
    subset = (gencodeID  %in% genes_to_keep)     &             #"subset" argument subsets from rowData(rse)
      (!grepl("^MT-",Symbol)        )                   #Gene symbol starting with 'MT-' (for mitochondrial genes)
  )
  
  
  
  exp_quant_normalised           = preprocessCore::normalize.quantiles(assays(rse_merged2)$logRPKM) #Quantile normalise samples. Ref: http://jtleek.com/genstats/inst/doc/02_05_normalization.html
  dimnames(exp_quant_normalised) = dimnames(assays(rse_merged2)$logRPKM)
  #exp_ranknorm                   = t(apply(exp_quant_normalised,1, RNOmni::rankNorm))  
  
  assays(rse_merged2)$quant.norm.logRPKM = exp_quant_normalised
  return(rse_merged2)
})


rse_neuronal = rse_preprocess_step3$gene
#rse_neuronal = rse_neuronal[,order(rse_neuronal$Line)]
exp          = assays(rse_neuronal)$logRPKM
#exp_RPKM = assays(rse_neuronal)$RPKM
sampleInfo = as.data.frame(colData(rse_neuronal)) %>% tibble::rownames_to_column("sampleID")
all.equal(colnames(exp), sampleInfo$SAMPLE_ID)

rm(rse_raw_list, rse_raw_list_updated, rse_preprocess_step1, rse_preprocess_step2, rse_preprocess_step3, rse_outlier_removed, rse_outlier_removed_big)

#Make longform data of expression + sampleInfo
fulldat = exp %>% as.data.frame %>% tibble::rownames_to_column("ensemblID") %>%  tidyr::pivot_longer(-ensemblID, names_to = "sampleID", values_to = "exp.logRPKM") %>% 
  left_join(sampleInfo,by = "sampleID") %>% mutate(across(any_of(c("DIVgroup", "Line", "RealGenome","Code", "RealCode", "RealDx")), ~factor(.x)))

#Select only specific columns from fulldat
subdat = fulldat %>% dplyr::select(any_of(c("ensemblID", "exp.logRPKM", "DIVgroup", "Line", "RealGenome", "RealCode", "RealDx")))
names(subdat$exp.logRPKM) = paste0(subdat$RealGenome,">",subdat$Line,">",subdat$DIVgroup,">",subdat$ensemblID)


###Get the mean expression by each RealGenome
subdat1 = subdat %>% dplyr::select(c(-RealCode))             %>%
  group_by(ensemblID,RealGenome)                             %>%
  mutate(across(RealDx, ~ .x)                                ,
            n.Lines       = n_distinct(Line)               ,
            n.DIVgroups   = n_distinct(DIVgroup)           ,
            n.all         = exp.logRPKM %>% length         ,
            mean          = exp.logRPKM %>% mean           ,
            var           = exp.logRPKM %>% var            ,
            values        = exp.logRPKM %>% list           ,
            scaled.values = exp.logRPKM %>% scale %>% list ,
            max.scale     = exp.logRPKM %>% scale %>% max  ,
            min.scale     = exp.logRPKM %>% scale %>% min   )                                                          %>%
  slice_head(n=1)                                            %>%
  ungroup()                                                  %>%
  dplyr::select(c(-Line,-DIVgroup,-exp.logRPKM))

saveRDS(subdat1, "redone/subdat1(redone).rds")
subdat1 = readRDS("redone/subdat1(redone).rds")
  

#########################################
##Make new expression assays with the averaged data
exp_logRPKM_new = subdat1 %>% pivot_wider(id_cols = ensemblID, values_from = mean, names_from = RealGenome) %>% tibble::column_to_rownames("ensemblID") %>% as.matrix()

#Make new colData
colData_new = as.data.frame(colData(rse_neuronal)) %>% 
  group_by(RealGenome)                                                                     %>% 
  summarise(across(c(RealDx,PRS, Experiment                                    ), ~.x ),
            across(c(h_mitoRate,combined_mitoRate, totalAssignedGene, rRNA_rate), mean))   %>%
  slice_head()                                                                             %>%
  ungroup()                                                                                %>%
  as.data.frame()                                                                          %>%
  set_rownames(.$RealGenome)
colData_new = colData_new[colnames(exp_logRPKM_new),]

#Make new rowData
rowData_new = rowRanges(rse_neuronal)
rowData_new = rowData_new[rownames(exp_logRPKM_new),]

rse_new = SummarizedExperiment(assays    = list(logRPKM = exp_logRPKM_new),
                               colData   = colData_new            ,
                               rowRanges = rowData_new            )

##Check normality of gene expression among samples 
li_exp = li_exp_new = vector()
for (i in 1:nrow(exp)){
  li_exp[i]     = shapiro.test(as.numeric(    exp[i,]))$p
  li_exp_new[i] = shapiro.test(as.numeric(exp_logRPKM_new[i,]))$p
}

# p > 0.05 support normality
table(li_exp     >0.05)
table(li_exp_new >0.05)


get.cells.proportions = function(rse,rse.name, assay.name, method, nmarkers){
  exp = assays(rse)[[assay.name]]       #Genes in rows
  exp = exp[!duplicated(rowData(rse)$Symbol) & !is.na(rowData(rse)$Symbol),]
  rownames(exp) = rowData(rse)$Symbol[!duplicated(rowData(rse)$Symbol) & !is.na(rowData(rse)$Symbol)]   #Convert ENSEMBL names to GeneSymbol
  zz1 = brainCells(exp, nMarker=nmarkers, species="human", method = method, scale = T)
  zz2 = brainCells(exp, nMarker=nmarkers, species="human", method = method, scale = F)
  colnames(zz2) = paste0(colnames(zz2),"_unscaled")
  zz = cbind(zz1,zz2)
  return(zz)
}



rse_new_preprocess_step4 = list(rse_new = rse_new) %>% purrr::imap(~{
  nmarkers = 50
  cell.prop.info_SVD = get.cells.proportions(.x, .y, assay.name = "logRPKM", method = "SVD", nmarkers)
  colData(.x) = cbind(colData(.x), cell.prop.info_SVD)
  return(.x)
})


setwd("redone\\results")
 saveRDS(rse_new_preprocess_step4      , "rse_preprocessed(redone).rds")         #List of rse object preprocessed (sample and gene outlier removed. Quantile)
