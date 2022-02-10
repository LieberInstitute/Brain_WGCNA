####################################################################################################
#                       Script to Preprocess samples match rse data for HP-DG                      #
#                               "Cell-population enrichment" studies                               #
#                                                                                                  #
####################################################################################################

options(java.parameters = "-Xmx8000m")  # Required for xlsx package
library(SummarizedExperiment)           # check this link: https://f1000research.com/articles/6-1558/v1 , https://www.bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
library(recount)                        # check this paper https://www.nature.com/articles/nbt.3838
library(jaffelab)                       # vignette http://127.0.0.1:13522/library/recount/doc/recount-quickstart.html
library(limma)
library(dendextend)
library(purrr)
library(dplyr)
library(preprocessCore)
library(pheatmap)
library(RNOmni)

set.seed(123)

setwd("C:/Users/Madhur/OneDrive - Johns Hopkins/")

#Read genomic eigenvariables
gnm = read.csv("data_source/genomicEigenVariables LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_Thrice_dropBrains_maf05_hwe6_geno10.csv",stringsAsFactors = F)
gnm = gnm[!duplicated(gnm$IID),]

getRSEdata = function(){
  rse = list()
  load("data_source/NIH_R21/Expression/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
  rse$hippo   = rse_gene
  load("data_source/NIH_R21/Expression/astellas_dg_hg38_rseGene_n263.rda")
  rse$dentate = rse_gene
  
  return(rse)
}

#Load all RSE data in a list
rse_raw_list = getRSEdata()


#Save gene_map/gene_info for later reference:
# gene_map           = as.data.frame(rowData(rse_raw_list$caudate), stringsAsFactors = F)
# gene_map           = gene_map[!duplicated(gene_map$ensemblID),]
# colnames(gene_map)[grep("ensembl",colnames(gene_map))] = "ensembl"
# colnames(gene_map)[grep("Symbol" ,colnames(gene_map))] = "hgnc"
# rownames(gene_map) = gene_map$ensembl
# saveRDS(gene_map, "C:/Users/mpariha1/Desktop/upload_to_google_drive/results/gene_map.rds")


#Combine Genomic Eigen variable information into the rse -colData:
rse_raw_list_updated = map(rse_raw_list,~{
  yy = gnm[match(.$BrNum, gnm$IID), c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")]
  colData(.) = cbind(colData(.),yy)
  return(.)
})


#Pre-process rse-data
rse_merged_list = rse_raw_list_updated %>% purrr::imap(~{
  rse_merged = if (.y !=  "dentate") merge_rse_metrics(.) else (.)
  assays(rse_merged)$RPKM     = recount::getRPKM(rse_merged, length_var = "Length", mapped_var = "numMapped")
  assays(rse_merged)$logRPKM  = log2(assays(rse_merged)$RPKM+1)
  return(rse_merged)
})


#Subset RSE objects for desired Samples:
minAge = -1;maxAge = 150; minRin = 6; diagnosis = c("Control"); race = c("AA","CAUC")
rse_preprocess_step1 = rse_merged_list %>% purrr::imap(~{
  rse_merged_small = SummarizedExperiment::subset(
    x      = .                                     ,
    subset = TRUE                                  ,                         #----"subset" argument subsets from rowData(rse) i.e. genes  ##(default = T, all rows)
    select =                                                                 #----"select" argument subsets from colData(rse) i.e. samples ##(default = T, all columns)
      (Age                >=   minAge       )   &
      (Age                <=   maxAge       )   &
      (Race              %in%  race         )   &
      (Dx                %in%  diagnosis    )   &
      (if (.y == "dentate") {TRUE} else {sapply(RIN, mean)  >=   minRin}) &  #RIN filter only if there are no missing values
      (!is.na(C1)                           )                                #Keep if C1 genomic information is not-missing
  )
  return(rse_merged_small)
})


commonBrNum_hippo_dentate = intersect(rse_preprocess_step1$hippo$BrNum, rse_preprocess_step1$dentate$BrNum)

rse_preprocess_step2 = rse_preprocess_step1 %>% purrr::imap(~{
  rse = SummarizedExperiment::subset(
    x      = .,
    select =(BrNum %in% commonBrNum_hippo_dentate)
  )
  return(rse)
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


rse_outlier_removed_big = map(rse_outlier_removed,"rse")
samples_removed         = map(rse_outlier_removed, "removed")

commonBrNum_hippo_dentate1 = intersect(rse_outlier_removed_big$hippo$BrNum, rse_outlier_removed_big$dentate$BrNum)

rse_outlier_removed_big1 = rse_outlier_removed_big %>% purrr::imap(~{
  rse = SummarizedExperiment::subset(
    x      = .,
    select =(BrNum %in% commonBrNum_hippo_dentate1)
  )
  return(rse)
})


#Outlier gene removal
rse_preprocess_step3 = imap(rse_outlier_removed_big1, ~{
  expMat            = assays(.)$RPKM
  median.gene.exprs = matrixStats::rowMedians(expMat)
  names(median.gene.exprs) = rownames(expMat)
  
  #Gene filters
  gene_filter0 = as.vector(rowSums(expMat ==0) <= ncol(expMat)*0.2)                #Genes to keep (not more than 20% samples are zero)
  gene_filter1 = as.vector(median.gene.exprs   >= .1)                              #Genes to keep (median exp > 0.1)
  
  
  median.gene.exprs_filtered = median.gene.exprs[gene_filter0 & gene_filter1]
  gene_filter2 = as.vector(abs(scale(median.gene.exprs_filtered)) <  3)            #Genes to keep (z-scoreed_median_exp < 3)
  
  
  genes_to_keep  = names(median.gene.exprs_filtered)[gene_filter2]
  
  print(.y)
  print(rbind(table(gene_filter0), table(gene_filter1), table(gene_filter2)))
  print(paste0(length(genes_to_keep),"/",dim(expMat)[1]))
  
  rse_merged2 = SummarizedExperiment::subset(
    x      = .                                     ,
    subset = (gencodeID  %in% genes_to_keep)     &          
             (!grepl("^MT-",Symbol)        )                                      #Gene symbol starting with 'MT-' (for mitochondrial genes)
  )
  return(rse_merged2)
})


commonGenes_hippo_dentate = intersect(rownames(rse_preprocess_step3$hippo),rownames(rse_preprocess_step3$dentate))

rse_preprocess_step4 = rse_preprocess_step3 %>% purrr::imap(~{
  rse = SummarizedExperiment::subset(
    x      = .,
    subset = (gencodeID %in% commonGenes_hippo_dentate)
  )
  return(rse)
})

rse_preprocess_step5 = map(rse_preprocess_step4,~{
  assays(., withDimnames=FALSE)$exp_quant_normalised  = preprocessCore::normalize.quantiles(assays(.)$logRPKM) #Quantile normalise samples. Ref: http://jtleek.com/genstats/inst/doc/02_05_normalization.html
  return(.)
})

setwd("C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\New_Project\\newResults")
saveRDS(rse_preprocess_step5      , "rse_preprocessed_newproject_redone.rds")         #List of rse object preprocessed (sample and gene outlier removed. Samples and genes matched)
