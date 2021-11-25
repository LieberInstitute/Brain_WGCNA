library(sva)
library(jaffelab)
set.seed(123)
library(SummarizedExperiment)
library(purrr)

setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Stem_Cell_Project_2021/")

rse_preprocessed = readRDS("redone/results/rse_preprocessed(redone).rds")        #Sample and gene outlier removed + samples-Quant-normalized
rse_preprocessed = map(rse_preprocessed,~ {assays(.x)$quant.norm.logRPKM = NULL; return(.x)})

###############
get.cleaned.obj = function(rse, rse.name, protect,  assay.to.clean){
  exp = assays(rse)[[assay.to.clean]]
 
  pcobj  = prcomp(t(exp),center = T, scale. = T)
  all.equal(rownames(as.data.frame(colData(rse))), rownames(pcobj$x))
  
  dat = cbind.data.frame(as.data.frame(colData(rse)),pcobj$x)
  mod = model.matrix(~  as.factor(RealDx) + h_mitoRate + rRNA_rate + totalAssignedGene + PC1
                     ,
                     data = dat)
  
  #Check of multi-colinearity between known covariates and PCs
  print(sum(alias(lm(t(exp)~mod))$Complete))
  
  
  cleaned_matrix_test = t(lm(t(exp) ~ mod)$residuals)
  cleaned_matrix = cleaningY(exp, mod = mod, P = protect)
  
  assays(rse)$svacleaned = cleaned_matrix
  
  check.shapiro1 = apply(cleaned_matrix,1, function(r) shapiro.test(r)$p.value)
  print(paste0("Non-normal genes (p <0.05) before ranknorm: ", sum(check.shapiro1 <0.05), " / ",length(check.shapiro1)))
  
  exp_ranknorm           = t(apply(cleaned_matrix,1, RNOmni::RankNorm))  #rankNorm each gene in the matrix
  assays(rse)$ranknorm   = exp_ranknorm
  
  check.shapiro2 = apply(exp_ranknorm,1, function(r) shapiro.test(r)$p.value)
  print(paste0("Non-normal genes (p <0.05) after ranknorm: ", sum(check.shapiro2 <0.05), " / ",length(check.shapiro2)))
  
  
  return(rse)
}

rse.sva.cleaned.big_logRPKM = list()
for (l in names(rse_preprocessed)){
  rse.sva.cleaned.big_logRPKM[[l]]    = get.cleaned.obj(rse_preprocessed[[l]], l, protect = 2,  assay.to.clean = "logRPKM")
}


setwd("redone/results/")
saveRDS(rse.sva.cleaned.big_logRPKM,"rse_sva_cleaned_protectDx(redone).rds")
