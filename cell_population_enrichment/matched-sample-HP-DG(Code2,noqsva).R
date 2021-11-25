library(sva)
library(jaffelab)
set.seed(123)
library(SummarizedExperiment)
library(purrr)
library(limma)

setwd("C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\New_Project\\newResults")
rse_preprocessed = readRDS("rse_preprocessed_newproject_redone.rds")        #Sample and gene outlier removed. Sample and genes matched.

get.cleaned.obj = function(rse, rse.name, protect,  assay.to.clean){
  rse$RIN = sapply(rse$RIN, mean)
  
  if(!grepl("dentate",rse.name)){
    mod = model.matrix(~ Age + as.factor(Sex) + mitoRate + rRNA_rate + totalAssignedGene + RIN +
                         C1 + C2 + C3 + C4+ C5+  C6 + C7 + C8 + C9 + C10,
                       data = as.data.frame(colData(rse)))
  } else {
    #Skipping RIN correction for dentate due to NAs
    mod = model.matrix(~ Age + as.factor(Sex) + mitoRate + rRNA_rate + totalAssignedGene +
                         C1 + C2 + C3 + C4+ C5+  C6 + C7 + C8 + C9 + C10,
                       data = as.data.frame(colData(rse)))
  }
  
  exp = assays(rse)[[assay.to.clean]]
  
  #Check of multi-colinearity between known covariates and PCs
  print(sum(alias(lm(t(exp)~mod))$Complete))
  
  
  cleaned_matrix = cleaningY(exp, mod = mod, P = protect)
  assays(rse)$svacleaned = cleaned_matrix
    
  check.shapiro = apply(cleaned_matrix,1, function(r) shapiro.test(r)$p.value)
  print(paste0("Non-normal genes (p <0.05) before ranknorm: ", sum(check.shapiro <0.05), " / ",length(check.shapiro)))
  
  exp_ranknorm           = t(apply(cleaned_matrix, 1, RNOmni::RankNorm))  #rankNorm each gene in the matrix
  assays(rse)$ranknorm   = exp_ranknorm
  
  return(rse)
}

##Regress confounders and apply genewise rank-normalisation
rse.sva.cleaned.big_quant.norm = list()
for (l in names(rse_preprocessed)){
  rse.sva.cleaned.big_quant.norm[[l]] = get.cleaned.obj(rse_preprocessed[[l]], l, protect = 1,  assay.to.clean = "exp_quant_normalised")
}

setwd("C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\New_Project\\newResults")
saveRDS(rse.sva.cleaned.big_quant.norm ,"rse_sva_cleaned_big_newproject_noAgeSqr_noPCAremoved_redone.rds")

