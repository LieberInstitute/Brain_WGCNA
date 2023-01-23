###############################################################################################
#                        Script to clean confounder effects                                   #
#                      "Regional-coexpression" and "Age-period" studies                       #
#                                   (NC samples only)                                         #
#                                  (Revision: Includes prenatal networks)                     #
###############################################################################################




library(sva)
library(jaffelab)
set.seed(123)
library(SummarizedExperiment)
library(purrr)

setwd("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/")


rse_preprocessed = readRDS("results/rse_preprocessed_GEclean_newcuts(dlpfc new groups).rds")        #Sample and gene outlier removed
samples_removed  = readRDS("results/outlier_samples_removed(dlpfc new groups).rds")

###############
get.cleaned.obj = function(rse, rse.name, protect,  assay.to.clean){
  rse$RIN = sapply(rse$RIN, mean)
  
  if(!grepl("dentate",rse.name)){
    mod = model.matrix(~ Age + I(scale(scale(Age)^2)) + as.factor(Sex) + mitoRate + rRNA_rate + totalAssignedGene + RIN +
                         C1 + C2 + C3 + C4+ C5+  C6 + C7 + C8 + C9 + C10 + neu,
                       data = as.data.frame(colData(rse)))
  } else {
    #Skipping RIN correction for dentate due to NAs
    mod = model.matrix(~ Age + I(scale(scale(Age)^2)) + as.factor(Sex) + mitoRate + rRNA_rate + totalAssignedGene +
                         C1 + C2 + C3 + C4+ C5+  C6 + C7 + C8 + C9 + C10 + neu,
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
rse.sva.cleaned.big_logRPKM = list()
for (l in names(rse_preprocessed)){
  rse.sva.cleaned.big_logRPKM[[l]]    = get.cleaned.obj(rse_preprocessed[[l]], l, protect = 3,  assay.to.clean = "logRPKM")
}

#Store parsed networks in a separate list()
rse.sva.cleaned.age.grps = samples_removed %>% purrr::imap(~{
  tissue = case_when(
    grepl("dlpfc.O",.y) ~ "dlpfc.O",
 #   grepl("dlpfc.A",.y) ~ "dlpfc.A",
    grepl("hippo"  ,.y) ~ "hippo",
    grepl("caudate",.y) ~ "caudate",
    grepl("dentate",.y) ~ "dentate"
  )
  
  if (tissue != "dentate"){
    age_min = as.numeric(strsplit2(strsplit2(.y, "__")[, 2],"_")[1])
    age_max = as.numeric(strsplit2(strsplit2(.y, "__")[, 2],"_")[2])
  } else {
    age_min = -1
    age_max = 150
  }
  
  temp_full = rse.sva.cleaned.big_logRPKM[[tissue]]
  
  temp_full_agerange = temp_full[,temp_full$Age>age_min & temp_full$Age<= age_max]
  
  temp      = temp_full_agerange[!colnames(temp_full_agerange) %in% .x]
  
  check.shapiro = apply(assays(temp)$ranknorm,1, function(r) shapiro.test(r)$p.value)
  print(paste0(.y,". # ",dim(temp)[2]))
  print(paste0("Non-normal genes (p <0.05) in ranknorm -> age-parsed sets: ", sum(check.shapiro <0.05), " / ",length(check.shapiro)))
  #if (dim(temp)[2] > 0) {gr[[i]] = temp}
  return(temp)
  
})



setwd("results/")
saveRDS(rse.sva.cleaned.age.grps,"rse_sva_cleaned_age_grps_noquantile_Neuclean_GEclean_newcuts_noPCA(dlpfc new groups).rds")
saveRDS(rse.sva.cleaned.big_logRPKM ,"rse_sva_cleaned_big_noquantile_Neuclean_GEclean_newcuts_noPCA(dlpfc new groups).rds")