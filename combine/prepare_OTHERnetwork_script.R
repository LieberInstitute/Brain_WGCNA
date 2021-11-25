####################################################################################################
#               Script to  collate and combine  network data from other published networks         #
#                                "Regional-coexpression" study                                     #
#                                                                                                  #
####################################################################################################


library(purrr)
library(SummarizedExperiment)
library(WGCNA)
library(limma)
library(tidyr)
library(gtools)
library(purrr)
library(dplyr)
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

WGCNAfn = function(.) {
  WGCNAcolors = c("grey",WGCNA::standardColors(200)) %>% set_names(0:200)
  return(WGCNAcolors[match(.,names(WGCNAcolors))])
  }


####################
####  prepare 'Radulescu2020' network

eu = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Eugenia (Radulescu2020)\\EugeniaSupplementary_table2_modules_revised_07022018.csv")
eu = eu %>% dplyr::select(-c(Order,Gene_Symbol))

eu = eu %>% tibble::column_to_rownames("ENSEMBL")
colnames(eu)[1]  = "modules"
colnames(eu)[-1] = substring(colnames(eu)[-1],4)
colnames(eu) = paste0("Radulescu2020|",colnames(eu))

saveRDS(eu, "C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Eugenia\\Eugenia(Madhur).rds")


####################
####  prepare 'Fromer2016 case' and 'Fromer2016 control' networks

case    <- read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Fromer 2016\\Case.csv")
control <- read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Fromer 2016\\Control.csv")

final_case    = case    %>% dplyr::select(-c(Gene.Symbol, Module)) %>% set_colnames(paste0("Fromer2016_case|"   ,c("Ensembl","modules","kTotal","kWithin","kOut"))) %>% tibble::column_to_rownames("Fromer2016_case|Ensembl")
final_control = control %>% dplyr::select(-c(Gene.Symbol, Module)) %>% set_colnames(paste0("Fromer2016_control|",c("Ensembl","modules","kTotal","kWithin","kOut"))) %>% tibble::column_to_rownames("Fromer2016_control|Ensembl")
final_Fromer  = merge(final_case, final_control, by = 0, all.x=T)  %>% tibble::column_to_rownames("Row.names")
saveRDS(final_Fromer, "C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Fromer 2016\\Fromer2016(Madhur).rds")

####################
####  Gandal1.science.aat8127
####  Transcriptome-wide isoform-level dysregulation in ASD, schizophrenia, and bipolar disorder
####  PsychENCODE
####  prepare 'Gandal2018PE' and 'Gandal2018PE_cs' networks

dat1 = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Gandal\\Gandal1.science.aat8127.csv")
dat2 = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Gandal\\Gandal1.science.aat8127_cs.csv")

dat1 = dat1 %>% tibble::column_to_rownames("ensembl_gene_id")
colnames(dat1)[1]  = "modules"
colnames(dat1)[-1] = unname(WGCNAfn(substring(colnames(dat1)[-1],4)))
colnames(dat1) = paste0("Gandal2018PE|",colnames(dat1))

dat2 = dat2 %>% tibble::column_to_rownames("ensembl_gene_id")
colnames(dat2)[1]  = "modules"
colnames(dat2)[-1] = unname(WGCNAfn(substring(colnames(dat2)[-1],4)))
colnames(dat2) = paste0("Gandal2018PE_cs|",colnames(dat2))

all.equal(rownames(dat1), rownames(dat2))
dat = cbind(dat1,dat2)

saveRDS(dat, "C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Gandal\\Gandal1.science.aat8127.rds")


####################
#### Gandal2.science.aad6469
#### Shared molecular neuropathology across major psychiatric disorders parallels polygenic overlap
####  prepare 'Gandal2018' network
dat = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Gandal\\Gandal2.science.aad6469.csv")

dat = dat %>% tibble::column_to_rownames("ensembl_gene_id")
colnames(dat)[1]  = "modules"
colnames(dat)[-1] = strsplit2(colnames(dat)[-1],"\\.")[,3]
colnames(dat) = paste0("Gandal2018|",colnames(dat))

saveRDS(dat, "C:\\Users\\Madhur\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Gandal\\Gandal2.science.aad6469.rds")


####################
#### Walker2019
#### Genetic Control of Expression and Splicing in Developing Human Brain Informs Disease Mechanisms
####  prepare 'Walker2019' network
dat = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Walker2019\\table3(gene-module assignment and module membership).csv")

dat1 = unique(dat) %>% set_rownames(.$"ENSG")#Removing 3 repeated rows
dat1[,c("ENSG","Gene","eQTL_nPval")] = NULL  
colnames(dat1)[1]  = "modules"
#colnames(dat1)[-1] = strsplit2(colnames(dat1)[-1],"\\.")[,3]
colnames(dat1) = paste0("Walker2019|",colnames(dat1))

saveRDS(dat1, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Walker2019\\Walker2019.rds")


####################
#### Werling2020:10.1016/j.celrep.2020.03.053
#### Whole-Genome and RNA Sequencing Reveal Variation and Transcriptomic Coordination in the Developing Human Prefrontal Cortex
#### prepare 'Werling2020' network
##   NIHMS1582939-supplement-4
##   "grey60" modules is named as silver in the files
dat = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Werling2020\\data.csv")
dat = dat %>% tibble::column_to_rownames("gene_id")
dat$gene_name = NULL
colnames(dat) = gsub("kME","",colnames(dat))
colnames(dat)[1] = "modules"
colnames(dat) = paste0("Werling2020|",colnames(dat))
saveRDS(dat1, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Werling2020\\Werling2020.rds")



####################
####  prepare 'Li2018' network
##    Modified from Table S10
##    No grey gene list
dat = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Li2018\\gene-module list (wide).csv", header = F)
dat = dat %>% tibble::column_to_rownames("V1")


dat1 = apply(dat,c(1,2), function(x) {
  #browser()
  if(x != "") {return(strsplit2(x,"\\|")[,1])} else {return(NA)}
  
})

dat2 = t(dat1)

dat3 = map(as.data.frame(dat2), ~ list(.x[!is.na(.x)])) %>% flatten %>% .[sort(names(.))]
names(dat3) = WGCNA::standardColors(length(dat3))

dat4 = tibble::enframe(dat3) %>% unnest_longer(value) %>% set_colnames(c("Li2018|modules","ensembl")) %>% tibble::column_to_rownames("ensembl")

saveRDS(dat4, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Walker2019\\Walker2019.rds")



####################
#### 10.1038/tp.2016.253
#### DRD2 co-expression network and a related polygenic index predict imaging, behavioral and clinical phenotypes linked to schizophrenia
#### prepare 'Pergola2017' network
##   Beta = 4
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\DRD2 (Pergola2017)\\data\\1. WGCNA_AIC199bis_Blom_b4_size40.RData")

acc = accession %>% tibble::rownames_to_column("probeID") %>% mutate(Accession = gsub(",","",Accession), Type = substr(probeID,1,3)) 

set.seed(07272020)
#Converting accession numbers to EnsemblID
bb1 = bitr(acc$Accession, fromType = "ACCNUM", toType = c("ENSEMBL"), OrgDb = "org.Hs.eg.db") %>%
  group_by(ACCNUM) %>% slice_sample()
#Converting gene symbol to EnsemblID
bb2 = bitr(acc$Gene, fromType = "SYMBOL", toType = c("ENSEMBL"), OrgDb = "org.Hs.eg.db") %>%
  group_by(SYMBOL) %>% slice_sample()

acc = acc %>% mutate( ENSEMBL1 = bb1$ENSEMBL[match(Accession, bb1$ACCNUM)],
                      ENSEMBL2 = bb2$ENSEMBL[match(Gene     , bb2$SYMBOL)],
                      ENSEMBL  = ifelse(is.na(ENSEMBL1), ENSEMBL2, ENSEMBL1)) %>%
  #dplyr::select(c(probeID,Type, ENSEMBL))                                     %>%
  filter(!is.na(ENSEMBL))

#Selecting only single probeID per EnsemblID in preferential order: hHC, hHA, hHR
acc = acc %>% group_by(ENSEMBL) %>% mutate(
  test = case_when(
    ("hHC" %in% Type)                                            ~ 1,
    !("hHC" %in% Type) & ("hHA" %in% Type)                       ~ 2,
    !("hHC" %in% Type) & !("hHA" %in% Type) & ("hHR" %in% Type)  ~ 3,
    TRUE                                                         ~ 0
  )) %>% ungroup() %>% 
  mutate(final = ifelse(((test == 1) & (Type == "hHC")) | ((test == 2) & (Type == "hHA")) | ((test == 3) & (Type == "hHR")), T, F)) %>%
  filter(final) %>% group_by(ENSEMBL) %>% dplyr::slice(1) %>% 
  dplyr::select(c(probeID, Gene, Accession, Type, ENSEMBL))


ss = IMk_module[,c("moduleColors", "kTotal", "kWithin", "kOut", "kDiff")]

eigen = eigen[2:nrow(eigen),]
all.equal(rownames(eigen), rownames(aic199bis.blom))
tt = signedKME(aic199bis.blom, datME = eigen)
colnames(tt) = substring(colnames(tt),4)

final = merge(ss,tt, by = 0, all = T)
colnames(final)[1:2] = c("gene","modules") 

final1 = final %>% mutate(ENSEMBL = acc$ENSEMBL[match(gene, acc$probeID)]) %>% filter(!is.na(ENSEMBL))
final2 = final1 %>% tibble::column_to_rownames("ENSEMBL") %>% dplyr::select(-gene) %>% set_colnames(paste0("Pergola2017|",colnames(.)))

saveRDS(final2, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Pergola2017.rds")

####################
#### Pergola2019 : https://doi.org/10.1016/j.biopsych.2019.03.981
#### DRD2 co-expression network and a related polygenic index predict imaging, behavioral and clinical phenotypes linked to schizophrenia
#### prepare 'Pergola2019' network
##   Unsigned, Pearson, beta: 5
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Pergola2019\\data\\2. WGCNA_k5_nu0_HK13_LIBD.RData")
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Pergola2019\\data\\2. LIBD.moduleMembership.RData")
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Pergola2019\\data\\2. LIBD.RUV.n343.g20993..k5.nu0.HK13.RData")

ss = IMk_module_LIBD[,c("moduleColors", "kTotal", "kWithin", "kOut", "kDiff")]
tt = MM.LIBD
colnames(tt) = substring(colnames(tt),3)

final = merge(ss,tt, by = 0, all = T)
colnames(final)[1:2] = c("gene","modules") 

Pergola2019 = final %>% tibble::column_to_rownames("gene") %>% set_colnames(paste0("Pergola2019|",colnames(.)))
saveRDS(Pergola2019, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\Pergola2019final(Madhur).rds")


####################
####  Pergola2020: 10.1101/2020.08.03.230227
####  mir137 :A miR-137-related biological pathway of risk for Schizophrenia is associated with human brain emotion processing
##    Unsigned, Spearman, beta: 7
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\MIR137 (Pergola2020)\\data\\3.  wgcna_cauc_net.RData")
ss = IMk_module[,c("moduleColors", "kTotal", "kWithin", "kOut", "kDiff")]

eigen = eigen[3:nrow(eigen),]
all.equal(rownames(eigen), rownames(wgcna_cauc))
tt = signedKME(wgcna_cauc, datME = eigen)
colnames(tt) = substring(colnames(tt),4)
final = merge(ss,tt, by = 0, all = T)
colnames(final)[1:2] = c("gene","modules") 

set.seed(07272020)
bb = bitr(final$gene, fromType = "SYMBOL", toType = c("ENSEMBL"), OrgDb = "org.Hs.eg.db")
bb = bb %>% group_by(SYMBOL) %>% slice_sample() %>% ungroup %>% group_by(ENSEMBL) %>% slice(1)
#bb = bb %>% group_by(SYMBOL) %>% summarise(ENSEMBL_all = list(ENSEMBL), n = lengths(ENSEMBL_all), ENSEMBL = unlist(ENSEMBL_all)[1],.groups = "keep") %>% slice(1) %>% arrange(-n) %>% ungroup %>% .[!duplicated(.$ENSEMBL),]

mir137 = final[final$gene %in% bb$SYMBOL,]
rownames(mir137) = bb$ENSEMBL[match(mir137$gene, bb$SYMBOL)]
mir137$gene = NULL
colnames(mir137) = paste0("Pergola2020|",colnames(.))
saveRDS(mir137, "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Other networks\\Other Network Data\\mir137final(Madhur).rds")
