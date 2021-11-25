####################################################################################################
#      Script to Combine WGCNA output from our networks and other previous published networks      #
#           "Regional-coexpression", "Age-period" and "Cell-population enrichment" studies         #
#                                                                                                  #
####################################################################################################



new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\old-to-new Network names all.csv") %>% tibble::deframe()

all_networks = list(
  WGCNAall      = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\WGCNAnetworks.rds"),
  
  Fromer2016    = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Fromer 2016/Fromer2016(Madhur).rds"),
  Gandal2018    = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Gandal/Gandal2.science.aad6469(Gandal2018).rds"),
  Gandal2018PE  = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Gandal/Gandal1.science.aat8127 (Gandal2018PE).rds"),
  Li2018        = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Li2018/Li2018.rds"),
  Pergola2017   = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/DRD2 (Pergola2017)/DRD2final(Madhur).rds"),
  Pergola2019   = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Pergola2019/Pergola2019final(Madhur).rds"),
  Pergola2020   = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/MIR137 (Pergola2020)/mir137final(Madhur).rds"),
  Radulescu2020 = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Eugenia (Radulescu2020)/Eugenia(Madhur).rds"),
  Walker2019    = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Walker2019/Walker2019.rds"),
  Werling2020   = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Other networks/Other Network Data/Werling2020/Werling2020.rds")
)


final = map(all_networks,~.x %>% tibble::rownames_to_column("EnsemblID")) %>% purrr::reduce(merge,by = "EnsemblID",all = T) %>% tibble::column_to_rownames("EnsemblID")

