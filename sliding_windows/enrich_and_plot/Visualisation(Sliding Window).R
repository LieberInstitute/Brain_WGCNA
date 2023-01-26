library(tidyverse)
#library(patchwork)
library(scales)
#(facetscales)  #For different scales for facets
library(qs)

width = list(
  PGC    = c("PGC"                             ),
  kb_000 = c("bp_0"     , "0"     , "kbp_0"    ),
  kb_020 = c("bp_20000" , "20000" , "kbp_20", "Ripke2014"),
  kb_040 = c("bp_40000" , "40000"              ),
  kb_050 = c("bp_50000" , "50000" , "kbp_50"   ),
  kb_060 = c("bp_60000" , "60000"              ),
  kb_080 = c("bp_80000" , "80000"              ),
  kb_100 = c("bp_1e+05" , "1e+05" , "kbp_100"  ),
  kb_150 = c("bp_150000", "150000", "kbp_150"  ),
  kb_200 = c("bp_2e+05" , "2e+05" , "kbp_200"  ),
  kb_250 = c("bp_250000", "250000", "kbp_250"  ),
  kb_300 = c("bp_3e+05" , "3e+05" , "kbp_300"  ),
  kb_350 = c("bp_350000", "350000", "kbp_350"  ),
  kb_400 = c("bp_4e+05" , "4e+05" , "kbp_400"  ),
  kb_450 = c("bp_450000", "450000", "kbp_450"  ),
  kb_500 = c("bp_5e+05" , "5e+05" , "kbp_500"  ),
  mb_001 = c("bp_1e+06" , "1e+06" , "Mbp_1"    ),
  mb_002 = c("bp_2e+06" , "1e+06" , "Mbp_2"    ),
  mb_003 = c("bp_3e+06" , "1e+06" , "Mbp_3"    ),
  mb_004 = c("bp_4e+06" , "1e+06" , "Mbp_4"    ),
  mb_005 = c("bp_5e+06" , "1e+06" , "Mbp_5"    ),
  mb_010 = c("bp_1e+07" , "1e+07" , "Mbp_10"   ),
  nn_200 = c("nn_200"   , "kbp_200_negative"   ),
  nn_500 = c("nn_500"   , "kbp_500_negative"   )
)

#Read NC/SCZ SummarizedExperiment objects for the expression matrices
##NC rse:
rse_NC = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Sliding_Window_Networks/results/rse_sva_cleaned_big_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")
##SCZ rse:
rse_SCZ = readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/results/rse_sva_SchizoNew_cleaned_big_noquantile_Neuclean_GEclean_newcuts_noPCA.rds")

age.df = rbind(
  rse_NC  %>% map_dfr(.id = "Tissue",~data.frame(Dx = "NC" , median.age = .x$Age)),
  rse_SCZ %>% map_dfr(.id = "Tissue",~data.frame(Dx = "SCZ", median.age = .x$Age))
) %>% 
  mutate(Tissue = case_when(
    Tissue == "dentate" ~ "DG",
    Tissue == "dlpfc"   ~ "DLPFC",
    Tissue == "hippo"   ~ "HP",
    Tissue == "caudate" ~ "CN"
  ))  %>% mutate(Tissue = factor(Tissue, levels = c("CN","DLPFC","HP","DG")))

#Read NC/SCZ fold change long form data 
ff_NC  <- readRDS("C:\\Users\\mpariha1\\Desktop\\grch38[PGC125new]\\sliding_windows\\enrichment_results (noquantile, Neuclean, sliding, NC)\\FoldEnrichment(NC).rds")
ff_SCZ <- readRDS("C:\\Users\\mpariha1\\Desktop\\grch38[PGC125new]\\sliding_windows\\enrichment_results (noquantile, Neuclean, sliding, SchizoNew)\\FoldEnrichment(SchizoNew).rds")

#Make a combined long form file for NC+SCZ
ff = bind_rows(list(NC = ff_NC, SCZ = ff_SCZ), .id = "Dx") %>% filter(list.type == "PGC3")
colnames(ff)[grep("bins",colnames(ff))] = c("bins.width")
ff$bins.width = as.character(stack(width)$ind[match(ff$bins.width,stack(width)$values)])

ff = ff %>% mutate(Tissue = case_when(
  tissue == "dlpfc"   ~ "DLPFC",
  tissue == "hippo"   ~ "HP",
  tissue == "caudate" ~ "CN",
  tissue == "dentate" ~ "DG",
))

#Caudate first window is removed from the plots [too far in age from other windows]
ff = ff[ff$newID != "NC.caudate__1" & ff$modules != "grey", ]##Caudate:218, dentate:46, dlpfc:221, hippo: 234

bins.width = "kb_200"

ff0 = ff %>% group_by(Dx, biotype, list.type, bins.width, newID) %>% 
  summarise(
    var.FoldChange      = var(Fold),
    sd.FoldChange       = sd(Fold),
    mean.FoldChange     = mean(Fold),
    median.FoldChange   = median(Fold),
    max.FoldChange      = max(Fold),
    min.FoldChange      = min(Fold),
    tissue_age          = sample(tissue_age,1),
    Tissue              = sample(Tissue,1),
    median.age          = sample(median.age,1),
    Fold                = list(Fold)
  ) %>% 
  ungroup() %>%
  arrange(list.type,biotype,bins.width, Tissue,median.age)

#Only retain SCZ AB enrichments
ff1 = ff0 %>% filter(list.type == "PGC3" & biotype == "all.biotypes" & bins.width %in% bin) %>% mutate(Tissue = factor(Tissue, levels = c("CN","DLPFC","HP","DG")))

p2 = ggplot(data = ff1, aes(x=median.age, color = interaction(Tissue,Dx), fill = interaction(Tissue,Dx))) +
  geom_smooth(data = ff1 %>% filter(Tissue == "DG"), method = "loess",aes(y = max.FoldChange, group = interaction(Dx,Tissue), linetype = Dx), alpha = 0.2, size = 1.1, span = 0.65, fullrange = F) +
  geom_smooth(data = ff1 %>% filter(Tissue != "DG"), method = "loess",aes(y = max.FoldChange, group = interaction(Dx,Tissue), linetype = Dx), alpha = 0.2, size = 1.1, span = 0.65, fullrange = F) + 
  facet_grid(Tissue~.) + 
  #geom_rug(data = ff1    %>% subset(Dx == "NC") , sides = "t",aes(color = interaction(Tissue,Dx)), alpha = 0.55, size = 0.75, outside = T, show.legend = F)+  #NC  windows
  #geom_rug(data = ff1    %>% subset(Dx == "SCZ"), sides = "t",aes(color = interaction(Tissue,Dx)), alpha = 0.55, size = 0.75, outside = F, show.legend = F)+  #SCZ windows
  geom_rug(data = age.df %>% subset(Dx == "NC") , sides = "t",aes(color = interaction(Tissue,Dx)), alpha = 0.5, length = unit(0.6,"lines"), size = 0.6, outside = F, show.legend = F)+  #NC  subjects
  geom_rug(data = age.df %>% subset(Dx == "SCZ"), sides = "b",aes(color = interaction(Tissue,Dx)), alpha = 0.5, length = unit(0.6,"lines"), size = 0.6, outside = F, show.legend = F)+  #SCZ subjects
  coord_cartesian(clip = "off")+
  #scale_x_continuous(expand = expansion(mult = c(0.005,NA)))+
  scale_y_continuous(expand = expansion(mult = c(0.08)))+
  theme_bw(base_size = 14)  + 
  theme(
    axis.text = element_text(size = 10.5),
    strip.text = element_text(size = 11, margin = margin(2,2,2,2)),
    axis.title = element_text(size = 11),
    plot.margin=grid::unit(c(3,3,1,3), "pt"),
    plot.title = element_text(hjust = 0.5),
    panel.spacing.y = unit(0.8, "lines"),
    panel.border    = element_rect(size = 0.35, color = prismatic::clr_lighten("black")),
    legend.position = "bottom",
    axis.ticks.length.x = unit(0.75,"lines"),
    axis.ticks.length.y = unit(0.2,"lines"),
    strip.background = element_rect(size = 0.25, fill = "grey90", color = "black"),
    legend.box.spacing = margin(0,0,0,0,"pt"),
    panel.grid.major = element_line(size = 0.4),
    panel.grid.minor = element_line(size = 0.4))+
  ylab(label = "Max  (Fold change)") +xlab("Median age")+
  ggtitle(paste0("SCZ PGC3 FoldChange @",bin," (no Grey)"))+
  guides(color = "none", fill = "none")+
  labs(linetype = "dx")
print(p2)
#ggsave(p2+ggtitle(label = NULL)            , filename = paste0("sliding_window_max.fold.change_", bin,".svg"),device = svg, width = 8, height = 6)

