library(dplyr)
library(purrr) ; library(gapminder)
library(SummarizedExperiment); library(qs);library(ggplot2); library(gridExtra);library(pryr); library(tibble); library(furrr); library(tidyr); library(readr)
library(dplyr); library(data.table); library(grid); library(limma); library(pheatmap); library(RColorBrewer); library(broom); library(MASS)
rm(list = ls())
sampleInfo_sliding_window <- readRDS("../Across age/sampleInfo_sliding_window.rds")
sw_ages = map_dbl(sampleInfo_sliding_window[7 : length(sampleInfo_sliding_window)], ~ {
  # print(which(..1 == 0))
  median(sampleInfo_sliding_window$Age[which(..1 == 1)])
})

st = "Z_stat"
sw_rlm_stats <- readRDS("results/sw_modenrMAGMA_rlm_stats.rds")
sw_stat_list = map(sw_rlm_stats, ~ ..1[[st]])
sw_stat_list[["median_age"]] = sw_ages
sw_stat_list[["network"]] = names(sw_ages)

sw_stat_wide = as.data.frame(sw_stat_list) 
sw_stat_wide = sw_stat_wide %>% as.data.frame()
sw_stat_wide["dx.tissue"] = gsub("__.*","",sw_stat_wide$network)
sw_stat_wide["dx"] = gsub("\\..*","",sw_stat_wide$network)
sw_stat_wide["dx"] = gsub("SchizoNew","SCZ",sw_stat_wide$dx)
sw_stat_wide["tissue"] = gsub(".*\\.","",sw_stat_wide$dx.tissue)
sw_stat_wide["tissue"] = gsub("caudate","CN",sw_stat_wide$tissue)
sw_stat_wide["tissue"] = gsub("dlpfc","DLPFC",sw_stat_wide$tissue)
sw_stat_wide["tissue"] = gsub("hippo","HP",sw_stat_wide$tissue)
sw_stat_wide["tissue"] = gsub("dentate","DG",sw_stat_wide$tissue)
sw_stat_wide = sw_stat_wide %>%
  # filter(!grepl('',network )) %>%
  filter(grepl('',network )) #%>%
  # select_if(~ any(. > 6))
# sw_modenrHMAGMA_corr_wide_sig = sw_modenrHMAGMA_corr_wide %>% select_if(~ any(. > 6))

sw_stat_long = data.table::melt(setDT(sw_stat_wide),
                                                  id.vars = c("network","median_age", "dx.tissue", "dx", "tissue"),
                                                  variable.name = "DRONC") %>%
  filter(grepl('',DRONC ) & grepl('',network ))
  # filter(grepl('neuron|nervous|synapse|dendritic',GO_term ) & grepl('',network ))
  # filter(grepl('immune',GO_term ) & grepl('',network ))
  # filter(grepl('chromatin|DNA',GO_term ) & grepl('',network ))
# filter(grepl('projection',GO_term ) & grepl('',network ))
# filter(grepl('_ion_|transport',GO_term ) & grepl('',network ))
# filter(grepl('development|differentiation|generation',GO_term ) & grepl('',network ))

# filter(grepl('neu|synapse',GO_term ) & grepl('',network ))

cc =c(c(
  "#1F78B4" , "cyan","#A6CEE3"  , "darkblue", 
  "#FB9A99" ,"#E31A1C", "#FDBF6F"  ,"darkgreen",
  "#6A3D9A" ,"#CAB2D6" ,
  "palegreen" , "pink" , "grey60", "black", "brown"
  , "darkmagenta", "gold"
)
)
sw_stat_long$group = rep("Cell type", times = nrow(sw_stat_long))
sw_stat_long$value = as.numeric(sw_stat_long$value)
sw_stat_long$tissue <- factor(sw_stat_long$tissue , levels = c("CN", "DLPFC", "HP", "DG"))
sw_stat_long$DRONC <- factor(sw_stat_long$DRONC , 
                             levels = c("exCA",  "exDG" , "exPFC" ,"GABA"  ,
                                        "ASC" ,  "ODC" , "MG" ,"END"  ,   
                                        "NSC"   ,  "OPC" ))
# # jpeg(file="saving_plot1.jpeg")
# gg1 = ggplot(data = sw_stat_long , aes(x=median_age, y= value,  group=interaction(dx.tissue, DRONC)))  +
#   scale_color_manual(values=cc) +
#   scale_fill_manual(values=cc) +
#   geom_smooth(aes(color=DRONC, fill =DRONC, linetype = dx),alpha = .05,method=loess,size=1.2, span=.0000001, level = .99)+
#   facet_grid(rows = vars(tissue))+ 
#   geom_hline(yintercept=c(3,-3), linetype="dashed", color = "red")+  
#   ylab("Z.stat for predicting MAGMA enrichment")+
#   xlab("Median age")+ 
#   # theme_bw()+
#   # geom_label()+
#   # geom_point(aes(fill  = GO_term, color = GO_term,shape = dx), alpha = 0.4) +
#   theme(
#     panel.background = element_rect(fill = "grey99",
#                                     # colour = "lightblue",
#                                     size = 0.5, linetype = "solid"
#                                     ),
#     legend.position="bottom",legend.text=element_text(size=9))
# # dev.off()
# x11(); gg1


jpeg("plots/DRONC_sw_2_samescale.jpeg",width =900, height = 1400)
gg1 = ggplot(data = sw_stat_long , aes(x=median_age, y= value,  group=interaction(dx.tissue, DRONC)))  +
  scale_color_manual(values=cc) +
  scale_fill_manual(values=cc) +
  geom_smooth(aes(color=DRONC, fill =DRONC, linetype = dx),alpha = .05,method=loess,size=2.1, span=.0000001, level = .99)+
  facet_grid(rows = vars(tissue)
             , cols = vars(group)
             )+
  geom_hline(yintercept=c(3,-3), linetype="dashed", color = "red", size =1.3)+ 
  coord_cartesian(ylim=c(-7, 23))+
  ylab("Z.stat for predicting MAGMA enrichment")+
  xlab("Median age")+ 
  theme_bw()+
  labs(col = "Cell type", fill = "Cell type")+
  # guides(fill=guide_legend(title="Cell type")) +
  # geom_label()+
  # geom_point(aes(fill  = DRONC, color = DRONC,shape = dx), alpha = 0.4) +
  theme(
    text = element_text(size = 50),
    # panel.background = element_rect(fill = "grey99",
    #                                 # colour = "lightblue",
    #                                 size = 0.5, linetype = "solid"
    # ),
    legend.position="right",legend.text=element_text(size=50), legend.key.size = unit(2, "cm"))
print(gg1)
dev.off()

