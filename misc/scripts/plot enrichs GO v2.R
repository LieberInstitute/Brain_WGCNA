library(dplyr)
library(purrr)
library(SummarizedExperiment); library(qs);library(ggplot2); library(gridExtra);library(pryr); library(tibble); library(furrr); library(tidyr); library(readr)
library(dplyr); library(data.table); library(grid); library(limma); library(pheatmap); library(RColorBrewer); library(broom); library(MASS)
library(RColorBrewer)

rm(list = ls())
sampleInfo_sliding_window <- readRDS("../Across age/sampleInfo_sliding_window.rds")
sw_ages = map_dbl(sampleInfo_sliding_window[7 : length(sampleInfo_sliding_window)], ~ {
  # print(which(..1 == 0))
  median(sampleInfo_sliding_window$Age[which(..1 == 1)])
})

st = "Z_stat"
sw_rlm_stats = readRDS("results/sw_modenrMAGMA_GOterms_CNS_rlm_stats.rds")
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
                                                  variable.name = "GO_term") #%>%
  # filter(grepl('',GO_term ) & grepl('',network ))
  # filter(grepl('synapse|dendritic',GO_term ) & grepl('',network ))
  # filter(grepl('immune',GO_term ) & grepl('',network ))
  # filter(grepl('projection',GO_term ) & grepl('',network ))
# filter(grepl('neuron|nervous',GO_term ) & grepl('',network ))%>%
# filter(grepl('_ion_|transport',GO_term ) & grepl('',network ))
# filter(grepl('development|differentiation|generation',GO_term ) & grepl('',network ))

sw_stat_long$GO_term<- gsub("^.{0,11}", "", sw_stat_long$GO_term)
sw_stat_long$GO_term_group = rep(NA, times = nrow(sw_stat_long))
sw_stat_long$GO_term_group[grepl('synapse|dendritic',sw_stat_long$GO_term )] = "Synapse and axon"
sw_stat_long$GO_term_group[grepl('projection',sw_stat_long$GO_term ) & 
                             grepl('neuron|nervous',sw_stat_long$GO_term )] = "Synapse and axon"
sw_stat_long$GO_term_group[grepl('development|differentiation|generation|process',sw_stat_long$GO_term )& 
                             grepl('neuron|nervous',sw_stat_long$GO_term )] = "NS development"
# sw_stat_long$GO_term_group[grepl('development|differentiation|generation',sw_stat_long$GO_term )& 
#                              grepl('neuron|nervous',sw_stat_long$GO_term )] = "Neuronal development"
sw_stat_long$GO_term_group[grepl('_ion_',sw_stat_long$GO_term )] = "Ions"
# sw_stat_long$GO_term_group[!grepl('Neuronal development|Synapse and axon',sw_stat_long$GO_term_group )&
#                              grepl('neuron|nervous',sw_stat_long$GO_term )] = "Other neuronal"
# sw_stat_long$GO_term_group[grepl('_ion_|transport',sw_stat_long$GO_term )& 
                             # grepl('',sw_stat_long$GO_term )] = "Ion transport"

sw_stat_long = sw_stat_long[!is.na(sw_stat_long$GO_term_group),]
sw_stat_long$GO_term <- factor(sw_stat_long$GO_term , levels = c(
  unique(sw_stat_long$GO_term[sw_stat_long$GO_term_group == "Synapse and axon"]), 
  unique(sw_stat_long$GO_term[sw_stat_long$GO_term_group == "NS development"]), 
  unique(sw_stat_long$GO_term[sw_stat_long$GO_term_group == "Ions"])
  ))
sw_stat_long$value = as.numeric(sw_stat_long$value)
sw_stat_long$tissue <- factor(sw_stat_long$tissue , levels = c("CN", "DLPFC", "HP", "DG"))
sw_stat_long$GO_term_group <- factor(sw_stat_long$GO_term_group , levels = c("Synapse and axon", 
                                                                             "NS development", 
                                                                             "Ions"))
cc = c(brewer.pal(n = 9, 
                   name = "Greens")[c(3,5,7)],
        brewer.pal(n = length(unique(sw_stat_long$GO_term[sw_stat_long$GO_term_group == "NS development"]))+1, 
                   name = "Purples")[2:9],
        brewer.pal(n = length(unique(sw_stat_long$GO_term[sw_stat_long$GO_term_group == "Ions"]))+1, 
                   name = "OrRd")[2:9]
        )


jpeg("plots/GO_sw_celltypescale.jpeg",width =2600, height = 1400)
gg1 = ggplot(data = sw_stat_long , aes(x=median_age, y= value,  group=interaction(dx.tissue, GO_term)))  +
  scale_color_manual(values=cc) +
  scale_fill_manual(values=cc) +
  geom_smooth(aes(color=GO_term, fill =GO_term, linetype = dx),alpha = .05,method=loess,size=2.1, span=.0000001, level = .99)+
  facet_grid(rows = vars(tissue), cols = vars(GO_term_group))+ 
  geom_hline(yintercept=c(-3,3), linetype="dashed", color = "red", size =1.3)+
  ylab("Z.stat for predicting MAGMA enrichment")+
  xlab("Median age")+ 
  coord_cartesian(ylim=c(-7, 23))+ 
  theme_bw()+
  # geom_label()+
  # geom_point(aes(fill  = GO_term, color = GO_term,shape = dx), alpha = 0.4) +
  theme(
    text = element_text(size = 50),
    # panel.background = element_rect(fill = "grey99",
    #                                 # colour = "lightblue",
    #                                 size = 0.5, linetype = "solid"
                                    # ),
    legend.position="right",legend.text=element_text(size=50), legend.key.size = unit(2, "cm"))
print(gg1)
dev.off()
x11(); gg1

