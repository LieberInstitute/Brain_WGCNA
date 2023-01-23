####################################################################################################
#             Script to plot main enrichment plots for SCZ risk modules for all networks           #
#        Script to plot GO, TF and Pathology enrichments for SCZ risk modules for all networks     #
#           "Regional-coexpression", "Age-period" and "Cell-population enrichment" studies         #
#                                                                                                  #
#                              Script to plot consensus hit/miss plot                              #
#                    (Revision: includes prenatal and replication networks as well                 #
####################################################################################################


library(limma)      
library(gtools)
library(stringr)
library(ggplot2)
library(WGCNA)
library(ggrepel)
library(patchwork)
library(magrittr)
library(gridtext)
library(ggh4x)       ##For the nested facets
library(colorspace)  ##For color scales
library(prismatic)   ## For simple color manipulation
library(tidyverse)

#A mapping of networks names standardized for the manuscript
new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\old-to-new Network names all.csv") %>% tibble::deframe()

#Bins of extensions applicable for SCZ enrichments
PGC.bins = c("PGC", "kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")

#Read list of module previously identified in publications for SCZ
priortised_modules = list(
  "Fromer2016_case"   = c("yellow","red","blue","lightyellow","greenyellow","cyan","grey60")              ,              
  "Fromer2016_control"= c("magenta","blue","greenyellow","salmon","turquoise","bisque")                   ,
  "Radulescu2020"     = c("blue","brown","green")                                                         , 
  "Walker2019"        = c("red","blue")                                                                   ,                                                                   
  "Werling2020"       = c("yellow","red","brown","magenta","cyan")                                        ,                                        
  "Li2018"            = c("green","skyblue3","blue","violet","royalblue","black")                         ,                         
  "Gandal2018PE"      = c("lightyellow","brown","red","black","blue","darkred","saddlebrown","turquoise") , 
  "Gandal2018PE_cs"   = c("pink","lightcyan")                                                             ,                                                             
  "Gandal2018a"        = c("turquoise","green","yellow","salmon","purple")                                 ,
  "Gandal2018b"      = c("lightyellow","brown","red","black","blue","darkred","saddlebrown","turquoise") , 
  "Gandal2018b_cs"   = c("pink","lightcyan")                                                             ,                                                             
  "Gandal2018"        = c("turquoise","green","yellow","salmon","purple")                                 ,
  "Pergola2017"       = c(NA)                                                                             ,
  "Pergola2019"       = c("darkgreen")                                                                    ,
  "Pergola2020"       = c("darkorange")                                                                   ,
  "BRNACC"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNAMY"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNCBH"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNCBL"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNCDT"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNCTX"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNCTXB24"         = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNCTXBA9"         = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNHIP"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNHYP"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNPUT"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNSNA"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "ALL"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BGA"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BROD"              = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "CEREBELLUM"        = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "CTX"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "NS.SCTX"           = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "STR"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "SUBCTX"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "WHOLE_BRAIN"       = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                 
)

#Networks generated in this manuscript
our_nets     = c("CN","CN Juvenile", "CN Adult", "CN Older Adult",
                 "DLPFC", "DLPFC Perinatal","DLPFC Juvenile", "DLPFC Adult", "DLPFC Older Adult",
                 "Replication Perinatal","Replication Juvenile",
                 "HP", "HP Perinatal", "HP Juvenile", "HP Adult", "HP Older Adult")
#Sample matched DG-HP networks
matched_nets = c("HP.QSVA","DG.noQSVA")
#Previously published networks considered in this manuscript
other_nets   = sort(c("Pergola2017","Pergola2019","Pergola2020",
                    "Fromer2016_case", "Fromer2016_control", "Gandal2018a", "Gandal2018b", "Gandal2018b_cs",
                    "BRNCTX","Radulescu2020","Walker2019","Werling2020"))
#Only Hartl2021_BRNCTX network was used in this manuscript for the downstream analysis
Hartl_nets = sort(c("BRNACC", "BRNAMY", "BRNCBH", "BRNCBL", "BRNCDT", "BRNCTX", "BRNCTXB24",
               "BRNCTXBA9","BRNHIP", "BRNHYP", "BRNPUT", "BRNSNA",  "BGA",    "STR",
               "BROD",   "CTX",    "SUBCTX", "NS.SCTX", "CEREBELLUM",    "ALL",    "WHOLE_BRAIN"))

# #Get the order of tissues for figures
# complete_tissue_order = c("CN","CN Juvenile", "CN Adult", "CN Older Adult",
#                           "DLPFC", "DLPFC Perinatal", "DLPFC Juvenile", "DLPFC Adult", "DLPFC Older Adult",
#                           "Replication Perinatal","Replication Juvenile",
#                           "HP", "HP Perinatal", "HP Juvenile", "HP Adult", "HP Older Adult",
#                           "HP.QSVA","DG.noQSVA",
#                           "Pergola2017","Pergola2019","Pergola2020",
#                           "Fromer2016_case", "Fromer2016_control", "Gandal2018a", "Gandal2018b", "Gandal2018b_cs",
#                           "BRNCTX","Radulescu2020","Walker2019", "Werling2020"
#                           )

#Read the long_form enrichments results (new with metanets and prenatal networks)
p_all = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\all.files_long_metanets_test.rds")

#Hartl2021 gene-module list
hm <- readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Hartl2021\\Hartl2021-gene-module-list.rds")
hm = hm %>% map(~{
  names(.x) = gsub("WHOLE_BRAIN\\.","BW.",names(.x))
  .x
})
hm = unique(unlist(map(hm,names)))
hm = gtools::mixedsort(hm)

#For Hartl2021 networks only, add psuedo colors for modules, for all other networks, module name is module color.
module.color = WGCNA::standardColors(length(hm)) %>% set_names(hm)
module.color["grey"] = "grey"

p_all$module.color[p_all$new_tissue_age %in% Hartl_nets] = module.color[p_all$module.color[p_all$new_tissue_age %in% Hartl_nets]]
p_all = p_all %>% mutate(
  modules = gsub("WHOLE_BRAIN\\.","BW.",modules),
  sig.modules = ifelse(-log10(bonf.vals)> -log10(0.05), modules,NA),
  tissue_age_sig.modules = ifelse(-log10(bonf.vals)> -log10(0.05), paste0(new_tissue_age,">",modules),NA))

#For panel A only consider SCZ enrichments
scz = p_all %>%
  filter(list.type == "PGC3" &
           !(modules %in% c("grey")) &  #Remove grey modules
           bins %in% PGC.bins &         #Only keeping PCG bins
           !(new_tissue_age %in% c("DG","HP.noQSVA","DG.QSVA"))) %>%
  mutate( 
    new_tissue_age = factor(new_tissue_age, levels = c("CN", "CN Juvenile", "CN Adult", "CN Older Adult",
                                                       "DLPFC", "DLPFC Prenatal", "DLPFC Perinatal", "Replication Perinatal","DLPFC Juvenile", "Replication Juvenile","DLPFC Adult", "DLPFC Older Adult",
                                                       "HP","HP Prenatal", "HP Perinatal", "HP Juvenile", "HP Adult", "HP Older Adult",
                                                       "HP.QSVA", "DG.noQSVA", "Fromer2016_case", "Fromer2016_control", "Gandal2018a", "Gandal2018b", "Gandal2018b_cs",
                                                       "Li2018", "Pergola2017", "Pergola2019", "Pergola2020", "Radulescu2020", "Walker2019", "Werling2020",
                                                       "ALL", "BGA", "BRNACC", "BRNAMY", "BRNCBH", "BRNCBL", "BRNCDT", "BRNCTX", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT","BRNSNA", "BROD", "CEREBELLUM", "CTX", "NS.SCTX", "STR", "SUBCTX", "WHOLE_BRAIN")))                              


scz = scz %>% rowwise() %>%
  mutate(priortised = factor(ifelse(modules %in% priortised_modules[[as.character(new_tissue_age)]],"Previously reported","Novel Report"),levels = c("Previously reported","Novel Report")),
         significant = ifelse(-log10(bonf.vals) > -log10(0.05),T,F)
  )

scz = scz %>% group_by(list.type, biotype, new_tissue_age, modules) %>% 
  mutate(sig.windows = sum(bonf.vals<0.05),
         everSig = ifelse(sig.windows>=1,T,F))


p = list()

#First make plot for each network separately, then arrange them using Patchwork
for (net in unique(scz$new_tissue_age)){
  pos = position_jitter(seed = 121,height = 0, width = 0.3)

  p[[net]] = scz %>% filter(new_tissue_age == net) %>% mutate(new_tissue_age = factor(gsub("BRNCTX$","Hartl2021_BRNCTX",new_tissue_age),levels = sort(gsub("BRNCTX$","Hartl2021_BRNCTX",levels(new_tissue_age)))), bins = factor(bins,levels = PGC.bins)) %>% 
    ggplot(aes(bins,-log10(bonf.vals), fill = I(module.color), color = significant, size = significant, shape = priortised)) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) + 
     geom_jitter(aes(alpha = significant), color = grDevices::adjustcolor("black",alpha.f = 0.5),position = pos, show.legend = F)+
     {if(net %in% Hartl_nets){ #For Hartl networks only write module names alongwith pseudo-colored modules
       geom_text_repel(aes(label = sig.modules,color = I(module.color)),alpha = 0.75, size = 3.5, direction = "y", min.segment.length = 10, force_pull = 1.4, force = 1.4,max.overlaps = 15, nudge_y = 0.7)
     }} +
    scale_size_manual( values = c(1.8,4.2))+guides(size = "none")+
    scale_alpha_manual(values = c(0.7,0.8))+guides(alpha = "none")+
    scale_shape_manual(guide = guide_legend(override.aes = list(size = 4.2, fill = "black", color = "black")),labels = c("Previously reported","Novel Report"), values = c(23,21), drop = F)+
    labs(shape = "")+
    scale_x_discrete(labels = c("PGC","kb_000" ,"kb_020","kb_050" ,"kb_100","kb_150" ,"kb_200","kb_250" ,"kb_500"))+
    scale_y_continuous(limits = c(0,15), breaks = c(0,5,10,15))+
    facet_wrap(new_tissue_age~.) + 
    ylab(label = NULL) + 
    xlab(label = NULL) +
    ggtitle(label = "")+
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_blank(),
          panel.grid.major.x = element_line(size = 0.2),
          panel.grid.minor.x = element_line(size = 0.),
          panel.grid.major.y = element_line(size = 0.3),
          panel.grid.minor.y = element_line(size = 0.3),
          panel.border    = element_rect(size = 0.1, color = prismatic::clr_lighten("black")),
          plot.margin = margin(0.5,0.5,0.5,0.5), #margin(0,0,0,0,unit = "pt"),
          legend.position = c(0,"Inf"),#"bottom",
          legend.margin=margin(-5,-5,-5,-5),
          legend.box.margin=margin(-5,-5,-5,-5),
          strip.text = element_text(margin = margin(0.3,0.3,0.3,0.3, "lines"), size = 14),
          strip.background = element_rect(size = 0.21, fill = "darkgrey", color = NA)
          
    )
  
  
  if (net %in% c("CTX","Pergola2020","HP.QSVA","HP")) { #Keep axis text for bottom left networks only
    p[[net]] = p[[net]] +
      # xlab(label = "PGC based reference bins")+
      # ylab(label = "-log10(bonf)") +
      theme(axis.text.x = element_text(angle = 90, size = 16, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 16),
            axis.ticks.x = element_line(size = 0.5, color = "white"),
            axis.ticks.length.x = unit(3,"pt"),
            axis.ticks.y = element_line(size = 0.5),
            axis.ticks.length.y = unit(3,"pt"),
            #axis.ticks.margin=unit(c(5,2),'pt')
            axis.title.x = element_text(size = 15, vjust = 0),
            axis.title.y = element_text(size = 15, vjust = 1)
      ) 
    
  }
}

#Correct the order of plots before applying design
p = p[c("CN", "CN Juvenile", "CN Adult", "CN Older Adult",
        "DLPFC", "DLPFC Perinatal", "DLPFC Juvenile", "DLPFC Adult", "DLPFC Older Adult",
        "Replication Perinatal","Replication Juvenile",
        "HP","HP Perinatal", "HP Juvenile", "HP Adult", "HP Older Adult",
        "HP.QSVA", "DG.noQSVA", "Fromer2016_case", "Fromer2016_control", "Gandal2018a", "Gandal2018b", "Gandal2018b_cs",
        "BRNCTX", "Pergola2017", "Pergola2019", "Pergola2020", "Radulescu2020", "Walker2019", "Werling2020",
        "ALL", "BGA", "BRNACC", "BRNAMY", "BRNCBH", "BRNCBL", "BRNCDT", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT","BRNSNA", "BROD", "CEREBELLUM", "CTX", "NS.SCTX", "STR", "SUBCTX","WHOLE_BRAIN")]

p_ours      = p[our_nets]
p_ours_AB   = map(p_ours,~ .x %+% subset(.x$data, biotype =="all.biotypes"))
p_ours_PC   = map(p_ours,~ .x %+% subset(.x$data, biotype =="protein.coding"))
p_others    = p[c("Fromer2016_case","Fromer2016_control","Gandal2018a","Gandal2018b","Gandal2018b_cs","BRNCTX","Pergola2017","Pergola2019","Pergola2020","Radulescu2020","Walker2019","Werling2020")]#;p[other_nets]
p_others_AB = map(p_others,~ .x %+% subset(.x$data, biotype =="all.biotypes"))
p_others_PC = map(p_others,~ .x %+% subset(.x$data, biotype =="protein.coding"))
p_matched     = p[matched_nets]
p_matched_AB  = map(p_matched,~ .x %+% subset(.x$data, biotype =="all.biotypes"))
p_matched_PC  = map(p_matched,~ .x %+% subset(.x$data, biotype =="protein.coding"))
p_Hartl     = p[Hartl_nets %>% .[.!= "BRNCTX"]]
p_Hartl_AB  = map(p_Hartl,~ .x %+% subset(.x$data, biotype =="all.biotypes"))
p_Hartl_PC  = map(p_Hartl,~ .x %+% subset(.x$data, biotype =="protein.coding"))

design_ours = c(
  area(t = 1,l =  1),                    area(t = 1,l =  3),area(t = 1,l =  4),area(t = 1,l =  5),  #"CN","CN Juvenile", "CN Adult", "CN Older Adult",
  area(t = 2,l =  1),area(t = 2,l =  2), area(t = 2,l =  3),area(t = 2,l =  4),area(t = 2,l =  5),  #"DLPFC","DLPFC Perinatal",  "DLPFC Juvenile", "DLPFC Adult", "DLPFC Older Adult",
                     area(t = 3,l =  2), area(t = 3,l =  3),                                        #"Replication Perinatal", "Replication juvenile"
  area(t = 4,l =  1),area(t = 4,l =  2), area(t = 4,l =  3),area(t = 4,l =  4),area(t = 4,l =  5)  #"HP","HP Perinatal", "HP Juvenile", "HP Adult", "HP Older Adult",
)


design_others = c(
  area(t = 1,l =  1), area(t = 1,l =  2),area(t = 1,l =  3),area(t = 1,l =  4),  #"Fromer2016_case", "Fromer2016_control", "Gandal2018", "Gandal2018PE"
  area(t = 2,l =  1), area(t = 2,l =  2),area(t = 2,l =  3),area(t = 2,l =  4),  #"Gandal2018PE_cs", "BRNCTX", "Pergola2017", "Pergola2019" 
  area(t = 3,l =  1), area(t = 3,l =  2),area(t = 3,l =  3),area(t = 3,l =  4)   #"Pergola2020", "Radulescu2020", "Walker2019", "Werling2020"
)

design_matched = c(
  area(t = 1,l =  1), area(t = 1,l =  2),area(t = 1,l =  3)  #"HP.QSVA","DG.noQSVA",blank
)


design_Hartl = c(
  area(t = 1,l =  1), area(t = 1,l =  2), area(t = 1,l =  3), area(t = 1,l =  4), area(t = 1,l =  5),
  area(t = 2,l =  1), area(t = 2,l =  2), area(t = 2,l =  3), area(t = 2,l =  4), area(t = 2,l =  5),
  area(t = 3,l =  1), area(t = 3,l =  2), area(t = 3,l =  3), area(t = 3,l =  4), area(t = 3,l =  5),
  area(t = 4,l =  1), area(t = 4,l =  2), area(t = 4,l =  3), area(t = 4,l =  4), area(t = 4,l =  5)
)

#Saving Panel A for figures 2 and 5
setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Figures")
pp_ours_AB = wrap_plots(p_ours_AB, design = design_ours, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(WGCNA,AB).svg"            ,width = 14, height = 7)
pp_ours_AB#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_ours_PC = wrap_plots(p_ours_PC, design = design_ours, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(WGCNA,PC).svg"            ,width = 14, height = 7.5)
pp_ours_PC#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_others_AB = wrap_plots(p_others_AB, design = design_others, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(Published,AB).svg"        ,width = 14, height = 7.5)
pp_others_AB#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_others_PC = wrap_plots(p_others_PC, design = design_others, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(Published,PC).svg"        ,width = 14, height = 7)
pp_others_PC#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_matched_AB = wrap_plots(p_matched_AB, design = design_matched, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(Matched,AB).svg"          ,width = 12, height = 4)
pp_matched_AB#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_matched_PC = wrap_plots(p_matched_PC, design = design_matched, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(Matched,PC).svg"          ,width = 12, height = 4)
pp_matched_PC#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_Hartl_AB = wrap_plots(p_Hartl_AB, design = design_Hartl, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(Hartl,AB).svg"          ,width = 16, height = 8)
pp_Hartl_AB#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()

pp_Hartl_PC = wrap_plots(p_Hartl_PC, design = design_Hartl, guides = "collect", heights = 5) + 
  plot_annotation(title = "Schizophrenia enrichment", theme = theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 20), legend.position = "bottom", plot.margin = margin(2,2,0,5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(-15,10,30,10)))
svg(filename = "panelwise(Hartl,PC).svg"          ,width = 16, height = 8)
pp_Hartl_PC#;grid::grid.draw(grid::textGrob("-log10(bonf)", x = 0.01, y =0.46, rot = 90, hjust = 0,  gp = gpar(fontsize = 20)));#grid::grid.draw(grid::textGrob("PGC3 based reference bins", x = 0.14, y = 0.11, rot = 0, hjust = 0.5, vjust = 0,  gp = gpar(fontsize = 21)))
dev.off()



###########################
###########################
###########################
####All Enrichments----####
mega.risk.modules = readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Hartl2021\\new_risk_modules (with Hartl).rds")
mega.risk.modules$risk.modules = append(mega.risk.modules$risk.modules[!grepl("Walker2019",mega.risk.modules$risk.modules)],c("Replication Juvenile>blue","Replication Perinatal>darkred","Replication Perinatal>red"))
mega.risk.modules$risk.modules.noGrey = append(mega.risk.modules$risk.modules.noGrey[!grepl("Walker2019",mega.risk.modules$risk.modules.noGrey)],c("Replication Juvenile>blue","Replication Perinatal>darkred","Replication Perinatal>red"))

tissue_order = c("CN","DLPFC", "HP",
                 "CN Juvenile", "CN Adult", "CN Older Adult",
                 "DLPFC Perinatal", "Replication Perinatal","DLPFC Juvenile",'Replication Juvenile', "DLPFC Adult", "DLPFC Older Adult",
                 "HP Perinatal", "HP Juvenile", "HP Adult", "HP Older Adult",
                 "HP.QSVA","DG.noQSVA",
                 "Fromer2016_case", "Fromer2016_control", "Gandal2018a", "Gandal2018b", "Gandal2018b_cs","BRNCTX",
                 "Pergola2017","Pergola2019","Pergola2020",
                 "Radulescu2020","Walker2019", "Werling2020",
                 "ALL", "BGA", "BRNACC", "BRNAMY", "BRNCBH", "BRNCBL", "BRNCDT", "BRNCTXB24", "BRNCTXBA9", "BRNHIP", "BRNHYP", "BRNPUT","BRNSNA", "BROD", "CEREBELLUM", "CTX", "NS.SCTX", "STR", "SUBCTX","WHOLE_BRAIN")

tissue_class = c("non-parsed","non-parsed","non-parsed",
                 "parsed","parsed","parsed",
                 "parsed","parsed","parsed","parsed","parsed","parsed",
                 "parsed","parsed","parsed","parsed",
                 "matched","matched",
                 "published","published","published","published","published",
                 "published","published","published",
                 "published","published","published","published",
                 rep("Hartl",20))

tissue_vec = tissue_class %>% set_names(tissue_order)


priortised_modules = list(
  "Fromer2016_case"   = c("yellow","red","blue","lightyellow","greenyellow","cyan","grey60")              ,              
  "Fromer2016_control"= c("magenta","blue","greenyellow","salmon","turquoise","bisque")                   ,
  "Radulescu2020"     = c("blue","brown","green")                                                         , 
  "Walker2019"        = c("red","blue")                                                                   ,                                                                   
  "Werling2020"       = c("yellow","red","brown","magenta","cyan")                                        ,                                        
  "Li2018"            = c("green","skyblue3","blue","violet","royalblue","black")                         ,                         
  "Gandal2018PE"      = c("lightyellow","brown","red","black","blue","darkred","saddlebrown","turquoise") , 
  "Gandal2018PE_cs"   = c("pink","lightcyan")                                                             ,                                                             
  "Gandal2018a"        = c("turquoise","green","yellow","salmon","purple")                                 ,
  "Gandal2018b"      = c("lightyellow","brown","red","black","blue","darkred","saddlebrown","turquoise") , 
  "Gandal2018b_cs"   = c("pink","lightcyan")                                                             ,                                                             
  "Gandal2018"        = c("turquoise","green","yellow","salmon","purple")                                 ,
  "Pergola2017"       = c(NA)                                                                             ,
  "Pergola2019"       = c("darkgreen")                                                                    ,
  "Pergola2020"       = c("darkorange")                                                                   ,
  "BRNACC"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNAMY"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNCBH"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNCBL"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNCDT"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNCTX"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNCTXB24"         = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNCTXBA9"         = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNHIP"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNHYP"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BRNPUT"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,                                                                    
  "BRNSNA"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "ALL"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BGA"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "BROD"              = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "CEREBELLUM"        = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "CTX"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "NS.SCTX"           = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "STR"               = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "SUBCTX"            = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                           ,
  "WHOLE_BRAIN"       = c("BW.M1", "BW.M4","CEREB.M3","CTX.M3","WHOLE_BRAIN.M1","WHOLE_BRAIN.M4")                                 
) %>% tibble::enframe() %>% unnest_longer(value) %>% unite(name,value,col = "ID", sep = ">") %>% pull(ID)


p_all$enrichment[p_all$enrichment == "pleiotropic_Lof_Denovo_CNVs"] = "LOF"
p_all$list.type [p_all$list.type  == "pleiotropic_Lof_Denovo_CNVs"] = "LOF"
p_all$bins      [p_all$bins       == "LoF"]                         = "LOF"

names(p_all$bonf.vals) = p_all$bins
p_all = p_all %>% mutate(new_tissue_age = as.character(new_tissue_age), biotype =case_when(biotype == "all.biotypes"~"AB",biotype=="protein.coding"~"PC")) %>% relocate(new_tissue_age, .after = tissue_age)


#Filter all enrichment files to keep only relevant bins depending upon the enrichment type 
p_all = p_all %>% mutate(bins_to_keep = ifelse(
  (list.type %in% c("PGC3","ad","adhd","als","asd","bip","cd","mdd","ms","ocd","pd","ptsd","ra","sa","uc") & bins %in% PGC.bins) |
    (list.type %in% c("PGC3.negative")                                                                       & bins %in% PGC.bins)        |
    (list.type %in% c("DEGs")                                                                                & bins %in% c("DEGs_human_Sousa", "DEGs SCZ"))                                                        | ##c("Apua_CAUDATE_sczd", "Clozapine_mouse", "DEGs_human_Sousa", "Fromer_nosva", "Haloperidol", "Jaffe_DLPFC_HIPPO_sczd", "Jaffe_DLPFC_sczd", "Jaffe_HIPPO_sczd", "Jaffe_sczd_2018", "Risperidone_aripripazole_blood")
    (list.type %in% c("DMGs")                                                                                & bins %in% c("DMGs general"))                                                                        | ##c("DMGs_Hannon", "DMGs_Jaffe", "DMGs_Kinoshita", "DMGs_Montano", "DMGs_Numata", "DMGs_Wockner")
    (list.type %in% c("Druggable_genes")                                                                     & bins %in% c("Drug general"))                                                                        | ##c("Finan", "IDG_Schizo_Tclin", "Santos", "Santos_Tclin","Wang")
    (list.type %in% c("TWAS")                                                                                & bins %in% c("TWAS general"))                                                                        | ##c("Gandal","Gusev", "Jaffe_DLPFC_pgc2.clozuk", "Jaffe_HIPPO_pgc2.clozuk", "Jaffe_DLPFC_pgc2", "Jaffe_HIPPO_pgc2", "Hall","Paquola_CAUDATE")
    (list.type %in% c("LOF")                                                                                 & bins %in% c("LOF"))                                                                                 , ##c("pleiotropic", "LoF", "Denovo", "CNVs")
  T, #Keep
  F  #Remove
)) %>% filter(bins_to_keep)

#Convert bonferroni corrected p-values to -log10(bonf.vals)
p_all$final.vals = -log10(p_all$bonf.vals)

##Load additional pre-computed enrichments NEW: Cell Specificity, GO (New) [Adds enrichment results for replication networks]
results = new.env()
load("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Paper new analysis/Revision/metanets/all_results(metanets).RData", envir = results)
results$res_AB_CellSpecificity = map_dfr(.id = "tissue_age",results$res_AB_CellSpecificity,~{ .x %>% as.data.frame() %>%  tibble::rownames_to_column("modules") %>% dplyr::select(c(modules,contains("bonf..")))}) %>% mutate(across(contains("bonf.."),~-log10(.x))) %>% set_rownames(paste0(.$tissue_age,">",.$modules)) %>% set_colnames(gsub("bonf..","",colnames(.)))

##Load additional pre-computed enrichments OLD: Cell Specificity, MAGMA, HMAGMA, GO, TF, binwiseMAGMA, PGC3.based.df = "PGC3 LOCI based AB", snp.based.df = "PGC3.loci.windows, top.gwas = "PGC3.perm.AB" (Old)
ne = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Shared Data\\Giulio's ML code\\Enrichments\\grch38[PGC125new]\\new_helper_data_grch38(with Hartl).RData", envir = ne)

#Rename networks and modules names for consistency
ne = map(as.list(ne), ~{ 
  xx = .x %>% set_rownames(gsub("\n"," ",rownames(.)))
  rownames(xx) = gsub("Gandal2018>"     ,"Gandal2018a>"   ,rownames(xx))
  rownames(xx) = gsub("Gandal2018PE>"   ,"Gandal2018b>"   ,rownames(xx))
  rownames(xx) = gsub("Gandal2018PE_cs>","Gandal2018b_cs>",rownames(xx))
  rownames(xx) = gsub("Gandal2018PE_cs>","Gandal2018b_cs>",rownames(xx))
  rownames(xx) = gsub("WHOLE_BRAIN\\."  ,"BW."            ,rownames(xx))
  xx
})                                

##Prepping enrichment names for final figures, keeping only relevant bins/columns for each enrichment
ne$MAGMA.df  = ne$MAGMA.df [,grep("tissue_age|modules|35.10kbp\\.PC"                 ,colnames(ne$MAGMA.df))]
ne$HMAGMA.df = ne$HMAGMA.df[,grep("tissue_age|modules|Adult_brain|Fetal_brain\\.PC"  ,colnames(ne$HMAGMA.df))]
colnames(ne$MAGMA.df)         = gsub("35.10kbp.PC","MAGMA.kb_35.10.PC",colnames(ne$MAGMA.df))
colnames(ne$HMAGMA.df)[3:4]   = c("H-MAGMA.Adult_brain.PC","H-MAGMA.Fetal_brain.PC")#paste0("HMAGMA.",colnames(ne$HMAGMA.df)[3])
colnames(ne$PGC3.based.df)[3] = "PGC3.loci.based" 
colnames(ne$snp.based.df)[3]  = "PGC3.loci.window+count" 
colnames(ne$top.gwas.df)[3]   = "PGC3.perm.AB+count"
ne$snp.based.df$modules[ne$snp.based.df$tissue_age == "BRNCTX"] = toupper(ne$snp.based.df$modules[ne$snp.based.df$tissue_age == "BRNCTX"])
ne$snp.based.df$modules = gsub("GREY","grey",ne$snp.based.df$modules)
rownames(ne$snp.based.df)[ne$snp.based.df$tissue_age == "BRNCTX"] = paste0(ne$snp.based.df$tissue_age[ne$snp.based.df$tissue_age == "BRNCTX"],">",ne$snp.based.df$modules[ne$snp.based.df$tissue_age == "BRNCTX"])
ne$top.gwas.df$modules[ne$top.gwas.df$tissue_age == "BRNCTX"] = toupper(ne$top.gwas.df$modules[ne$top.gwas.df$tissue_age == "BRNCTX"])
ne$top.gwas.df$modules = gsub("GREY","grey",ne$top.gwas.df$modules)
rownames(ne$top.gwas.df)[ne$top.gwas.df$tissue_age == "BRNCTX"] = paste0(ne$top.gwas.df$tissue_age[ne$top.gwas.df$tissue_age == "BRNCTX"],">",ne$top.gwas.df$modules[ne$top.gwas.df$tissue_age == "BRNCTX"])

#Where enrichment results are not available for replication network, add rows with NA values
ne$PGC3.based.df = rbind(ne$PGC3.based.df, data.frame(tissue_age = c("Replication Perinatal","Replication Perinatal", "Replication Juvenile"), modules = c("red","darkred","blue"), PGC3.loci.based = NA) %>% set_rownames(c("Replication Perinatal>red","Replication Perinatal>darkred","Replication Juvenile>blue")))
ne$HMAGMA.df     = rbind(ne$HMAGMA.df    , data.frame(tissue_age = c("Replication Perinatal","Replication Perinatal", "Replication Juvenile"), modules = c("red","darkred","blue"), `H-MAGMA.Adult_brain.PC` = NA, `H-MAGMA.Fetal_brain.PC` = NA, `Adult_brain.ALL` = NA,check.names = F) %>% set_rownames(c("Replication Perinatal>red","Replication Perinatal>darkred","Replication Juvenile>blue")) )
ne$MAGMA.df      = rbind(ne$MAGMA.df     , data.frame(tissue_age = c("Replication Perinatal","Replication Perinatal", "Replication Juvenile"), modules = c("red","darkred","blue"), "MAGMA.kb_35.10.PC" = c(0,0,3.48), check.names = F) %>% set_rownames(c("Replication Perinatal>red","Replication Perinatal>darkred","Replication Juvenile>blue")) )
ne$top.gwas.df   = rbind(ne$top.gwas.df  , data.frame(tissue_age = c("Replication Perinatal","Replication Perinatal", "Replication Juvenile"), modules = c("red","darkred","blue"), "PGC3.perm.AB+count" = NA, check.names = F) %>% set_rownames(c("Replication Perinatal>red","Replication Perinatal>darkred","Replication Juvenile>blue")) )

#Make a long form df of these extra enrichments 
extra_enrichments = ne[c("PGC3.based.df","HMAGMA.df","snp.based.df","top.gwas.df","MAGMA.df")] %>% purrr::reduce(left_join) %>%
  pivot_longer(-c(tissue_age,modules), names_to = "list.type", values_to = "final.vals") %>% mutate(tissue_age = gsub("\n"," ",tissue_age), enrichment = list.type, bins = list.type, ID = paste0(tissue_age,">",modules))
extra_enrichments$ID = gsub("Gandal2018>"     ,"Gandal2018a>"   ,extra_enrichments$ID)
extra_enrichments$ID = gsub("Gandal2018PE>"   ,"Gandal2018b>"   ,extra_enrichments$ID)
extra_enrichments$ID = gsub("Gandal2018PE_cs>","Gandal2018b_cs>",extra_enrichments$ID)
extra_enrichments$ID = gsub("Gandal2018PE_cs>","Gandal2018b_cs>",extra_enrichments$ID)
extra_enrichments$ID = gsub("WHOLE_BRAIN\\."  ,"BW."            ,extra_enrichments$ID)

#For long form df from enrichment results calculated earlier and the extra enrichment results
p_all1 = full_join(p_all, results$res_AB_CellSpecificity %>% pivot_longer(cols = c(-tissue_age,-modules), values_to = "final.vals", names_to = "bins") %>% mutate(list.type = "Cell", enrichment = "Cell")) %>% 
  mutate(ID = paste0(tissue_age,">",modules)) %>% filter(ID %in% mega.risk.modules$risk.modules.noGrey)
#For pathology enrichments, get the count of significant windows (number of PC bins with bonf.vals <0.05)
p_patho_only = p_all1 %>% filter(ID %in% mega.risk.modules$risk.modules.noGrey & biotype %in% "AB" & list.type %in% c("ad","adhd","als","asd","bip","cd","mdd","ms","ocd","pd","ptsd","ra","sa","uc") & bins %in% c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500")) %>% 
  group_by(list.type,ID) %>% summarise(n.sig.windows = sum(final.vals > -log10(0.05)))


p_all2 = full_join(p_all1, extra_enrichments) %>% filter(ID %in% mega.risk.modules$risk.modules.noGrey)

#Make ggplot heatmap of -log10(bonf.vals)
pl = p_all2 %>% mutate(`-log10(bonf.vals)` = ifelse(final.vals > -log10(0.05), final.vals, 0)) %>% 
  ggplot(aes(bins, ID, fill = `-log10(bonf.vals)`)) + 
  geom_tile(width = 1, height = rel(200), color = prismatic::clr_lighten("black"), size = 0.2, position = "identity")+
  #geom_text(data = extra_enrichments %>% mutate(`-log10(bonf.vals)` = ifelse(final.vals > -log10(0.05), final.vals, 0)) %>% filter(bins %in% c("PGC3.loci.window+count","PGC3.perm.AB+count")),aes(label = final.vals))+
  facet_grid(ID~.,scale = "free", switch = "y")+ 
  scale_x_discrete(position = "top")+
  theme_minimal() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.text.x.top = element_text(angle = 46, vjust = 0, hjust = 0, size = 15),
    axis.title = element_blank(),
    axis.ticks.length = unit(0,"pt"),
    panel.spacing = unit(2,"pt")
  ) +
  guides(fill  = guide_legend(ncol = 1, title = paste0("-log10(bonf)"), title.position = "top", title.theme = element_text(face = "bold", hjust = 0.2), label.position = "right", label.theme = element_text(angle = 0, hjust = 0.5), keywidth = unit(20,"pt"), keyheight = unit(20,"pt")),
         color = guide_legend(ncol = 1, title = paste0("-log10(bonf)"), title.position = "top", title.theme = element_text(face = "bold", hjust = 0.2), label.position = "right", label.theme = element_text(angle = 0, hjust = 0.5), keywidth = unit(20,"pt"), keyheight = unit(20,"pt")),
         size  = "none")

#Make ggplot heatmap of count values (pathologies)
pl.patho = p_patho_only %>% 
  ggplot(aes(list.type, ID, fill = n.sig.windows)) + 
  geom_tile(width = rel(1), height = rel(200), color = prismatic::clr_lighten("black"), size = 0.2, position = "identity")+
  geom_text(aes(label = n.sig.windows))+  
  facet_grid(ID~.,scale = "free", switch = "y")+ 
  scale_x_discrete(position = "top")+
  theme_minimal() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.text.x.top = element_text(angle = 46, vjust = 0, hjust = 0, size = 15),
    axis.title = element_blank(),
    axis.ticks.length = unit(0,"pt"),
    panel.spacing = unit(2,"pt")
  ) 

#Plotting following enrichments for all SCZ risk modules: SCZ, Cell specificity, DEGs, DMGs, Druggable Genes, LOF, TWAS
gp = list(
  pl       %+% subset(pl$data       %>% mutate(bins = factor(bins, levels = c("PGC","kbp_0","kbp_20","kbp_50","kbp_100","kbp_150","kbp_200","kbp_250","kbp_500"))), list.type %in%  "PGC3" & biotype %in% "PC") + scale_fill_gradient(low = prismatic::clr_lighten("white"), high = "#4B0092", limits = c(0,10), oob =scales::oob_squish, na.value = 'grey')                          + theme(strip.background = element_rect(color = "lightgrey", fill = "lightgrey")),  #SCZ bins
  pl       %+% subset(pl$data       %>% mutate(bins = factor(bins, levels = c("ASC","END","exCA","exDG","exPFC","GABA","MG","NSC","ODC","OPC")))                  , list.type %in% "Cell")                      + scale_fill_gradient(low = prismatic::clr_lighten("white"), high = "#2E8B57", limits = c(0,15), oob =scales::oob_squish, na.value = 'grey')                          + theme(strip.background = element_blank(), strip.text.y.left =  element_blank()),  #Cell specificity bins
  pl       %+% subset(pl$data       %>% mutate(bins = factor(bins, levels = c("DEGs SCZ","DEGs_human_Sousa","DMGs general","Drug general","LOF","TWAS general"))) , !is.na(bins))                               + scale_fill_gradient(low = prismatic::clr_lighten("white"), high = "#2E8B57", limits = c(0,15), oob =scales::oob_squish, na.value = 'grey')                          + theme(strip.background = element_blank(), strip.text.y.left = element_blank()),    #Other enrichments
  pl       %+% subset(pl$data       %>% mutate(bins = factor(bins, levels = c("PGC3.loci.based","H-MAGMA.Adult_brain.PC","H-MAGMA.Fetal_brain.PC","Adult_brain.ALL","MAGMA.kb_35.10.PC"))) , !is.na(bins))      + scale_fill_gradient(low = prismatic::clr_lighten("white"), high = "#2E8B57", limits = c(0,15), oob =scales::oob_squish, na.value = 'grey') + theme(strip.background = element_blank(), strip.text.y.left = element_blank()),    #Other enrichments (-log10(bonf))
  pl       %+% subset(pl$data       %>% mutate(bins = factor(bins, levels = c("PGC3.loci.window+count","PGC3.perm.AB+count"))) , !is.na(bins))                                                                  + scale_fill_binned_sequential(limits = c(0,9), breaks = c(0:9), na.value = 'grey', guide = "none") + theme(strip.background = element_blank(), strip.text.y.left = element_blank()) + guides(fill ="none"),    #Other enrichments (count)
  pl.patho + scale_fill_binned_sequential(limits = c(0,9), breaks = c(0:9)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank())        #Pathology enrichments (number of significant windows count)
) %>% wrap_plots(guides = "collect",widths = c(9,10,6,5,2,14),ncol=6)
#ggsave(gp          , filename = "enrichments_plots_all.svg"          ,device = svg, width = 17, height = 12)

# #####################
# ######################
# #### Go plot --------
library(ggh4x)
#Format GO results into a df
GO.df = map_dfr(.id= "enrich_type",transpose(results$res_AB_GO),~{
  ont = .x
  map_dfr(.id = "network",ont,~{ .x@compareClusterResult} )
}) %>% mutate(enrich_type = strsplit2(enrich_type,"_")[,4],
              risk.module = paste0(network,">",Cluster), .before = network)

#Filter to identified risk modules only
GO.df = GO.df %>% filter(risk.module %in% mega.risk.modules$risk.modules.noGrey & p.adjust <0.05)
#Further filter to keep only top-3 significant processes per module per enrichment type (BP/MF/CC/KEGG/Reactome)  
GO.top3 = GO.df %>% group_by(enrich_type,risk.module) %>% slice_min(p.adjust,n=3) %>% split(.,.$enrich_type)

#GO ID ordered by frequency
GO.top3.ordered = imap(GO.top3 ,~{
  GO.df %>% filter(ID %in% unique(.x$ID)) %>% group_by(ID,Description) %>% dplyr::count() %>% arrange(desc(n)) %>% pull(ID)
})


GO.plot = GO.df %>% 
  filter(ID %in% unlist(GO.top3.ordered)  & enrich_type %in% c("BP","KEGG","Pathway")) %>% add_count(ID,name = "n1") %>% filter(n1>4) %>% add_count(risk.module, name = "n2") %>% filter(n2>3) %>%
  complete(risk.module,nesting(Description,ID,enrich_type))  %>%
  mutate(risk.module         = gsub("\n"," ",risk.module),
         Tissue              = strsplit2(risk.module,">")[,1],
         Modules             = strsplit2(risk.module,">")[,2],
         Modules_highlighted = ifelse(risk.module %in% priortised_modules, paste0("underline(italic(",tolower(Modules),"))") , tolower(Modules)),
         Tissue              = factor(gsub("BRNCTX","Hartl2021_BRNCTX",gsub(" ","~",Tissue)), levels = gsub("BRNCTX","Hartl2021_BRNCTX",gsub(" ","~",tissue_order))),
         ID                  = factor(ID,levels = unlist(GO.top3.ordered), labels = Description[match(unlist(GO.top3.ordered),ID)])) %>% 
  ggplot(aes(x = ID, y = interaction(Modules_highlighted,Tissue,sep=">"), color = -log10(p.adjust), fill = -log10(p.adjust), size = Count)) + 
  geom_tile(width = 1, height = rel(200), color = prismatic::clr_lighten("black",shift = 0.75), size = 0.1, position = "identity") +
  facet_nested(Tissue + Modules_highlighted~enrich_type, scale = "free", space = "free", switch = "y", nest_line = element_line(colour = "darkgrey", linetype = 1), strip = strip_nested(clip = "off", size = "variable", background_y = elem_list_rect(fill = c("darkgrey", "grey80")), by_layer_y = T), labeller = label_parsed, drop = T) +
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 90), na.translate = FALSE, expand = expansion(mult = c(.01)))+
  scale_fill_gradientn(colours = c("black","yellow","darkblue"), limits = c(0,15), oob =scales::oob_squish, na.value = "white")+
  scale_color_gradientn(colours = c("black","yellow","darkblue"), limits = c(0,15), oob =scales::oob_squish, na.value = "white")+
  coord_cartesian(clip = "off")+
  xlab(label = NULL) + ylab(label = NULL) + #labs(fill = .y) +
  guides(fill  = guide_colorbar(nrow = 1, title = paste0("-log10(bonf)"), title.position = "top", title.theme = element_text(face = "bold", hjust = 0.2), label.position = "bottom", label.theme = element_text(angle = 0, hjust = 0.5), keywidth = unit(20,"pt"), keyheight = unit(20,"pt")),
         color = guide_colorbar(nrow = 1, title = paste0("-log10(bonf)"), title.position = "top", title.theme = element_text(face = "bold", hjust = 0.2), label.position = "bottom", label.theme = element_text(angle = 0, hjust = 0.5), keywidth = unit(20,"pt"), keyheight = unit(20,"pt")),
         size  = "none")+
  theme_light(base_size = 11) + 
  theme(
    strip.text               = element_text(color = "black"), 
    axis.title               = element_blank(),
    axis.text.x              = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.ticks.length.x      = unit(0,"pt"),
    axis.ticks.length.y.left = unit(0,"pt"),
    axis.text.y              = element_blank(),
    strip.text.y.left        = element_text(angle = 0, margin = margin(5,5,5,5,"pt"), vjust = 0.5, size = 14),
    strip.text.x             = element_text(angle = 0, margin = margin(5,5,5,5,"pt"), vjust = 0.5, size = 14),
    strip.placement          = "outside",
    legend.position          = "bottom",
    legend.direction         = "horizontal",
    legend.spacing.x         = unit(1.0, 'pt'),
    plot.title               = element_text(hjust = 0.5),
    panel.background         = element_blank(),
    panel.border             = element_rect(size = 0.1, color = prismatic::clr_lighten("black",shift = 0.75)),
    #panel.border             = element_blank(),
    plot.margin              = margin(3,3,3,3), 
    panel.grid.major.y       = element_line(color = "grey", size = 0.0),
    panel.grid.major.x       = element_line(color = "grey", size = 0.0)
  )

# setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Figures")
# ggsave(GO.plot, filename = "GO.plots.new.svg"            ,width = 17, height = 11)


#########Final TF plots##########
library(ggh4x)

#Format TF enrichment results into df
temp = new.env()
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Enrichments\\enrich_TF.gprofiler_META_Perinatal.RData",envir = temp)
ne$TF.df = rbind(ne$TF.df,map_dfr(.id = "modules",temp$result,~.x) %>% mutate(tissue_age = "Replication Perinatal",) %>% group_nest(tissue_age,modules))
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Enrichments\\enrich_TF.gprofiler_META_Juvenile.RData",envir = temp)
ne$TF.df = rbind(ne$TF.df,map_dfr(.id = "modules",temp$result,~.x) %>% mutate(tissue_age = "Replication Juvenile",) %>% group_nest(tissue_age,modules))
names(ne$TF.df$data) = paste0(ne$TF.df$tissue_age,">",ne$TF.df$modules)

TF = ne$TF.df
names(TF$data) = map(names(TF$data), ~{ 
  .x = gsub("Gandal2018>"     ,"Gandal2018a>"   ,.x)
  .x = gsub("Gandal2018PE>"   ,"Gandal2018b>"   ,.x)
  .x = gsub("Gandal2018PE_cs>","Gandal2018b_cs>",.x)
  .x = gsub(">WHOLE_BRAIN"    ,">BW"            ,.x)
  .x
})
names(TF$data) = gsub("\n"," ",names(TF$data))

ss = strsplit2(mega.risk.modules$risk.modules.noGrey,">")
#Keep only significant TF results
TF.df = TF$data[unique(c(paste0(ss[,1],">",ss[,2]),paste0(ss[,1],">",tolower(ss[,2]))))] %>% .[lengths(.)>0] %>% bind_rows(.id = "risk.module") %>% filter(p.adjust <0.05) %>% mutate(Term_ID = strsplit2(Term_ID,"_")[,1], Description = as.character(Description))
#Further filter to keep only top-3 significant processes
TF.top3 = TF.df %>% group_by(risk.module) %>% slice_min(p.adjust,n=3) %>% pull(Description) %>% unique
TF.top3.ordered = TF.df %>% filter(Description %in% unique(TF.top3)) %>% group_by(Description) %>% dplyr::count() %>% arrange(desc(n)) %>% pull(Description)

  
TF.plot = TF.df %>%
  filter(Description %in% TF.top3.ordered) %>% add_count(Description,name = "n1") %>% filter(n1>4) %>% add_count(risk.module, name = "n2") %>% filter(n2>3) %>%
  complete(risk.module,nesting(Description)) %>%
  mutate(risk.module         = gsub("\n"," ",risk.module),
         Tissue              = strsplit2(risk.module,">")[,1],
         Modules             = strsplit2(risk.module,">")[,2],
         Modules_highlighted = ifelse(risk.module %in% (priortised_modules %>% strsplit2(.,">") %>% c(paste0(.[,1],">",.[,2]),paste0(.[,1],">",tolower(.[,2]))) %>% unique), paste0("underline(italic(",tolower(Modules),"))") , tolower(Modules)),
         Tissue              = factor(gsub(" ","~",Tissue), levels = gsub(" ","~",tissue_order)),
         ID                  = factor(Description,levels = unique(TF.top3.ordered))) %>%
  ggplot(aes(x = ID, y = interaction(Modules_highlighted,Tissue,sep=">"), color = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_tile(width = 1, height = rel(200), color = prismatic::clr_lighten("black",shift = 0.75), size = 0.1, position = "identity") +
  facet_nested(Tissue + Modules_highlighted~., scale = "free", space = "free", switch = "y", nest_line = element_line(colour = "darkgrey", linetype = 1), strip = strip_nested(clip = "off", size = "variable", background_y = elem_list_rect(fill = c("darkgrey", "grey80")), by_layer_y = T), labeller = label_parsed, drop = T) +
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 40), na.translate = FALSE, expand = expansion(mult = c(-0.01)))+
  scale_fill_gradientn(colours = c("black","yellow","darkblue"), limits = c(0,30), oob =scales::oob_squish, na.value = "white")+
  scale_color_gradientn(colours = c("black","yellow","darkblue"), limits = c(0,30), oob =scales::oob_squish, na.value = "white")+
  coord_cartesian(clip = "off")+
  xlab(label = NULL) + ylab(label = NULL) + #labs(fill = .y) +
  guides(fill  = guide_colorbar(nrow = 1, title = paste0("-log10(bonf)"), title.position = "top", title.theme = element_text(face = "bold", hjust = 0.2), label.position = "bottom", label.theme = element_text(angle = 0, hjust = 0.5), keywidth = unit(20,"pt"), keyheight = unit(20,"pt")),
         color = guide_colorbar(nrow = 1, title = paste0("-log10(bonf)"), title.position = "top", title.theme = element_text(face = "bold", hjust = 0.2), label.position = "bottom", label.theme = element_text(angle = 0, hjust = 0.5), keywidth = unit(20,"pt"), keyheight = unit(20,"pt")),
         size  = "none")+
  theme_light(base_size = 11) +
  theme(
    strip.text               = element_text(color = "black"),
    axis.title               = element_blank(),
    axis.text.x              = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13),
    axis.ticks.length.x      = unit(0,"pt"),
    axis.ticks.length.y.left = unit(0,"pt"),
    axis.text.y              = element_blank(),
    strip.text.y.left        = element_text(angle = 0, margin = margin(5,5,5,5,"pt"), vjust = 0.5, size = 14),
    strip.text.x             = element_text(angle = 0, margin = margin(5,5,5,5,"pt"), vjust = 0.5, size = 14),
    strip.placement          = "outside",
    legend.position          = "bottom",
    legend.direction         = "horizontal",
    legend.spacing.x         = unit(1.0, 'pt'),
    plot.title               = element_text(hjust = 0.5),
    panel.background         = element_blank(),
    panel.border             = element_rect(size = 0.1, color = prismatic::clr_lighten("black",shift = 0.75)),
    plot.margin              = margin(3,3,3,3),
    panel.grid.major.y       = element_line(color = "grey", size = 0.0),
    panel.grid.major.x       = element_line(color = "grey", size = 0.0)
  )

setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Figures")
#ggsave(TF.plot, filename = "TF.plots.svg"            ,width = 15, height = 9)

#############################################################
################### Consensus genes plot#####################
#############################################################

#Read the consensus list
cc <- readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/SNPextension/consensus_genes/consensus_gene_list_new.withHartl.rds")

#Read consensus genes list full annotation with genename
cc_full <- readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/SNPextension/consensus_genes/consensus_gene_list_full_annotation.withHartl.rds")
cc = map(cc,~{names(.x) = cc_full$external_gene_name[match(.x,cc_full$ensembl_gene_id)];.x})

#Read the regular gene list
gm <- readRDS("C:/Users/mpariha1/Desktop/upload_to_google_drive/OneDrive - Johns Hopkins/Shared Data/Giulio's ML code/Enrichments/grch38[PGC125new]/gene-module list (wide_form_test) (all networks)[grch38].rds")
names(gm) = new.names[names(gm)]
names(gm) = gsub("Gandal2018$"     , "Gandal2018a"   , names(gm))
names(gm) = gsub("Gandal2018PE$"   , "Gandal2018b"   , names(gm))
names(gm) = gsub("Gandal2018PE_cs$", "Gandal2018b_cs", names(gm))

#Read Hartl gene-module list
hm <- readRDS("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Hartl2021\\Hartl2021-gene-module-list.rds")
hm = map(hm,~{
  names(.x) = gsub("WHOLE_BRAIN","BW"   ,names(.x))
  names(.x) = gsub("CEREBELLUM" ,"CEREB",names(.x))
  .x
})

consensus.postnatal.tissue.order = c("DLPFC Juvenile","DLPFC Adult","DLPFC Older Adult",
                                      "Fromer2016_case","Fromer2016_control","Gandal2018a","BRNCTX",
                                      "Pergola2017","Pergola2019","Pergola2020","Radulescu2020")

consensus.prenatal.tissue.order = c("DLPFC Perinatal","Li2018","Walker2019","Werling2020") 
 
#Simple 0/1 heatmap (blue/grey)
gm1 = c(gm,list(BRNCTX = hm$BRNCTX)) %>% tibble::enframe(name = "network",value = "list1") %>%
  unnest_longer(list1,values_to = "list2", indices_to = "modules") %>%
  unnest_longer(list2,values_to = "ensembl") %>%
  unite(network,modules,col = "ID", sep = ">",remove = F) %>%
  complete(ensembl,network) %>%
  mutate(gene = cc_full$external_gene_name[match(ensembl,cc_full$ensembl_gene_id)],
         isExpressed          = ifelse(is.na(ID),0,1)                                                                ,
         inPostnatalConsensus = ifelse(ensembl %in% cc$SCZ_postnatal_positive.new,1,0)                               ,
         inPrenatalConsensus  = ifelse(ensembl %in% cc$SCZ_prenatal_positive,1,0)                                    ,
         inConsensus          = ifelse(ensembl %in% c(cc$SCZ_prenatal_positive,cc$SCZ_postnatal_positive.new),1,0)   ,
         isRiskModules        = ifelse(tolower(ID) %in% tolower(unlist(mega.risk.modules$consensus.risk.modules)),1,0)                      ,
         final                = case_when(
           inConsensus ==1 & isRiskModules ==1 & isExpressed==1 ~ "skyblue",
           inConsensus ==1 & isRiskModules ==0 & isExpressed==1 ~ "#F083FAFF",
           inConsensus ==1 & isRiskModules ==0 & isExpressed==0 ~ "#474747FF"     ,
           TRUE                                                 ~ "white"    
         ),
         networkType = case_when(
           network %in% consensus.postnatal.tissue.order ~ "Postnatal\nNetwork",
           network %in% consensus.prenatal.tissue.order  ~ "Perinatal\nNetwork",
           TRUE                                          ~ "Other"
           ),
         geneType = case_when(
           ensembl %in% cc$SCZ_prenatal_positive      ~ "Perinatal\nConsensus\nGenes",
           ensembl %in% cc$SCZ_postnatal_positive.new ~ "Postnatal\nConsensus\nGenes",
           TRUE                                       ~ "Other Genes"
           ))%>% 
  filter(inConsensus==1 & network %in% c(consensus.postnatal.tissue.order,consensus.prenatal.tissue.order))


gg.consensus= gm1 %>% 
  mutate(network = factor(gsub("BRNCTX","Hartl2021_BRNCTX",network), levels = c(gsub("BRNCTX","Hartl2021_BRNCTX",consensus.postnatal.tissue.order),consensus.prenatal.tissue.order)),
         gene    = factor(gene   , levels = names(sort(table(gm1$gene))))) %>% 
  ggplot(aes(x = gene, y = network)) +
  geom_tile(aes(fill = final, height = ifelse(isRiskModules,0.85,0.85),width = ifelse(isRiskModules,0.85,0.85)), color = "grey")+
  facet_nested(networkType+network~geneType+gene, scale = "free", space = "free", switch = "y", 
               nest_line = element_line(colour = "darkgrey", linetype = 1), 
               strip = strip_nested(clip = "off", size = "variable", 
                                    background_y = elem_list_rect(fill = c("darkgrey", "grey80")),
                                    text_y = elem_list_text(angle = c(90,0)), by_layer_y = T,
                                    background_x = elem_list_rect(fill = c("darkgrey", "grey80")),
                                    text_x = elem_list_text(angle = c(0,90)), by_layer_x = T),drop = T) +
  scale_x_discrete(position = "top",label = function(x) stringr::str_trunc(x, 40), na.translate = FALSE, expand = expansion(mult = c(.01)))+
  scale_fill_identity(guide = "legend",labels = c("In a SCZ module","Not in a SCZ module", "Not in the network"), breaks = c("skyblue","#F083FAFF","#474747FF"))+#scale_fill_manual(values = c("lightblue","grey"))+
  scale_size_manual(values = c(0.75,0.45),breaks = c("a1","a2"))+
  coord_cartesian(clip = "off")+
  xlab(label = NULL) + ylab(label = NULL) + #labs(fill = .y) +
  guides(size  = "legend")+
  theme_light(base_size = 11) + 
  theme(
    strip.text               = element_text(color = "black"), 
    axis.title               = element_blank(),
    axis.text.x.top          = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0, size = 15),
    axis.ticks.length.x      = unit(0,"pt"),
    axis.ticks.length.y.left = unit(0,"pt"),
    axis.text.y              = element_blank(),
    strip.text.y.left        = element_text(angle = 0, margin = margin(2,4,2,4,"pt"), vjust = 0.5, size = 14),
    strip.text.x             = element_text(angle = 0, margin = margin(2,4,2,4,"pt"), vjust = 0.5, size = 14),
    strip.placement          = "outside",
    legend.position          = "bottom",
    legend.direction         = "horizontal",
    legend.spacing.x         = unit(1.0, 'pt'),
    plot.title               = element_text(hjust = 0.5),
    panel.background         = element_blank(),
    panel.border             = element_rect(size = 0.0, color = prismatic::clr_lighten("black",shift = 0.75)),
    #panel.border             = element_blank(),
    plot.margin              = margin(6,6,6,6), 
    panel.grid.major.y       = element_blank(),#element_line(color = prismatic::clr_darken("grey",shift = 0.05)),
    panel.grid.major.x       = element_blank(),#element_line(color = prismatic::clr_darken("grey",shift = 0.05)),
    panel.spacing            = unit(3,"pt") 
  )

setwd("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Paper new analysis\\Revision\\metanets\\Figures")
#ggsave(gg.consensus, filename = "consensus_postnatal_heatmap_new3.svg"            ,width = 11, height = 6.5)
