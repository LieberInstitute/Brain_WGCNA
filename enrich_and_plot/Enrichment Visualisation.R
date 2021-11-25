####################################################################################################
#             Script to plot main enrichment plots for SCZ risk modules for all networks           #
#        Script to plot GO, TF and Pathology enrichments for SCZ risk modules for all networks     #
#           "Regional-coexpression", "Age-period" and "Cell-population enrichment" studies         #
#                                                                                                  #
####################################################################################################

library(limma)
library(gtools)

library(stringr)
library(ggplot2)

library(WGCNA)
library(dplyr)
library(purrr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(magrittr)
library(pheatmap)
library(metap)

##Read updated network names
new.names = read.csv("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\old-to-new Network names all.csv") %>% tibble::deframe()

##Read list of SCZ risk modules:
SCZ_risk_modules_df = readRDS(file =  "C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\SCZ.risk.modules.df.rds") 

##Read enrichment results prepared for figures
load("C:\\Users\\mpariha1\\Desktop\\upload_to_google_drive\\OneDrive - Johns Hopkins\\Code for WGCNA paper\\all.figure.df.RData",figure.df <-new.env())

## Sort each enrichment df to match order of SCZ_risk_modules_df
figure.df.li = imap(as.list(figure.df), ~ {
  if ("new.network.module" %in% colnames(.x)) {
    .x[match(.x$new.network.module,SCZ_risk_modules_df$new.network.module),] %>% tibble::column_to_rownames("new.network.module")
  } else{
    .x[rownames(.x) %in% SCZ_risk_modules_df$new.network.module,]
  }
})

# tissue.levels_a = c(
#   "CN","HP","DLPFC",
#   "Fromer2016_case","Fromer2016_control",
#   "Gandal2018",
#   "Li2018",
#   "Pergola2017","Pergola2019","Pergola2020",
#   "Radulescu2020",
#   "Walker2019",
#   "Werling2020"
#   )
# 
# tissue.levels_b = c(
#   "DG.noQSVA",
#   "HP.QSVA", "HP.noQSVA"
# )
# 
# tissue.levels_c = c(
#   "CN Juvneile"    ,"CN Adult"      , "CN Older Adult",
#   "HP Perinatal"   ,"HP Juvneile"   ,"HP Adult"       , "HP Older Adult",
#   "DLPFC Perinatal","DLPFC Juvneile","DLPFC Adult"    , "DLPFC Older Adult"
# )


hl = HeatmapAnnotation(which = "row", 
                       Network  = anno_block(labels = unique(SCZ_risk_modules_df$new.network),  
                                             gp = gpar(fill = 0, border = "black", lwd = 0.3, lty = 1),
                                             labels_gp = gpar(col = "black",fontsize = 13, fontface = "bold", just = c("right","top")),
                                             labels_rot = 0, show_name = F, width = unit(45,"mm"), height = unit(25,"mm")
                       ),
                       Modules  = anno_text(x = SCZ_risk_modules_df$module, just = "center", location = 0.5,  #both just and Location to move labels up/down
                                            gp = gpar(fill = SCZ_risk_modules_df$module, col = "white", border = "black", fontsize = 12, lwd = 1.2, fontface = "bold"),
                                            height = max_text_height(SCZ_risk_modules_df$module)*3, width = max_text_width(SCZ_risk_modules_df$module)*1.5
                       ),
                       
                       gap = unit(2, "mm"),
                       show_legend = F
)


col_fun_count           = colorRamp2(c(0, 5, 9 ), c("black","violet", "blue"))                ## Col scale: a
col_fun_cont_very_small = colorRamp2(c(0 ,2, 5 ), c("black", "yellow", "red"))                ## Col scale: b
col_fun_cont_big        = colorRamp2(c(0, 10, 80, 100), c("black", "yellow", "red", "green")) ## Col scale: c


numbers_flag = TRUE
## Main heatmap (PGC HYPER.ALL, PGC HYPER.PC, PGC3 HYPER.NEG.ALL, PGC3 HYPER.NEG.PC, PGC3 PERM.ALL, PGC3 LOCI.WINDOWS)
if (numbers_flag) {my_cell_fun = function(j, i, x, y, w, h, col="white") {
  mat = cbind(figure.df.li$dat.block1,figure.df.li$dat.block2) %>% round(1)
  grid.text(mat[i, j], x, y)
  # if (mat[i,j] > 2) {
  #   grid.text("*", x, y, hjust = 0.5, vjust = 0.8, gp = gpar(fontsize=20, col = "green"))
  # }
}
}else {my_cell_fun = NA}

ht0 = Heatmap(cbind(figure.df.li$dat.block1,figure.df.li$dat.block2), na_col = "grey", name = "ht0", col = col_fun_count, cluster_rows = F, cluster_columns = F,
              row_title    = NULL ,row_title_gp    = gpar(fontsize = 0, fontface = "bold"), row_title_side    = "left", row_title_rot    = 90, row_labels = SCZ_risk_modules_df$new.network,
              column_title = NULL ,column_title_gp = gpar(fontsize = 0, fontface = "bold"), column_title_side = "top" , column_title_rot = 90,
              
              row_names_side    = "left", row_names_rot    = 0 , row_names_gp    = gpar(fontsize = 0, fontface = "bold"     ),
              column_names_side = "bottom" , column_names_rot = 90, column_names_gp = gpar(fontsize = 10, fontface = 2, lty = 1),
              rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
              show_heatmap_legend = F,
              cell_fun = my_cell_fun,
)



#HMAGMA and MAGMA
if (numbers_flag) {my_cell_fun = function(j, i, x, y, w, h, col="white") {
  mat = cbind(figure.df.li$HMAGMA.df,figure.df.li$MAGMA.df) %>% round(1)
  grid.text(mat[i, j], x, y)}
}else {my_cell_fun = NA}

ht1 = Heatmap(cbind(figure.df.li$HMAGMA.df,figure.df.li$MAGMA.df), name = "ht1", col = col_fun_cont_very_small, cluster_rows = F, cluster_columns = F,
              column_names_rot = 90, column_names_side = "bottom" , column_names_gp = gpar(fontsize = 10, fontface = 2, lty = 1),
              row_names_gp = gpar(fontsize = 0), row_labels = SCZ_risk_modules_df$new.network,
              rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
              
              show_heatmap_legend = F,
              cell_fun = my_cell_fun,
)


#TWAS SCZ, DEGS SCZ, DMGs SCZ, DEGS HUMAN.APES, LOF, DRUG
if (numbers_flag) {my_cell_fun = function(j, i, x, y, w, h, col="white") {
  mat = cbind(figure.df.li$TWAS.SCZ.df, figure.df.li$DEGs.SCZ.df, figure.df.li$DMGs.SCZ.df, figure.df.li$DEGs.HUMAN.APES.df, figure.df.li$LOF.df, figure.df.li$DRUG.df) %>% round(1)
  grid.text(mat[i, j], x, y)}
}else {my_cell_fun = NA}

ht2 = Heatmap(cbind(figure.df.li$TWAS.SCZ.df, figure.df.li$DEGs.SCZ.df, figure.df.li$DMGs.SCZ.df, figure.df.li$DEGs.HUMAN.APES.df, figure.df.li$LOF.df, figure.df.li$DRUG.df), name = "ht2", col = col_fun_cont_big, cluster_rows = F, cluster_columns = F,
              column_names_rot = 90, column_names_side = "bottom" , column_names_gp = gpar(fontsize = 10, fontface = 2, lty = 1),
              row_names_gp = gpar(fontsize = 0), row_labels = SCZ_risk_modules_df$new.network,
              rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
              
              show_heatmap_legend = F,
              cell_fun = my_cell_fun,
)



#Cell Specificity
if (numbers_flag) {my_cell_fun = function(j, i, x, y, w, h, col="white") {
  mat = figure.df.li$cell.specificity.df %>% round(1)
  grid.text(mat[i, j], x, y)}
}else {my_cell_fun = NA}

ht3 = Heatmap(figure.df.li$cell.specificity.df, name = "ht3", col = col_fun_cont_big,
              column_names_rot = 90,column_names_side = "bottom" , column_names_gp = gpar(fontsize = 10, fontface = 2, lty = 1), cluster_rows = F, cluster_columns = F,
              row_names_gp = gpar(fontsize = 0), row_labels = SCZ_risk_modules_df$new.network,
              rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
              
              show_heatmap_legend = F,
              cell_fun = my_cell_fun,
)


#Bin count
lgd_a = Legend(title = gt_render("  <sup><sup>a </sup></sup> #sig.bins   "), title_position = "topcenter",
              labels_gp = gpar(col = "red", font = 3), title_gp = gpar(fontsize = 12, fontface = "bold", just = 0, box_fill = "#FF000015", box_lty = 0),
              legend_height = unit(30, "mm"), legend_width = unit(40, "mm"), direction = "horizontal", border = "black",
              at = c(0,5,9), labels = c(0,5,9), col_fun = col_fun_count)

#PGC3 loci and HMAGMA
lgd_b = Legend(title = gt_render(" <sup><sup>b </sup></sup> -log10(bonf)  "), title_position = "topcenter",
               labels_gp = gpar(col = "red", font = 3), title_gp = gpar(fontsize = 12, fontface = "bold", just = 0, box_fill = "#0000FF15", box_lty = 0),
               legend_height = unit(30, "mm"), legend_width = unit(40, "mm"), direction = "horizontal", border = "black",
               col_fun = col_fun_cont_very_small)

#DEG,DMG,TWAS,LOF Drug etc (misc) and Cell specificity
lgd_c = Legend(title = gt_render(" <sup><sup>c </sup></sup> -log10(bonf)  "), title_position = "topcenter",
              labels_gp = gpar(col = "red", font = 3), title_gp = gpar(fontsize = 12, fontface = "bold", just = 0, box_fill = "#00FF0015", box_lty = 0),
              legend_height = unit(30, "mm"), legend_width = unit(40, "mm"), direction = "horizontal", border = "black",
              col_fun = col_fun_cont_big)

pd = packLegend(lgd_a, lgd_b, lgd_c, max_width = unit(20, "cm"), 
                direction = "horizontal", column_gap = unit(10, "mm"), row_gap = unit(1, "cm"))

# #Global options
ht = hl + ht0 + ht1 + ht2 + ht3

#setwd("C:/Users/mpariha1/Documents")
pdf("grch38.testComplexHeatmap(test).pdf", width = 15, height = 30, useDingbats = T, compress = F, pointsize = 5)

ht_opt(RESET = TRUE)


d = draw(ht, annotation_legend_list = pd,
         auto_adjust = T, merge_legends = T, main_heatmap = "ht0",
         ht_gap = unit(5, "mm"),
         padding = unit(c(6, 5, 6, 8), "mm"), #bottom, left, top, right paddings
         row_split = factor(SCZ_risk_modules_df$new.network, levels = unique(SCZ_risk_modules_df$new.network)), row_gap = unit(1.5,"lines"), cluster_row_slices = FALSE, cluster_column_slices = FALSE,
         heatmap_legend_side = "bottom"#,
) 

#Set background color for legend blocks
fnn = function(viewport,ww = unit(1.01, "npc"), hh = unit(1.01, "npc"), ff){
  seekViewport(viewport)
  grid.rect(width = ww, height = hh,  gp = gpar(fill = ff, lwd = 0, lty = 0))#, just = c("left","bottom"))
  seekViewport("global")
}

fnn("ht0_column_names_1", ff = "#FF000015")
fnn("ht1_column_names_1", ff = "#0000FF15")
fnn("ht2_column_names_1", ff = "#00FF0015")
fnn("ht3_column_names_1", ff = "#00FF0015")


ht_opt(RESET = TRUE)
dev.off()


rm(hl,ht0,ht1,ht2,ht3,lgd_a,lgd_b,lgd_c,d, pd, ht)

###################################################
###################################################

## Pathology enrichment plots ----
col_fun_count           = colorRamp2(c(0, 5, 9 ), c("black","violet", "blue"))        ## Col scale: patho
hl.patho = HeatmapAnnotation(which = "row", 
                             Network  = anno_block(labels = unique(SCZ_risk_modules_df$new.network),  
                                                   gp = gpar(fill = 0, border = "black", lwd = 0.3, lty = 1),
                                                   labels_gp = gpar(col = "black",fontsize = 13, fontface = "bold", just = c("right","top")),
                                                   labels_rot = 0, show_name = F, width = unit(45,"mm"), height = unit(25,"mm")
                             ),
                             Modules  = anno_text(x = SCZ_risk_modules_df$module, just = "center", location = 0.5,  #both just and Location to move labels up/down
                                                  gp = gpar(fill = SCZ_risk_modules_df$module, col = "white", border = "black", fontsize = 12, lwd = 1.2, fontface = "bold"),
                                                  height = max_text_height(SCZ_risk_modules_df$module)*3, width = max_text_width(SCZ_risk_modules_df$module)*1.5
                             ),
                             
                             gap = unit(2, "mm"),
                             show_legend = F
                             
)

figure.df.li$Patho.df[figure.df.li$Patho.df==0] = NA

ht.patho = Heatmap(figure.df.li$Patho.df, na_col = "grey", name = "ht.patho", col = col_fun_count, cluster_rows = F, cluster_columns = F, #cluster_columns = cluster_within_group(zz1, ann_row$tissue),
                   row_title    = NULL ,row_title_gp    = gpar(fontsize = 0, fontface = "bold"), row_title_side    = "left", row_title_rot    = 90,
                   column_title = "Disease pathology enrichments" ,column_title_gp = gpar(fontsize = 16, fontface = "bold"), column_title_side = "top" , column_title_rot = 0,
                   
                   row_names_side    = "left", row_names_rot    = 0 , row_names_gp    = gpar(fontsize = 10, fontface = "bold"     ), row_labels = SCZ_risk_modules_df$new.network,
                   column_names_side = "top" , column_names_rot = 90, column_names_gp = gpar(fontsize = 12, fontface = 2, lty = 1),
                   rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
                   
                   show_heatmap_legend = F,
                   cell_fun = NA,
)

lgd.patho = Legend(title = gt_render("#sig.bins"), title_position = "topcenter",
                   labels_gp = gpar(col = "red", font = 3), title_gp = gpar(fontsize = 12, fontface = "bold", just = 0, box_fill = "#FF000015", box_lty = 0),
                   legend_height = unit(30, "mm"), legend_width = unit(50, "mm"), direction = "horizontal", border = "black",
                   at = c(0,5,9), labels = c(0,5,9), col_fun = col_fun_count)

pd.patho = packLegend(lgd.patho, max_width = unit(25, "cm"), 
                      direction = "horizontal", column_gap = unit(10, "mm"), row_gap = unit(1, "cm"))

ht.patho = hl.patho + ht.patho

pdf("grch38.testPatho_ComplexHeatmap.pdf", width = 8, height = 16, useDingbats = T, compress = F, pointsize = 5)

d = draw(ht.patho, annotation_legend_list = pd.patho,
         auto_adjust = T, merge_legends = T,
         ht_gap = unit(c(1), "mm"),
         padding = unit(c(6, 5, 6, 8), "mm"), #bottom, left, top, right paddings)
         #legend_grouping = "original", 
         row_split = factor(SCZ_risk_modules_df$new.network, levels = unique(SCZ_risk_modules_df$new.network)), row_gap = unit(1.5,"lines"), cluster_row_slices = FALSE, cluster_column_slices = TRUE, cluster_rows = FALSE,
         heatmap_legend_side = "bottom",
         #width = unit(200,"mm"),
) 


dev.off()
rm(d, ht.patho, hl.patho, lgd.patho, pd.patho)

## GO enrichment plots ----
col_fun_cont_big        = colorRamp2(c(0, 15, 30), c("black", "yellow", "blue"))   #Enrichment

for (h in c("BP","MF","CC")){
  
  dat = figure.df.li$GO.df[,grepl(paste0(h,">"),colnames(figure.df.li$GO.df))]
  ann_col = strsplit2(colnames(dat),">") %>% as.data.frame %>% set_colnames(c("GO","ID","term"))
  ann_row = strsplit2(rownames(dat),">") %>% as.data.frame %>% set_colnames(c("new.network","module"))
  
  ht0 = Heatmap(as.matrix(dat), na_col = "grey", name = paste0("ht.",h), col = col_fun_cont_big, cluster_rows = F, cluster_columns = F,
                row_title    = NULL ,row_title_gp    = gpar(fontsize = 0, fontface = "bold"), row_title_side    = "left", row_title_rot    = 90,
                column_title = paste0("GO:",h," enrichments") ,column_title_gp = gpar(fontsize = 16, fontface = "bold"), column_title_side = "top" , column_title_rot = 0,
                
                row_names_side    = "left", row_names_rot    = 0 , row_names_gp    = gpar(fontsize = 0, fontface = "bold"     ),
                column_names_side = "bottom" , column_names_rot = 90, column_names_gp = gpar(fontsize = 9, fontface = 2, lty = 1), column_labels =strsplit2(str_wrap(ann_col$term, width = 50),"\n")[,1], 
                rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
                
                show_heatmap_legend = F,
                cell_fun = NA
  )
  
  hl = HeatmapAnnotation(which = "row", 
                         Network  = anno_block(labels = unique(ann_row$new.network), 
                                               gp = gpar(fill = 0, border = "black", lwd = 0.3, lty = 1),
                                               labels_gp = gpar(col = "black",fontsize = 13, fontface = "bold", just = c("right","top")),
                                               labels_rot = 0, show_name = F, width = unit(45,"mm"), height = unit(25,"mm")
                         ),
                         Modules  = anno_text(x = ann_row$module,  just = "center", location = 0.5,  #both just and Location to move labels up/down
                                              gp = gpar(fill = ann_row$module, col = "white", border = "black", fontsize = 14, lwd = 1.2, fontface = "bold"),
                                              height = max_text_height(ann_row$module)*3, width = max_text_width(ann_row$module)*1.5
                         ),
                         
                         gap = unit(2, "mm"),
                         show_legend = F
  )
  
  
  lgd1 = Legend(title = gt_render("-log10(bonf.val)"), title_position = "topcenter",
                labels_gp = gpar(col = "red", font = 3), title_gp = gpar(fontsize = 12, fontface = "bold", just = 0, box_fill = "#FF000015", box_lty = 0),
                legend_height = unit(30, "mm"), legend_width = unit(50, "mm"), direction = "horizontal", border = "black",
                at = c(0,15,30), labels = c(0,15,30), col_fun = col_fun_cont_big)
  
  pd = packLegend(lgd1, max_width = unit(25, "cm"), 
                  direction = "horizontal", column_gap = unit(10, "mm"), row_gap = unit(1, "cm"))
  
  
  
  pdf(paste0("grch38.testGO_ComplexHeatmap_",h,".pdf"), width = 12, height = 13, useDingbats = T, compress = F, pointsize = 5)
  
  ht = hl + ht0
  
  d = draw(ht, annotation_legend_list = pd,
           auto_adjust = T, merge_legends = T,
           ht_gap = unit(c(1), "mm"),
           padding = unit(c(6, 5, 6, 8), "mm"), #bottom, left, top, right paddings)
           #legend_grouping = "original", 
           row_split = factor(ann_row$new.network, levels = unique(ann_row$new.network)), row_gap = unit(1.5,"lines"), cluster_row_slices = TRUE, cluster_column_slices = TRUE, cluster_rows = FALSE,
           heatmap_legend_side = "bottom",
           #width = unit(200,"mm"),
  ) 
  dev.off()
}

rm(d, dat, h, ht, ht0, lgd1, pd, ann_col, ann_row)


## TF enrichment plots ----
col_fun_cont_big        = colorRamp2(c(0, 25, 50), c("black", "yellow", "blue"))   #Enrichments

dat = figure.df.li$TF.df
ann_row = strsplit2(rownames(dat),">") %>% as.data.frame %>% set_colnames(c("new.network","module"))


ht.TF = Heatmap(dat, na_col = "grey", name = "ht.patho", col = col_fun_cont_big, cluster_rows = F, cluster_columns = F, #cluster_columns = cluster_within_group(zz1, ann_row$tissue),
                   row_title    = NULL ,row_title_gp    = gpar(fontsize = 0, fontface = "bold"), row_title_side    = "left", row_title_rot    = 90,
                   column_title = "TF enrichments" ,column_title_gp = gpar(fontsize = 16, fontface = "bold"), column_title_side = "top" , column_title_rot = 0,
                   
                   row_names_side    = "left", row_names_rot    = 0 , row_names_gp    = gpar(fontsize = 0, fontface = "bold"     ),
                   column_names_side = "top" , column_names_rot = 90, column_names_gp = gpar(fontsize = 12, fontface = 2, lty = 1),
                   rect_gp = gpar(col = "black", lty = 1, lwd = 0.5),
                   
                   show_heatmap_legend = F,
                   cell_fun = NA,
)

hl = HeatmapAnnotation(which = "row", 
                       Network  = anno_block(labels = unique(ann_row$new.network),  
                                             gp = gpar(fill = 0, border = "black", lwd = 0.3, lty = 1),
                                             labels_gp = gpar(col = "black",fontsize = 13, fontface = "bold", just = c("right","top")),
                                             labels_rot = 0, show_name = F, width = unit(45,"mm"), height = unit(25,"mm")
                       ),
                       Modules  = anno_text(x = ann_row$module, just = "center", location = 0.5,  #both just and Location to move labels up/down
                                            gp = gpar(fill = ann_row$module, col = "white", border = "black", fontsize = 12, lwd = 1.2, fontface = "bold"),
                                            height = max_text_height(ann_row$module)*3, width = max_text_width(ann_row$module)*1.5
                       ),
                       
                       gap = unit(2, "mm"),
                       show_legend = F
)

lgd.TF = Legend(title = gt_render("-log10(bonf.val)"), title_position = "topcenter",
                labels_gp = gpar(col = "red", font = 3), title_gp = gpar(fontsize = 12, fontface = "bold", just = 0, box_fill = "#FF000015", box_lty = 0),
                legend_height = unit(30, "mm"), legend_width = unit(50, "mm"), direction = "horizontal", border = "black",
                at = c(0,25,50), labels = c(0,25,50), col_fun = col_fun_cont_big)

pd = packLegend(lgd.TF, max_width = unit(25, "cm"), 
                direction = "horizontal", column_gap = unit(10, "mm"), row_gap = unit(1, "cm"))



pdf("grch38.testTF_ComplexHeatmap(top-3).pdf", width = 12, height = 13, useDingbats = T, compress = F, pointsize = 5)

ht = hl + ht.TF
d = draw(ht, annotation_legend_list = pd,
         auto_adjust = T, merge_legends = T,
         ht_gap = unit(c(1), "mm"),
         padding = unit(c(6, 5, 6, 8), "mm"), #bottom, left, top, right paddings)
         #legend_grouping = "original", 
         row_split = factor(ann_row$new.network, levels = unique(ann_row$new.network)), row_gap = unit(1.5,"lines"), cluster_row_slices = TRUE, cluster_column_slices = TRUE, cluster_rows = FALSE,
         heatmap_legend_side = "bottom",
         #width = unit(200,"mm"),
) 

dev.off()
