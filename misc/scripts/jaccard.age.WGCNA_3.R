##### multiplot ####
library(ggplot2)
library(grid)
library(gridExtra)
require(ggrepel)

data.null = rbind(Caudate.null,Dlpfc.null,Hippo.null)
data.null$age_window = factor(data.null$age_window,levels = c("perinatal-juvenile","juvenile-adult","adult-older adult"))
data.null$type = as.character(data.null$type)
data.null$type[data.null$type %in% "SCZ"] = "SCZ risk genes"; data.null$type[data.null$type %in% "ALL"] = "ALL genes"
data.null$type = factor(data.null$type,levels = c("SCZ risk genes","ALL genes"))
data$type = as.character(data$type)
data$type[data$type %in% "SCZ"] = "SCZ risk genes"; data$type[data$type %in% "ALL"] = "ALL genes"
data$type = factor(data$type,levels = c("SCZ risk genes","ALL genes"))
data$tissue = gsub("caudate","CN",
                   gsub("dlpfc","DLPFC",
                        gsub("hippocampus","HP",data$tissue)))
data.null$tissue = gsub("caudate","CN",
                        gsub("dlpfc","DLPFC",
                             gsub("hippocampus","HP",data.null$tissue)))


svg("Caudate.dlpfc.hippo_jaccard.index_boxplot2_different.universe.svg", width=6.5, height=13)
# ggplot(data, aes(fill = type, y = JI, x = age_window)) + 
#     geom_bar(position="dodge",stat = "identity") + facet_grid(. ~ tissue) + xlab("age stages") + ylab("Jaccard Index") +
#         scale_fill_discrete(name = "type of genes") + 
#         theme(axis.title = element_text(size = 25,face = "bold"),
#               axis.text.x = element_text(size = 16),
#               axis.text.y = element_text(size = 19),
#               strip.text = element_text(size = 22),
#               legend.text = element_text(size = 18),
#               legend.title = element_text(size = 22,face = "bold"))

ggplot(data.null, aes(y = JI, x = age_window)) + 
  geom_boxplot(outlier.shape = 21,width = 0.4,aes(fill = age_window),data = data.null) + 
  geom_point(data = data,color = "black",size = 3, shape = 23,fill = "red") +
  facet_grid(tissue ~ type,scales = "free_x") + xlab("") + ylab("Jaccard Index") +
  scale_fill_discrete(name = "Age transitions") + 
  theme_bw() +
  theme(axis.title = element_text(size = 25,face = "bold"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(angle = 90,size = 19,vjust = 0.5),
        strip.text = element_text(size = 19),
        legend.text = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        
        # legend.title = element_text(size = 16,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(2.6, 3, 2.6, 3), "lines"))

dev.off()


save(data,data.null, file = "grey.age.transitions.jaccard.RData")
