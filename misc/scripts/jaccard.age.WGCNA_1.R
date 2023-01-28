library(ggplot2)

Caudate = data.frame(tissue = rep("caudate",4),
                     age_window = c(rep("juvenile-adult",2),rep("adult-older adult",2)),
                     type = rep(c("SCZ","ALL"),2),
                     JI = c(0.62,0.65,0.68,0.74))

DLPFC = data.frame(tissue = rep("dlpfc",6),
                     age_window = c(rep("perinatal-juvenile",2),rep("juvenile-adult",2),rep("adult-older adult",2)),
                     type = rep(c("SCZ","ALL"),3),
                     JI = c(0.44,0.49,0.68,0.72,0.70,0.75))

HIPPO = data.frame(tissue = rep("hippocampus",6),
                     age_window = c(rep("perinatal-juvenile",2),rep("juvenile-adult",2),rep("adult-older adult",2)),
                     type = rep(c("SCZ","ALL"),3),
                     JI = c(0.52,0.57,0.65,0.71,0.60,0.69))

data = rbind(Caudate,DLPFC,HIPPO)
data$age_window = factor(data$age_window,levels = c("perinatal-juvenile","juvenile-adult","adult-older adult"))
data$type = factor(data$type,levels = c("SCZ","ALL"))

jpeg("Caudate.dlpfc.hippo_jaccard.index_plot1.jpeg",units="in", width=18, height=8, res=250)
# ggplot(data, aes(fill = type, y = JI, x = age_window)) + 
#     geom_bar(position="dodge",stat = "identity") + facet_grid(. ~ tissue) + xlab("age stages") + ylab("Jaccard Index") +
#         scale_fill_discrete(name = "type of genes") + 
#         theme(axis.title = element_text(size = 25,face = "bold"),
#               axis.text.x = element_text(size = 16),
#               axis.text.y = element_text(size = 19),
#               strip.text = element_text(size = 22),
#               legend.text = element_text(size = 18),
#               legend.title = element_text(size = 22,face = "bold"))

ggplot(data, aes(fill = age_window, y = JI, x = type)) + 
    geom_bar(position="dodge",stat = "identity") + facet_grid(. ~ tissue) + xlab("") + ylab("Jaccard Index") +
        scale_fill_discrete(name = "Age periods") + 
        theme(axis.title = element_text(size = 25,face = "bold"),
              axis.text = element_text(size = 19),
              strip.text = element_text(size = 22),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 22,face = "bold"))

dev.off()
