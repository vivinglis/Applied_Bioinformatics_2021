# Pearson correlation values plots

rm(list = ls())

# Packages
library(ggplot2)
library(cowplot)

# Read in tables
corr_avgs <- read.table("Applied_Bioinformatics_2021/data_tables/1_pearson_corr_avgs.tsv", header = TRUE) 
corr_avgs_pro <- read.table("Applied_Bioinformatics_2021/data_tables/1_pearson_corr_avgs_with_pro.tsv", header = TRUE) 
corr_rep <- read.table("Applied_Bioinformatics_2021/data_tables/1_2_pearson_corr_rep.tsv", header = TRUE) 
corr_rep_pro <- read.table("Applied_Bioinformatics_2021/data_tables/1_2_pearson_corr_rep_with_pro.tsv", header = TRUE) 
 
# All genes
corr_avgs$Tissue <- as.factor(corr_avgs$Tissue)
avgs <- ggplot(corr_avgs, aes(x=Tissue, y=Value, fill=Tissue))+
  geom_violin() + 
  labs(title="Using average of replicas", y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")

corr_avgs_pro$Tissue <- as.factor(corr_avgs_pro$Tissue)
avgs_p <- ggplot(corr_avgs_pro, aes(x=Tissue, y=Value, fill=Tissue))+
          geom_violin() + 
          labs(title="Using average of replicas with promotor", y = "Correlation value")+
          geom_boxplot(width=0.1, fill = "white")

corr_rep$Tissue <- as.factor(corr_rep$Tissue)
rep <- ggplot(corr_rep, aes(x=Tissue, y=Value, fill=Tissue))+
       geom_violin() + 
       labs(title="Using replicas", y = "Correlation value")+
       geom_boxplot(width=0.1, fill = "white")

corr_rep_pro$Tissue <- as.factor(corr_rep_pro$Tissue)
rep_p <- ggplot(corr_rep, aes(x=Tissue, y=Value, fill=Tissue))+
         geom_violin() + 
         geom_boxplot(width=0.1, fill = "white") +
         labs(title="Using replicas with promotor", y = "Correlation value")

plot_row <- plot_grid(avgs, avgs_p, rep, rep_p, labels = c("A", "B", "C", "D"))

title <- ggdraw() + 
  draw_label(
    "Pearson correation values for all genes",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(
  title, plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

#---------------------------------
# Top correlation values
top_num <- 100

veg <- subset(corr_avgs, Tissue == "veg")
sex <- subset(corr_avgs, Tissue == "sex")
top_veg <- head(veg[order(-veg$Value, decreasing =TRUE),], top_num)
top_sex <- head(sex[order(-sex$Value, decreasing =TRUE),], top_num)
top <- rbind(top_veg, top_sex)
avgs <- ggplot(top, aes(x=Tissue, y=Value, fill=Tissue))+
        geom_violin() + 
        labs(title="Using average of replicas ", y = "Correlation value")+
        geom_boxplot(width=0.1, fill = "white") +
        ylim(c(-1,-0.88))

veg <- subset(corr_avgs_pro, Tissue == "veg")
sex <- subset(corr_avgs_pro, Tissue == "sex")
top_veg <- head(veg[order(-veg$Value, decreasing =TRUE),], top_num)
top_sex <- head(sex[order(-sex$Value, decreasing =TRUE),], top_num)
top <- rbind(top_veg, top_sex)
avgs_p <- ggplot(top, aes(x=Tissue, y=Value, fill=Tissue))+
          geom_violin() + 
          labs(title="Using average of replicas  with promotor ", y = "Correlation value")+
          geom_boxplot(width=0.1, fill = "white")+
          ylim(c(-1,-0.88))

veg <- subset(corr_rep, Tissue == "veg")
sex <- subset(corr_rep, Tissue == "sex")
top_veg <- head(veg[order(-veg$Value, decreasing =TRUE),], top_num)
top_sex <- head(sex[order(-sex$Value, decreasing =TRUE),], top_num)
top <- rbind(top_veg, top_sex)
rep <- ggplot(top, aes(x=Tissue, y=Value, fill=Tissue))+
       geom_violin() + 
       labs(title="Using replicas", y = "Correlation value")+
       geom_boxplot(width=0.1, fill = "white")+
       ylim(c(-1,-0.88))

veg <- subset(corr_rep_pro, Tissue == "veg")
sex <- subset(corr_rep_pro, Tissue == "sex")
top_veg <- head(veg[order(-veg$Value, decreasing =TRUE),], top_num)
top_sex <- head(sex[order(-sex$Value, decreasing =TRUE),], top_num)
top <- rbind(top_veg, top_sex)
rep_p <- ggplot(top, aes(x=Tissue, y=Value, fill=Tissue))+
         geom_violin() + 
         labs(title="Using replicas  with promotor ", y = "Correlation value")+
         geom_boxplot(width=0.1, fill = "white")+
         ylim(c(-1,-0.88))

plot_row <- plot_grid(avgs, avgs_p, rep, rep_p, labels = c("A", "B", "C", "D"))

title <- ggdraw() + 
  draw_label(
    "Pearson correation values for top 100 genes in each tissue",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(
  title, plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)
