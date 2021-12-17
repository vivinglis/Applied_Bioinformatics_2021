# Intersection analyses 
# Spearman

rm(list = ls())

# Packages
library(ggplot2)
library(cowplot)
library(ggVennDiagram)

# Read in tables
corr_avgs <- read.table("Applied_Bioinformatics_2021/data_tables/1_spearman_corr_avgs.tsv", header = TRUE) 
corr_avgs_pro <- read.table("Applied_Bioinformatics_2021/data_tables/1_spearman_corr_avgs_with_pro.tsv", header = TRUE) 
corr_rep <- read.table("Applied_Bioinformatics_2021/data_tables/1_spearman_corr_rep.tsv", header = TRUE) 
corr_rep_pro <- read.table("Applied_Bioinformatics_2021/data_tables/1_spearman_corr_rep_with_pro.tsv", header = TRUE)  

top_num <- 100

# Avgs without promotor
veg <- subset(corr_avgs, Tissue == "veg")
sex <- subset(corr_avgs, Tissue == "sex")

veg_sa <- subset(veg, Allel == "sa")
veg_ba <- subset(veg, Allel == "BA")
sex_sa <- subset(sex, Allel == "sa")
sex_ba <- subset(sex, Allel == "BA")

top_veg_sa <- head(veg_sa[order(-veg_sa$Value, decreasing =TRUE),], top_num)
top_veg_ba <- head(veg_ba[order(-veg_ba$Value, decreasing =TRUE),], top_num)
top_sex_ba <- head(sex_ba[order(-sex_ba$Value, decreasing =TRUE),], top_num)
top_sex_sa <- head(sex_sa[order(-sex_sa$Value, decreasing =TRUE),], top_num)

top_genes <- data.frame(cbind(top_veg_sa$Geneid, top_veg_ba$Geneid, top_sex_sa$Geneid, top_sex_ba$Geneid))
names(top_genes) <- c("veg_sa", "veg_ba", "sex_sa", "sex_ba")

a <- ggVennDiagram(top_genes) + labs(title = "Avgs without promotor")

# Avgs with promotor, 
veg <- subset(corr_avgs_pro, Tissue == "veg")
sex <- subset(corr_avgs_pro, Tissue == "sex")

veg_sa <- subset(veg, Allel == "sa")
veg_ba <- subset(veg, Allel == "BA")
sex_sa <- subset(sex, Allel == "sa")
sex_ba <- subset(sex, Allel == "BA")

top_veg_sa <- head(veg_sa[order(-veg_sa$Value, decreasing =TRUE),], top_num)
top_veg_ba <- head(veg_ba[order(-veg_ba$Value, decreasing =TRUE),], top_num)
top_sex_ba <- head(sex_ba[order(-sex_ba$Value, decreasing =TRUE),], top_num)
top_sex_sa <- head(sex_sa[order(-sex_sa$Value, decreasing =TRUE),], top_num)

top_genes <- data.frame(cbind(top_veg_sa$Geneid, top_veg_ba$Geneid, top_sex_sa$Geneid, top_sex_ba$Geneid))
names(top_genes) <- c("veg_sa", "veg_ba", "sex_sa", "sex_ba")

b <- ggVennDiagram(top_genes) + labs(title = "Avgs with promotor")

#Rep without promotor
veg <- subset(corr_rep, Tissue == "veg")
sex <- subset(corr_rep, Tissue == "sex")

veg_sa <- subset(veg, Allel == "sa")
veg_ba <- subset(veg, Allel == "BA")
sex_sa <- subset(sex, Allel == "sa")
sex_ba <- subset(sex, Allel == "BA")

top_veg_sa <- head(veg_sa[order(-veg_sa$Value, decreasing =TRUE),], top_num)
top_veg_ba <- head(veg_ba[order(-veg_ba$Value, decreasing =TRUE),], top_num)
top_sex_ba <- head(sex_ba[order(-sex_ba$Value, decreasing =TRUE),], top_num)
top_sex_sa <- head(sex_sa[order(-sex_sa$Value, decreasing =TRUE),], top_num)

top_genes <- data.frame(cbind(top_veg_sa$Geneid, top_veg_ba$Geneid, top_sex_sa$Geneid, top_sex_ba$Geneid))
names(top_genes) <- c("veg_sa", "veg_ba", "sex_sa", "sex_ba")

c <- ggVennDiagram(top_genes) + labs(title = "Rep without promotor")


# Rep with promotor, 
veg <- subset(corr_rep_pro, Tissue == "veg")
sex <- subset(corr_rep_pro, Tissue == "sex")

veg_sa <- subset(veg, Allel == "sa")
veg_ba <- subset(veg, Allel == "BA")
sex_sa <- subset(sex, Allel == "sa")
sex_ba <- subset(sex, Allel == "BA")

top_veg_sa <- head(veg_sa[order(-veg_sa$Value, decreasing =TRUE),], top_num)
top_veg_ba <- head(veg_ba[order(-veg_ba$Value, decreasing =TRUE),], top_num)
top_sex_ba <- head(sex_ba[order(-sex_ba$Value, decreasing =TRUE),], top_num)
top_sex_sa <- head(sex_sa[order(-sex_sa$Value, decreasing =TRUE),], top_num)

top_genes <- data.frame(cbind(top_veg_sa$Geneid, top_veg_ba$Geneid, top_sex_sa$Geneid, top_sex_ba$Geneid))
names(top_genes) <- c("veg_sa", "veg_ba", "sex_sa", "sex_ba")

d <- ggVennDiagram(top_genes) + labs(title = "Rep with promotor")

# Plot all in same plot
plot_row <- plot_grid(a, b, c, d, labels = c("A", "B", "C", "D"))

title <- ggdraw() + 
  draw_label(
    "Comparison of top 100 genes between samples, Spearman correlation values",
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
