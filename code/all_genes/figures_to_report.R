# Plots using no mating type separated correlation values

rm(list = ls())

# Packages
library(ggplot2)
library(cowplot)
library(ggVennDiagram)

corr_avgs_pro <- read.table("../../data_tables/all_genes/1_pearson_corr_avgs_with_pro_non_allel_sep.tsv", header = TRUE) 
corr_avgs <- read.table("../../data_tables/all_genes/1_pearson_corr_avgs_non_allel_sep.tsv", header = TRUE)
corr_rep <- read.table("../../data_tables/all_genes/1_pearson_corr_rep_non_allel_sep.tsv", header = TRUE) 
corr_rep_pro <- read.table("../../data_tables/all_genes/1_pearson_corr_rep_with_pro_non_allel_sep.tsv", header = TRUE)


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

s <- 15 

plot_row <- plot_grid(avgs + theme(text = element_text(size = s)),
                      avgs_p + theme(text = element_text(size = s)),
                      rep + theme(text = element_text(size = s)),
                      rep_p + theme(text = element_text(size = s)),
                      labels = c("A", "B", "C", "D"))

title <- ggdraw() + theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# Plot with avgs pro top 100 
top_num <- 100

veg <- subset(corr_avgs_pro, Tissue == "veg")
sex <- subset(corr_avgs_pro, Tissue == "sex")
top_veg <- head(veg[order(-veg$Value, decreasing =TRUE),], top_num)
top_sex <- head(sex[order(-sex$Value, decreasing =TRUE),], top_num)
top <- rbind(top_veg, top_sex)
avgs_p <- ggplot(top, aes(x=Tissue, y=Value, fill=Tissue))+
  geom_violin() + 
  labs(y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")+
  ylim(c(-1,-0.88))

avgs_p <- avgs_p + theme(text = element_text(size = s),
                         axis.text.x = element_text(size = 15), 
                         axis.text.y = element_text(size = 15))
avgs_p





#------------------------Spread among top genes--------------------------

top_num <- 100

# Load dataset avgs pro, created from correlation value where mating type was separated
corr_avgs_pro <- read.table("../../data_tables/all_genes/1_pearson_corr_avgs_with_pro.tsv", header = TRUE) 

# Avgs with promotor, 
veg <- subset(corr_avgs_pro, Tissue == "veg")
sex <- subset(corr_avgs_pro, Tissue == "sex")

veg_sa <- subset(veg, Allel == "sa")
veg_ba <- subset(veg, Allel == "BA")
sex_sa <- subset(sex, Allel == "sa")
sex_ba <- subset(sex, Allel == "BA")

#Add gene id as index 
rownames(veg_sa) <- veg_sa$Geneid
rownames(veg_ba) <- veg_ba$Geneid
rownames(sex_sa) <- sex_sa$Geneid
rownames(sex_ba) <- sex_ba$Geneid

top_veg_sa <- head(veg_sa[order(-veg_sa$Value, decreasing =TRUE),], top_num)
rownames(top_veg_sa) <- top_veg_sa$Geneid

spread <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Geneid", "Sample", "Value"))

# veg sa as top
for (gene in top_veg_sa$Geneid) {
  row <- data.frame(Geneid = gene, Sample = "veg_ba", Value = veg_ba[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "sex_sa", Value = sex_sa[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "sex_ba", Value = sex_ba[gene, "Value"])
  spread <- rbind(spread, row)
  # Top 
  row <- data.frame(Geneid = gene, Sample = "veg_sa", Value = top_veg_sa[gene, "Value"])
  spread <- rbind(spread, row)
}

a <- ggplot(spread, aes(x=Sample, y=Value, fill=Sample))+
  geom_violin() + 
  labs(title="Top 100 in veg sa", y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")


# veg ba as top
top_veg_ba <- head(veg_ba[order(-veg_ba$Value, decreasing =TRUE),], top_num)
rownames(top_veg_ba) <- top_veg_ba$Geneid

spread <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Geneid", "Sample", "Value"))

for (gene in top_veg_ba$Geneid) {
  row <- data.frame(Geneid = gene, Sample = "veg_sa", Value = veg_sa[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "sex_sa", Value = sex_sa[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "sex_ba", Value = sex_ba[gene, "Value"])
  spread <- rbind(spread, row)
  # Top 
  row <- data.frame(Geneid = gene, Sample = "veg_ba", Value = top_veg_ba[gene, "Value"])
  spread <- rbind(spread, row)
}

b <- ggplot(spread, aes(x=Sample, y=Value, fill=Sample))+
  geom_violin() + 
  labs(title="Top 100 in veg ba ", y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")


# sex sa as top
top_sex_sa <- head(sex_sa[order(-sex_sa$Value, decreasing =TRUE),], top_num)
rownames(top_sex_sa) <- top_sex_sa$Geneid

spread <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Geneid", "Sample", "Value"))

for (gene in top_sex_sa$Geneid) {
  row <- data.frame(Geneid = gene, Sample = "veg_sa", Value = veg_sa[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "veg_ba", Value = veg_ba[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "sex_ba", Value = sex_ba[gene, "Value"])
  spread <- rbind(spread, row)
  # Top 
  row <- data.frame(Geneid = gene, Sample = "sex_sa", Value = top_sex_sa[gene, "Value"])
  spread <- rbind(spread, row)
}

c <- ggplot(spread, aes(x=Sample, y=Value, fill=Sample))+
  geom_violin() + 
  labs(title="Top 100 in sex sa", y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")

# sex ba as top
top_sex_ba <- head(sex_ba[order(-sex_ba$Value, decreasing =TRUE),], top_num)
rownames(top_sex_ba) <- top_sex_ba$Geneid

spread <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Geneid", "Sample", "Value"))

for (gene in top_sex_ba$Geneid) {
  row <- data.frame(Geneid = gene, Sample = "veg_sa", Value = veg_sa[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "veg_ba", Value = veg_ba[gene, "Value"])
  spread <- rbind(spread, row)
  row <- data.frame(Geneid = gene, Sample = "sex_sa", Value = sex_sa[gene, "Value"])
  spread <- rbind(spread, row)
  # Top 
  row <- data.frame(Geneid = gene, Sample = "sex_ba", Value = top_sex_ba[gene, "Value"])
  spread <- rbind(spread, row)
}

d <- ggplot(spread, aes(x=Sample, y=Value, fill=Sample))+
  geom_violin() + 
  labs(title="Top 100 in sex sa", y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")

# Plot all in same plot
s <- 15
plot_row <- plot_grid(a + theme(text = element_text(size = s)),
                      b + theme(text = element_text(size = s)),
                      c + theme(text = element_text(size = s)),
                      d + theme(text = element_text(size = s)), labels = c("A", "B", "C", "D"))


plot_grid(plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# --------------------------intersection --------------------------------
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

b <- ggVennDiagram(top_genes, cex = 100)




