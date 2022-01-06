#Correlation analysis
library("stringi")
library("ggplot2")
library("tidyverse")
# Read in tables
#corr_avgs <- read.table("1_pear_corr_GO_terms_avgs.tsv", header = TRUE) 
corr_avgs <- read.table("1_data_for_corr_avg_and_std_GO_terms_filtered_4_avgs.tsv", header = TRUE) 
corr_avgs_pro <- read.table("1_data_for_corr_avg_and_std_GO_terms_filtered_4_avgs_pro.tsv", header = TRUE) 
corr_reps <- read.table("1_data_for_corr_avg_and_std_GO_terms_filtered_4_reps.tsv", header = TRUE) 
corr_reps_pro <- read.table("1_data_for_corr_avg_and_std_GO_terms_filtered_4_reps_pro.tsv", header = TRUE) 

top_num <- 100
std_threshold <- 0.25
corr_threshold <- -0.25
# #sample_name <- "sa_sex"
# top_corr_avgs <- head(corr_avgs[order(corr_avgs$std_corr),], top_num)
# top_corr_avgs_pro <- head(corr_avgs_pro[order(corr_avgs_pro$std_corr),], top_num)
# top_corr_reps <- head(corr_reps[order(corr_reps$std_corr),], top_num)
# top_corr_reps_pro <- head(corr_reps_pro[order(corr_reps_pro$std_corr),], top_num)

top_corr_avgs <- corr_avgs[corr_avgs$Avg_corr<corr_threshold ,]
top_corr_avgs_pro <- corr_avgs_pro[corr_avgs_pro$Avg_corr<corr_threshold ,]
top_corr_reps <- corr_reps[corr_reps$Avg_corr<corr_threshold ,]
top_corr_reps_pro <- corr_reps_pro[corr_reps_pro$Avg_corr<corr_threshold ,]

top_corr_avgs <- top_corr_avgs[top_corr_avgs$std_corr<std_threshold ,]
top_corr_avgs_pro <- top_corr_avgs_pro[top_corr_avgs_pro$std_corr<std_threshold ,]
top_corr_reps <- top_corr_reps[top_corr_reps$std_corr<std_threshold ,]
top_corr_reps_pro <- top_corr_reps_pro[top_corr_reps_pro$std_corr<std_threshold ,]

sample_BA_veg <- grep("_BA_veg", colnames(gene_expr), value = TRUE)
sample_BA_sex <- grep("_BA_sex", colnames(gene_expr), value = TRUE)
sample_sa_veg <- grep("_sa_veg", colnames(gene_expr), value = TRUE)
sample_sa_sex <- grep("_sa_sex", colnames(gene_expr), value = TRUE)

samples <- list(sample_BA_veg, sample_BA_sex, sample_sa_veg, sample_sa_sex)
all_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Gene", "lineage", "met", "expr"))
for (gene in GO_terms$peptidyl.prolyl.cis.trans.isomerase.activity){
  if (is.na(gene)){
    break
  }
  for (sample in samples){
    per_gene <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Gene", "lineage","met", "expr"))
  
    for (lineage in sample){
      new <- data.frame(Gene = gene, lineage = lineage, met=met[gene,lineage], expr = gene_expr[gene, lineage])
      per_gene <- rbind(per_gene, new)
    }
  
    all_df <- rbind(all_df, per_gene)
  }
}

# 
# p_avgs <- ggplot(data=top_corr_avgs, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
#   labs(title="Gene expression values averaged")
# 
# 
# 
# p_avgs_pro <- ggplot(data=top_corr_avgs_pro, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
#   labs(title="Gene expression values averaged with promotor region")
# 
# 
# 
# p_reps <- ggplot(data=top_corr_reps, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg"))+
#   labs(title="All gene expression replicates")
# 
# 
# 
# p_reps_pro <- ggplot(data=top_corr_reps_pro, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
#   labs(title="All gene expression replicates with promotor region")

library(cowplot)

p_row <- plot_grid(p_avgs, p_avgs_pro, p_reps, p_reps_pro, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste0("Average correlation between methylation level and gene expression for the genes in GO terms containing 4 or more genes where SDev < ",std_threshold),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_grid <- plot_grid(
  title, p_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

p_avgs <- ggplot(data=top_corr_avgs, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged")



p_avgs_pro <- ggplot(data=top_corr_avgs_pro, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged with promotor region")



p_reps <- ggplot(data=top_corr_reps, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg"))+
  labs(title="All gene expression replicates")



p_reps_pro <- ggplot(data=top_corr_reps_pro, aes(x=Sample, y=Avg_corr, fill=Sample, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="All gene expression replicates with promotor region")
