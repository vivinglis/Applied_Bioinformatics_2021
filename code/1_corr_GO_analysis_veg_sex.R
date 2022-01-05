#Correlation analysis
library("stringi")
library("ggplot2")
library("tidyverse")

####Read in tables
pearson_corr_avgs <- read.table("../Applied_Bioinformatics_2021/k_data_tables/1_pearson_corr_avgs_non_allel_sep.tsv", header = TRUE)
pearson_corr_avgs_with_pro <- read.table("../Applied_Bioinformatics_2021/k_data_tables/1_pearson_corr_avgs_with_pronon_allel_sep.tsv", header = TRUE)
pearson_corr_rep <- read.table("../Applied_Bioinformatics_2021/k_data_tables/1_pearson_corr_rep_non_allel_sep.tsv", header = TRUE)
pearson_corr_rep_with_pro <- read.table("../Applied_Bioinformatics_2021/k_data_tables/1_pearson_corr_rep_with_pronon_allel_sep.tsv", header = TRUE)
GO_terms <- read.delim("1_filtered_GO_terms_per_gene_4_or_more_genes_per_GO.tsv")

corrs <- list(pearson_corr_avgs, pearson_corr_avgs_with_pro, pearson_corr_rep, pearson_corr_rep_with_pro)
avg_df <- as.data.frame(matrix(nrow=0, ncol = 7))
top_df <- list()
for (corr in corrs){
  # Dividing data into subsets
  pearson_corr_veg <- corr %>% filter(Tissue == "veg")
  pearson_corr_sex <- corr %>% filter(Tissue == "sex")
  rownames(pearson_corr_sex) <- pearson_corr_sex$Geneid
  rownames(pearson_corr_veg) <- pearson_corr_veg$Geneid
  GO_sex <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("GO_term","Avg_Corr","SDev_Corr", "Avg_k_val", "SDev_k_val"))
  GO_veg <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("GO_term","Avg_Corr","SDev_Corr", "Avg_k_val", "SDev_k_val"))
  for (i in 1:ncol(GO_terms)){
    per_GO_sex <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("GO_term","Corr", "k_val"))
    per_GO_veg <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("GO_term","Corr", "k_val"))

    for (j in 1:nrow(GO_terms)) {
      if (is.na(GO_terms[j, i])) {
        break
      }
      if (!is.na(pearson_corr_sex[GO_terms[j, i],"Value"])){
        sex <- data.frame(colnames(GO_terms)[i], pearson_corr_sex[GO_terms[j, i],"Value"], pearson_corr_sex[GO_terms[j, i], "k_value"])
        names(sex) <- c("GO_term","Corr", "k_val")
        per_GO_sex <- rbind(per_GO_sex, sex)
        
      }
      if (!is.na(pearson_corr_veg[GO_terms[j, i], "Value"])){
        veg <- data.frame(colnames(GO_terms)[i], pearson_corr_veg[GO_terms[j, i], "Value"], pearson_corr_veg[GO_terms[j, i], "k_value"])
        names(veg) <- c("GO_term","Corr", "k_val")
        per_GO_veg <- rbind(per_GO_veg, veg)
      }
     
    }
    temp_GO_sex <- data.frame("Sex", colnames(GO_terms)[i], mean(per_GO_sex$Corr), sqrt(var(per_GO_sex$Corr)), mean(per_GO_sex$k_val), sqrt(var(per_GO_sex$k_val)))
    temp_GO_veg <- data.frame("Veg", colnames(GO_terms)[i], mean(per_GO_veg$Corr), sqrt(var(per_GO_veg$Corr)), mean(per_GO_veg$k_val), sqrt(var(per_GO_veg$k_val)))
    names(temp_GO_sex) <- c("Tissue", "GO_term", "Avg_Corr", "SDev_Corr", "Avg_k_val", "SDev_k_val")
    names(temp_GO_veg) <- c("Tissue", "GO_term", "Avg_Corr", "SDev_Corr", "Avg_k_val", "SDev_k_val")
    GO_sex <- rbind(GO_sex, temp_GO_sex, temp_GO_veg)
  }
  # Plotting

  top_df <- list(top_df, GO_sex)
}


corr_avgs <- as.data.frame(top_df[[1]][[1]][[1]][[2]])
corr_avgs_pro <- as.data.frame(top_df[[1]][[1]][[2]])
corr_reps <- as.data.frame(top_df[[1]][[2]])
corr_reps_pro <- as.data.frame(top_df[[2]])

corr_veg_avgs <- corr_avgs %>% filter(Tissue == "Veg")
corr_sex_avgs <- corr_avgs %>% filter(Tissue == "Sex")
corr_veg_avgs_pro <- corr_avgs_pro %>% filter(Tissue == "Veg")
corr_sex_avgs_pro <- corr_avgs_pro %>% filter(Tissue == "Sex")
corr_veg_reps <- corr_reps %>% filter(Tissue == "Veg")
corr_sex_reps <- corr_reps %>% filter(Tissue == "Sex")
corr_veg_reps_pro <- corr_reps_pro %>% filter(Tissue == "Veg")
corr_sex_reps_pro <- corr_reps_pro %>% filter(Tissue == "Sex")

top_num <- 10
std_threshold <- 1
corr_threshold <- 1
# #sample_name <- "sa_sex"
# 
# top_corr_avgs <- head(corr_avgs[order(corr_avgs$Avg_Corr),], top_num)
# top_corr_avgs_pro <- head(corr_avgs_pro[order(corr_avgs_pro$Avg_Corr),], top_num)
# top_corr_reps <- head(corr_reps[order(corr_reps$Avg_Corr),], top_num)
# top_corr_reps_pro <- head(corr_reps_pro[order(corr_reps_pro$Avg_Corr),], top_num)

# top_corr_avgs <- head(corr_avgs[order(corr_avgs$SDev_Corr),], top_num)
# top_corr_avgs_pro <- head(corr_avgs_pro[order(corr_avgs_pro$SDev_Corr),], top_num)
# top_corr_reps <- head(corr_reps[order(corr_reps$SDev_Corr),], top_num)
# top_corr_reps_pro <- head(corr_reps_pro[order(corr_reps_pro$SDev_Corr),], top_num)

corr_avgs <- corr_avgs[corr_avgs$Avg_Corr<corr_threshold ,]
corr_avgs_pro <- corr_avgs_pro[corr_avgs_pro$Avg_Corr<corr_threshold ,]
corr_reps <- corr_reps[corr_reps$Avg_Corr<corr_threshold ,]
corr_reps_pro <- corr_reps_pro[corr_reps_pro$Avg_Corr<corr_threshold ,]

top_corr_avgs <- corr_avgs[corr_avgs$SDev_Corr<std_threshold ,]
top_corr_avgs_pro <- corr_avgs_pro[corr_avgs_pro$SDev_Corr<std_threshold ,]
top_corr_reps <- corr_reps[corr_reps$SDev_Corr<std_threshold ,]
top_corr_reps_pro <- corr_reps_pro[corr_reps_pro$SDev_Corr<std_threshold ,]

top_corr_avgs <- na.omit(top_corr_avgs)
top_corr_avgs_pro <- na.omit(top_corr_avgs_pro)
top_corr_reps <- na.omit(top_corr_reps)
top_corr_reps_pro <- na.omit(top_corr_reps_pro)

# sample_BA_veg <- grep("_BA_veg", colnames(gene_expr), value = TRUE)
# sample_BA_sex <- grep("_BA_sex", colnames(gene_expr), value = TRUE)
# sample_sa_veg <- grep("_sa_veg", colnames(gene_expr), value = TRUE)
# sample_sa_sex <- grep("_sa_sex", colnames(gene_expr), value = TRUE)
# 
# samples <- list(sample_BA_veg, sample_BA_sex, sample_sa_veg, sample_sa_sex)
# all_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Gene", "lineage", "met", "expr"))
# for (gene in GO_terms$peptidyl.prolyl.cis.trans.isomerase.activity){
#   if (is.na(gene)){
#     break
#   }
#   for (sample in samples){
#     per_gene <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Gene", "lineage","met", "expr"))
#     
#     for (lineage in sample){
#       new <- data.frame(Gene = gene, lineage = lineage, met=met[gene,lineage], expr = gene_expr[gene, lineage])
#       per_gene <- rbind(per_gene, new)
#     }
#     
#     all_df <- rbind(all_df, per_gene)
#   }
# }


p_avgs <- ggplot(data=top_corr_avgs, aes(x=Tissue, y=SDev_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  
  labs(title="Gene expression values averaged")



p_avgs_pro <- ggplot(data=top_corr_avgs_pro, aes(x=Tissue, y=SDev_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")
#+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  
  #labs(title="Gene expression values averaged with promotor region")



p_reps <- ggplot(data=top_corr_reps, aes(x=Tissue, y=SDev_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  labs(title="All gene expression replicates")



p_reps_pro <- ggplot(data=top_corr_reps_pro, aes(x=Tissue, y=SDev_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  
  labs(title="All gene expression replicates with promotor region")

library(cowplot)

p_row <- plot_grid(p_avgs, p_avgs_pro, p_reps, p_reps_pro, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste0("Standard deviation of pearson correlation between methylation level and gene expression for the genes in GO terms containing 4 or more genes"),# where SDev < ",std_threshold),
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

p_avgs_avg <- ggplot(data=top_corr_avgs, aes(x=Tissue, y=Avg_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +

  labs(title="Gene expression values averaged")



p_avgs_pro_avg <- ggplot(data=top_corr_avgs_pro, aes(x=Tissue, y=Avg_Corr, fill=Tissue, label=GO_term)) +
  #geom_violin() +
  #geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  #geom_point()
  #geom_text(size=6) +
  geom_label(size=10) 
  #+

  #labs(title="Gene expression values averaged with promotor region")



p_reps_avg <- ggplot(data=top_corr_reps, aes(x=Tissue, y=Avg_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +
  labs(title="All gene expression replicates")



p_reps_pro_avg <- ggplot(data=top_corr_reps_pro, aes(x=Tissue, y=Avg_Corr, fill=Tissue, label=GO_term)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  geom_text(size=3) +
  geom_label() +

  labs(title="All gene expression replicates with promotor region")


p_row_avg <- plot_grid(p_avgs_avg, p_avgs_pro_avg, p_reps_avg, p_reps_pro_avg, labels = "AUTO")

title_avg <- ggdraw() +
  draw_label(
    paste0("Average pearson correlation between methylation level and gene expression for the genes in GO terms containing 4 or more genes"),# where SDev < ",std_threshold),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_grid_avg <- plot_grid(
  title_avg, p_row_avg,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)