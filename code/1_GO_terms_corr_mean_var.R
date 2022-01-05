#Correlation analysis
library("stringi")
library("ggplot2")
library("tidyverse")
# Read in tables
#pearson_corr_avgs <- read.table("1_pear_corr_GO_terms_avgs.tsv", header = TRUE) 
pearson_corr_avgs <- read.table("1_2_pearson_corr_avgs.tsv", header = TRUE) 
pearson_corr_avgs_with_pro <- read.table("1_2_pearson_corr_avgs_with_pro.tsv", header = TRUE) 
pearson_corr_rep <- read.table("../Applied_Bioinformatics_2021/data_tables/1_2_pearson_corr_rep.tsv", header = TRUE) 
pearson_corr_rep_with_pro <- read.table("../Applied_Bioinformatics_2021/data_tables/1_2_pearson_corr_rep_with_pro.tsv", header = TRUE) 

#NEED TO REMOVE NAs
GO_terms  <- read.delim("1_filtered_GO_terms_per_gene_4_or_more_genes_per_GO.tsv")

corrs <- list(pearson_corr_avgs, pearson_corr_avgs_with_pro, pearson_corr_rep, pearson_corr_rep_with_pro)

sample_names <- c("BA_veg", "sa_veg", "sa_sex", "BA_sex")
res_list <- list()

for (corr in corrs){
  corr_complete <- na.omit(corr)
  corr_BA <- corr_complete %>% filter(Allel == "BA")
  corr_sa <- corr_complete %>% filter(Allel== "sa")
  corr_BA_veg <- corr_BA %>% filter(Tissue== "veg")
  corr_sa_sex <- corr_sa %>% filter(Tissue == "sex")
  corr_sa_veg <- corr_sa %>% filter(Tissue== "veg")
  corr_BA_sex <- corr_BA %>% filter(Tissue == "sex")

  rownames(corr_BA_veg) <- corr_BA_veg$Geneid
  rownames(corr_sa_veg) <- corr_sa_veg$Geneid
  rownames(corr_sa_sex) <- corr_sa_sex$Geneid
  rownames(corr_BA_sex) <- corr_BA_sex$Geneid
  samples <- list(BA_veg = corr_BA_veg, sa_veg = corr_sa_veg, BA_sex = corr_BA_sex, sa_sex = corr_sa_sex)
  res_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Sample", "GO_term", "Avg_corr", "std_corr"))
  num = 1
  for (sample in samples){
      for (i in 1:ncol(GO_terms)) {
        per_GO <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),  "corr")
        
        for (j in 1:nrow(GO_terms)) {
          if (is.na(GO_terms[j, i])) {
            break
          }
          new <- data.frame( sample[GO_terms[j, i], 4])
          per_GO <- rbind(per_GO, new)                
        }
        mean_val <- mean(per_GO[,1], na.rm = TRUE)
        std_val <- sqrt(var(per_GO[,1], y=NULL, na.rm = TRUE))
        corr_new <- data.frame(Sample = sample_names[num], GO_term = colnames(GO_terms)[i], Avg_corr = mean_val, std_corr = std_val)
        res_df <- rbind(res_df, corr_new)
      }
    num = num + 1
  }
  res_list <- list(res_list, res_df)
}


data_avgs <- as.data.frame(res_list[[1]][[1]][[1]][[2]])
data_avgs$Sample <- as.factor(data_avgs$Sample)
data_avgs$std_corr <- as.double(data_avgs$std_corr)

p_avgs <- ggplot(data=data_avgs, aes(x=Sample, y=std_corr, fill=Sample)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged")


data_avgs_pro <- as.data.frame(res_list[[1]][[1]][[2]])
#data_avgs_pro$GO_terms <- as.data.frame(res_list[[1]][[1]][[3]])



data_avgs_pro$Sample <- as.factor(data_avgs_pro$Sample)
data_avgs_pro$std_corr <- as.double(data_avgs_pro$std_corr)

p_avgs_pro <- ggplot(data=data_avgs_pro, aes(x=Sample, y=std_corr, fill=Sample)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged with promotor region")


data_reps <- as.data.frame(res_list[[1]][[2]])
#data_reps$GO_terms <- as.data.frame(res_list[[1]][[3]])



data_reps$Sample <- as.factor(data_reps$Sample)
data_reps$std_corr <- as.double(data_reps$std_corr)

p_reps <- ggplot(data=data_reps, aes(x=Sample, y=std_corr, fill=Sample)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg"))+
  labs(title="All gene expression replicates")

data_reps_pro <- as.data.frame(res_list[[2]])
#data_reps_pro$GO_terms <- as.data.frame(res_list[[3]])


data_reps_pro$Sample <- as.factor(data_reps_pro$Sample)
data_reps_pro$std_corr <- as.double(data_reps_pro$std_corr)

p_reps_pro <- ggplot(data=data_reps_pro, aes(x=Sample, y=std_corr, fill=Sample)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="All gene expression replicates with promotor region")

library(cowplot)

p_row <- plot_grid(p_avgs, p_avgs_pro, p_reps, p_reps_pro, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste("Standard deviation of the correlation between methylation level and gene expression for the genes in GO terms containing 4 or more genes"),
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
# 
# data_avgs <- as.data.frame(res_list[[1]][[1]][[1]][[2]])
# #data_avgs$GO_terms <- as.data.frame(res_list[[1]][[1]][[1]][[3]])
# 
# 
# 
# data_avgs$Sample <- as.factor(data_avgs$Sample)
# data_avgs$Avg_corr <- as.double(data_avgs$Avg_corr)
# 
# p_avgs <- ggplot(data=data_avgs, aes(x=Sample, y=Avg_corr, fill=Sample)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
#   labs(title="Gene expression values averaged")
# 
# 
# data_avgs_pro <- as.data.frame(res_list[[1]][[1]][[2]])
# #data_avgs_pro$GO_terms <- as.data.frame(res_list[[1]][[1]][[3]])
# 
# 
# 
# data_avgs_pro$Sample <- as.factor(data_avgs_pro$Sample)
# data_avgs_pro$Avg_corr <- as.double(data_avgs_pro$Avg_corr)
# 
# p_avgs_pro <- ggplot(data=data_avgs_pro, aes(x=Sample, y=Avg_corr, fill=Sample)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
#   labs(title="Gene expression values averaged with promotor region")
# 
# 
# data_reps <- as.data.frame(res_list[[1]][[2]])
# #data_reps$GO_terms <- as.data.frame(res_list[[1]][[3]])
# 
# 
# 
# data_reps$Sample <- as.factor(data_reps$Sample)
# data_reps$Avg_corr <- as.double(data_reps$Avg_corr)
# 
# p_reps <- ggplot(data=data_reps, aes(x=Sample, y=Avg_corr, fill=Sample)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg"))+
#   labs(title="All gene expression replicates")
# 
# data_reps_pro <- as.data.frame(res_list[[2]])
# #data_reps_pro$GO_terms <- as.data.frame(res_list[[3]])
# 
# 
# data_reps_pro$Sample <- as.factor(data_reps_pro$Sample)
# data_reps_pro$Avg_corr <- as.double(data_reps_pro$Avg_corr)
# 
# p_reps_pro <- ggplot(data=data_reps_pro, aes(x=Sample, y=Avg_corr, fill=Sample)) +
#   geom_violin() +
#   geom_boxplot(width=0.1, fill="white")+
#   scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
#   labs(title="All gene expression replicates with promotor region")
# 
# library(cowplot)
# 
# p_row <- plot_grid(p_avgs, p_avgs_pro, p_reps, p_reps_pro, labels = "AUTO")
# 
# title <- ggdraw() +
#   draw_label(
#     paste("Average correlation between methylation level and gene expression for the genes in GO terms containing 4 or more genes"),
#     fontface = 'bold',
#     x = 0,
#     hjust = 0
#   ) +
#   theme(
#     # add margin on the left of the drawing canvas,
#     # so title is aligned with left edge of first plot
#     plot.margin = margin(0, 0, 0, 17)
#   )
# p_grid <- plot_grid(
#   title, p_row,
#   ncol = 1,
#   # rel_heights values control vertical title margins
#   rel_heights = c(0.1, 1)
# )
# 



# num = "avg_and_std"
# sample_name = "4"
# 
# write.table(data_avgs, paste0("1_data_for_corr_",num,"_GO_terms_filtered_",sample_name,"_avgs.tsv"), quote=FALSE, sep='\t', row.names = FALSE)
# write.table(data_avgs_pro, paste0("1_data_for_corr_",num,"_GO_terms_filtered_",sample_name,"_avgs_pro.tsv"), quote=FALSE, sep='\t', row.names = FALSE)
# write.table(data_reps, paste0("1_data_for_corr_",num,"_GO_terms_filtered_",sample_name,"_reps.tsv"), quote=FALSE, sep='\t', row.names = FALSE)
# write.table(data_reps_pro, paste0("1_data_for_corr_",num,"_GO_terms_filtered_",sample_name,"_reps_pro.tsv"), quote=FALSE, sep='\t', row.names = FALSE)
