# Filter correlation tables based on k value 

# Packages 
library("dplyr")
library(ggplot2)
library(cowplot)
library(ggVennDiagram)

# Run topgenes.R here first 
# Run midgenes.R here first 
# Run botgenes.R here first 
# Run topgenes_veg.R here first 
# Run midgenes_veg.R here first 
# Run botgenes_veg.R here first

# sex tissue
top_genes_sex_sa <- Reduce(intersect, list(top_crassa_sa_sex, 
                                           top_L1_sa_sex,
                                           top_L10_sa_sex, 
                                           top_L6_sa_sex, 
                                           top_sito_sa_sex))

top_genes_sex_ba <- Reduce(intersect, list(top_crassa_BA_sex, 
                                           top_L1_BA_sex,
                                           top_L10_BA_sex, 
                                           top_L6_BA_sex,  
                                           top_sito_BA_sex))

top_genes_sex <- union(top_genes_sex_ba, top_genes_sex_sa)


mid_genes_sex_sa <- Reduce(intersect, list(mid_crassa_sa_sex, 
                                           mid_L1_sa_sex,
                                           mid_L10_sa_sex, 
                                           mid_L6_sa_sex, 
                                           mid_sito_sa_sex))

mid_genes_sex_ba <- Reduce(intersect, list(mid_crassa_BA_sex, 
                                           mid_L1_BA_sex,
                                           mid_L10_BA_sex, 
                                           mid_L6_BA_sex, 
                                           mid_sito_BA_sex))

mid_genes_sex <- union(mid_genes_sex_ba, mid_genes_sex_sa)

bot_genes_sex_sa <- Reduce(intersect, list(bot_crassa_sa_sex, 
                                           bot_L1_sa_sex,
                                           bot_L10_sa_sex, 
                                           bot_L6_sa_sex, 
                                           bot_sito_sa_sex))

bot_genes_sex_ba <- Reduce(intersect, list(bot_crassa_BA_sex, 
                                           bot_L1_BA_sex,
                                           bot_L10_BA_sex,  
                                           bot_L6_BA_sex,  
                                           bot_sito_BA_sex))

bot_genes_sex <- union(bot_genes_sex_ba, bot_genes_sex_sa)

# veg tissue
top_genes_veg <- Reduce(intersect, list(top_crassa_BA_veg, top_crassa_sa_veg, 
                                        top_L1_BA_veg, top_L1_BA_veg,
                                        top_L10_BA_veg, top_L10_BA_veg, 
                                        top_L6_BA_veg, top_L6_sa_veg, 
                                        top_sito_BA_veg, top_sito_sa_veg))


mid_genes_veg <- Reduce(intersect, list(mid_crassa_BA_veg, mid_crassa_sa_veg, 
                                        mid_L1_BA_veg, mid_L1_BA_veg,
                                        mid_L10_BA_veg, mid_L10_BA_veg, 
                                        mid_L6_BA_veg, mid_L6_sa_veg, 
                                        mid_sito_BA_veg, mid_sito_sa_veg))

bot_genes_veg <- Reduce(intersect, list(bot_crassa_BA_veg, bot_crassa_sa_veg, 
                                        bot_L1_BA_veg, bot_L1_BA_veg,
                                        bot_L10_BA_veg, bot_L10_BA_veg, 
                                        bot_L6_BA_veg, bot_L6_sa_veg, 
                                        bot_sito_BA_veg, bot_sito_sa_veg))





# CAN only filtrate based  on the tabels wheer i havent seperated on mating type ??? 
# Mating type problem matching geneids 

# Pearson -----------------------------------------------
df <- read.table("new/1_pearson_corr_avgs.tsv", header = TRUE)

mean_val <- mean(df$k_value, na.rm = TRUE)
sd_val <- sd(df$k_value, na.rm = TRUE)

min_val <- mean_val - sd_val
max_val <- mean_val + sd_val

df_filt1 <- subset(df, k_value < min_val)
df_filt2 <- subset (df, k_value > max_val)

df_filtered <- rbind(df_filt1, df_filt2)
rownames(df_filtered) <- df_filtered$Geneid

# bot genes is empty 

#top genes
top_genes_df <- rbind(df_filtered["Ntsc3293",], df_filtered["Ntsc7137",]) # dont exist in filtered (?) 

#mid genes

mid_genes_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Geneid", "Allel", "Tissue", "Value", "k_value"))
for (geneid in df_filtered$Geneid) {
  if (geneid %in% mid_genes) {
    mid_genes_df <- rbind(mid_genes_df, df_filtered[geneid, ])
    
  }
}




ggplot(df_filtered, aes(x=Tissue, y=Value, fill=Tissue))+
  geom_violin() + 
  labs(title="Using average of replicas", y = "Correlation value")+
  geom_boxplot(width=0.1, fill = "white")


# All diff tables 
read.table("new/1_pearson_corr_avgs.tsv", header = TRUE)
read.table("new/1_pearson_corr_avgs_non_allel_sep.tsv", header = TRUE)
read.table("new/1_pearson_corr_avgs_with_pro.tsv", header = TRUE)
read.table("new/1_pearson_corr_avgs_with_pro_non_allel_sep.tsv", header = TRUE)

read.table("new/1_pearson_corr_rep.tsv", header = TRUE)
read.table("new/1_pearson_corr_rep_non_allel_sep.tsv", header = TRUE)
read.table("new/1_pearson_corr_rep_with_pro.tsv", header = TRUE)
read.table("new/1_pearson_corr_rep_with_pro_non_allel_sep.tsv", header = TRUE)

read.table("new/1_spearman_corr_avgs.tsv", header = TRUE)
read.table("new/1_spearman_corr_avgs_non_allel_sep.tsv", header = TRUE)
read.table("new/1_spearman_corr_avgs_with_pro.tsv", header = TRUE)
read.table("new/1_spearman_corr_avgs_with_pro_non_allel_sep.tsv", header = TRUE)

read.table("new/1_spearman_corr_rep.tsv", header = TRUE)
read.table("new/1_spearman_corr_rep_non_allel_sep.tsv", header = TRUE)
read.table("new/1_spearman_corr_rep_with_pro.tsv", header = TRUE)
read.table("new/1_spearman_corr_rep_with_pro_non_allel_sep.tsv", header = TRUE)




