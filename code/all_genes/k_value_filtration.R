# 1. Divide between sex and veg
# 2. Use avg met for all 10 samples 
# 3. Do a distribution curve over methylation 
# 4. Divide genes between high, low and avg methylation
# 5. Divide k value filtration between high, low and middle levels. 

met <- read.table("../../data_tables/1_orthologous_genes_methylation.tsv", header = TRUE)

# With promotor region 
#met <- read.table("../../data_tables/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)

rownames(met) <- met$Geneid

samples_sex <- grep("_sex", colnames(met), value = TRUE)
samples_veg <- grep("_veg", colnames(met), value = TRUE)

avg_met <- setNames(data.frame(matrix(ncol = 2, nrow = 7109)), c("avg_met_sex", "avg_met_veg"))
avg_met <- cbind(met$Geneid, avg_met)

colnames(avg_met)[1] <- "Geneid"
rownames(avg_met) <- avg_met$Geneid

#gene <- "Ntsc3922"
for (geneid in met$Geneid) {
  sum_value <- 0 
  for (sample in samples_sex) {
    sum_value <- sum_value + met[geneid, sample]
  }
  avg_value <- sum_value/10
  avg_met[geneid, "avg_met_sex"] <- avg_value
  
  # veg 
  sum_value <- 0 
  for (sample in samples_veg) {
    sum_value <- sum_value + met[geneid, sample]
  }
  avg_value <- sum_value/10
  avg_met[geneid, "avg_met_veg"] <- avg_value
}

# Distributions

# sex
mean_sex <- mean(avg_met$avg_met_sex, na.rm=TRUE)
sd_sex <- sd(avg_met$avg_met_sex, na.rm=TRUE)

two_sd <- mean_sex + 2*sd_sex
mtwo_sd <- mean_sex - 2*sd_sex

top_genes_sex <- list()
mid_genes_sex <- list()
bot_genes_sex <- list()

for (geneid in avg_met$Geneid) {
  avg_met_value <- avg_met[geneid, "avg_met_sex"]
  
  if (is.na(avg_met_value) == FALSE ) {
    # top gene
    if (avg_met_value > two_sd) {
      top_genes_sex <- append(top_genes_sex, geneid)
    }
    # Bottom gene
    else if (avg_met_value < mtwo_sd) {
      bot_genes_sex <- append(bot_genes_sex, geneid)
    }
    else if (avg_met_value < two_sd & avg_met_value > mtwo_sd) {
      mid_genes_sex <- append(mid_genes_sex, geneid)
    }
    else {
      print("woaw")
    }
  }
}


# veg
mean_veg <- mean(avg_met$avg_met_veg, na.rm=TRUE)
sd_veg <- sd(avg_met$avg_met_veg, na.rm=TRUE)

two_sd <- mean_veg + 2*sd_veg
mtwo_sd <- mean_veg - 2*sd_veg

top_genes_veg <- list()
mid_genes_veg <- list()
bot_genes_veg <- list()

for (geneid in avg_met$Geneid) {
  avg_met_value <- avg_met[geneid, "avg_met_veg"]
  
  if (is.na(avg_met_value) == FALSE ) {
    # top gene
    if (avg_met_value > two_sd) {
      top_genes_veg <- append(top_genes_veg, geneid)
    }
    # Bottom gene
    else if (avg_met_value < mtwo_sd) {
      bot_genes_veg <- append(bot_genes_veg, geneid)
    }
    else if (avg_met_value < two_sd & avg_met_value > mtwo_sd) {
      mid_genes_veg <- append(mid_genes_veg, geneid)
    }
    else {
      print("woaw")
    }
  }
}


# K value filtration 

# Without promotor region
df <- read.table("../../data_tables/all_genes/1_pearson_corr_avgs_non_allel_sep.tsv", header = TRUE)
#df <- read.table("../../data_tables/all_genes/1_pearson_corr_rep_non_allel_sep.tsv", header = TRUE)
#df <- read.table("../../data_tables/all_genes/1_spearman_corr_avgs_non_allel_sep.tsv", header = TRUE)
#df <- read.table("../../data_tables/all_genes/1_spearman_corr_rep_non_allel_sep.tsv", header = TRUE)

# Diff tables with promotor region, obs change met table to row 11. 
#df <- read.table("../../data_tables/all_genes/1_pearson_corr_avgs_with_pro_non_allel_sep.tsv", header = TRUE)
#df <- read.table("../../data_tables/all_genes/1_pearson_corr_rep_with_pro_non_allel_sep.tsv", header = TRUE)
#df <- read.table("../../data_tables/all_genes/1_spearman_corr_avgs_with_pro_non_allel_sep.tsv", header = TRUE)
#df <- read.table("../../data_tables/all_genes/1_spearman_corr_rep_with_pro_non_allel_sep.tsv", header = TRUE)

mean_val <- mean(df$k_value, na.rm = TRUE)
sd_val <- sd(df$k_value, na.rm = TRUE)

min_val <- mean_val - sd_val
max_val <- mean_val + sd_val

df_filt1 <- subset(df, k_value < min_val)
df_filt2 <- subset (df, k_value > max_val)

df_filtered <- rbind(df_filt1, df_filt2)
rownames(df_filtered) <- df_filtered$Geneid

#sex
df_filtered_sex <- subset(df_filtered, Tissue == "sex" )
df_filtered_top_sex <- subset(df_filtered_sex, Geneid %in% top_genes_sex)
df_filtered_mid_sex <- subset(df_filtered_sex, Geneid %in% mid_genes_sex)
df_filtered_bot_sex <- subset(df_filtered_sex, Geneid %in% bot_genes_sex)

#veg
df_filtered_veg <- subset(df_filtered, Tissue == "veg" )
df_filtered_top_veg <- subset(df_filtered_veg, Geneid %in% top_genes_veg)
df_filtered_mid_veg <- subset(df_filtered_veg, Geneid %in% mid_genes_veg)
df_filtered_bot_veg <- subset(df_filtered_veg, Geneid %in% bot_genes_veg)

df_filtered_top <- rbind(df_filtered_top_sex, df_filtered_top_veg)
df_filtered_mid <- rbind(df_filtered_mid_sex, df_filtered_mid_veg)
df_filtered_bot <- rbind(df_filtered_bot_sex, df_filtered_bot_veg)


a <- ggplot(df_filtered_top, aes(x=Tissue, y=Value, fill=Tissue))+
            geom_violin() + 
            labs(title="High methylated genes", y = "Correlation value")+
            geom_boxplot(width=0.1, fill = "white") + ylim(c(-1,1)) 
  

b <- ggplot(df_filtered_mid, aes(x=Tissue, y=Value, fill=Tissue))+
            geom_violin() + 
            labs(title="Avgerage methylated genes", y = "Correlation value")+
            geom_boxplot(width=0.1, fill = "white") + + ylim(c(-1,1))

c <- ggplot(df_filtered_bot, aes(x=Tissue, y=Value, fill=Tissue))+
            geom_violin() + 
            labs(title="Low methylated genes", y = "Correlation value")+
            geom_boxplot(width=0.1, fill = "white") + ylim(c(-1,1))

plot_row <- plot_grid(a, b, c, labels = c("A", "B", "C"))

title <- ggdraw() + 
  draw_label(
    "Pearson correation values for filtered avgs genes",
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






