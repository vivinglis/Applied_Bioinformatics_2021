
#Reading data
pearson_corr_avgs <- read.table('1_2_pearson_corr_avgs.tsv', header = TRUE) 
pearson_corr_avgs_with_pro <- read.table("1_2_pearson_corr_avgs_with_pro.tsv", header = TRUE)
pearson_corr_rep <- read.table("../Applied_Bioinformatics_2021/data_tables/1_2_pearson_corr_rep.tsv", header = TRUE)
pearson_corr_rep_with_pro <- read.table("../Applied_Bioinformatics_2021/data_tables/1_2_pearson_corr_rep_with_pro.tsv", header = TRUE)

corrs <- list(pearson_corr_avgs, pearson_corr_avgs_with_pro, pearson_corr_rep, pearson_corr_rep_with_pro)
top_df <- list()

for (corr in corrs){
  #Filtering data 
  corr_BA <- corr %>% filter(Allel == "BA")
  corr_sa <- corr %>% filter(Allel== "sa")
  corr_BA_veg <- corr_BA %>% filter(Tissue== "veg")
  corr_sa_sex <- corr_sa %>% filter(Tissue == "sex")
  corr_sa_veg <- corr_sa %>% filter(Tissue== "veg")
  corr_BA_sex <- corr_BA %>% filter(Tissue == "sex")
  
  
  #Collecting all data in top_df
  data_df <- as.data.frame(corr_BA_sex$Value)
  data_df$sa_sex <- corr_sa_sex$Value
  data_df$BA_veg <- corr_BA_veg$Value
  data_df$sa_veg <- corr_sa_veg$Value
  colnames(data_df)[1] <- "BA_sex"
  rownames(data_df) <-corr_BA_veg$Geneid
  top_df <- list(top_df, data_df)
}

library(ggplot2)
library(ggpmisc)
my.formula <- y ~ x

#All individual plots for avgs
data_avgs <- as.data.frame(top_df[[1]][[1]][[1]][[2]])

p_avgs_sex <- ggplot(data=data_avgs, aes(x=BA_sex, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  labs(title="Sex plot")

p_avgs_BA <- ggplot(data=data_avgs, aes(x=BA_sex, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="BA plot")

p_avgs_sa <- ggplot(data=data_avgs, aes(x=sa_veg, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="sa plot")

p_avgs_veg <- ggplot(data=data_avgs, aes(x=sa_veg, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="Veg plot")



#All individual plots for avgs_pro
data_avgs_pro <- as.data.frame(top_df[[1]][[1]][[2]])

p_avgs_pro_sex <- ggplot(data=data_avgs_pro, aes(x=BA_sex, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  labs(title="Sex plot")

p_avgs_pro_BA <- ggplot(data=data_avgs_pro, aes(x=BA_sex, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="BA plot")

p_avgs_pro_sa <- ggplot(data=data_avgs_pro, aes(x=sa_veg, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="sa plot")

p_avgs_pro_veg <- ggplot(data=data_avgs_pro, aes(x=sa_veg, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="Veg plot")


#All individual plots for reps
data_reps <- as.data.frame(top_df[[1]][[2]])
p_reps_sex <- ggplot(data=data_reps, aes(x=BA_sex, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  labs(title="Sex plot")

p_reps_BA <- ggplot(data=data_reps, aes(x=BA_sex, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="BA plot")

p_reps_sa <- ggplot(data=data_reps, aes(x=sa_veg, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="sa plot")

p_reps_veg <- ggplot(data=data_reps, aes(x=sa_veg, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="Veg plot")



#All individual plots for reps_pro
data_reps_pro <- as.data.frame(top_df[[2]])
p_reps_pro_sex <- ggplot(data=data_reps_pro, aes(x=BA_sex, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  labs(title="Sex plot")

p_reps_pro_BA <- ggplot(data=data_reps_pro, aes(x=BA_sex, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="BA plot")

p_reps_pro_sa <- ggplot(data=data_reps_pro, aes(x=sa_veg, y=sa_sex)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="sa plot")

p_reps_pro_veg <- ggplot(data=data_reps_pro, aes(x=sa_veg, y=BA_veg)) +
  geom_violin() +
  geom_smooth(method=lm)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title="Veg plot")




library(cowplot)
 #Doing all of the grid plots
p_row <- plot_grid(p_avgs_sex, p_avgs_BA, p_avgs_sa, p_avgs_veg, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste("Correlation of correlation between methylation and gene expression for all genes where expression is averaged"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_grid_avgs <- plot_grid(
  title, p_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

p_row <- plot_grid(p_avgs_pro_sex, p_avgs_pro_BA, p_avgs_pro_sa, p_avgs_pro_veg, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste("Correlation of correlation between methylation and gene expression for all genes where expression is averaged and the promotor region is included"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_grid_avgs_pro <- plot_grid(
  title, p_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)


p_row <- plot_grid(p_reps_sex, p_reps_BA, p_reps_sa, p_reps_veg, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste("Correlation of correlation between methylation and gene expression for all genes where all replicates of the gene expression is included"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_grid_reps <- plot_grid(
  title, p_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)



p_row <- plot_grid(p_reps_pro_sex, p_reps_pro_BA, p_reps_pro_sa, p_reps_pro_veg, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    paste("Correlation of correlation between methylation and gene expression for all genes where all replicates of the gene expression and the promotor region is included"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_grid_reps_pro <- plot_grid(
  title, p_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)