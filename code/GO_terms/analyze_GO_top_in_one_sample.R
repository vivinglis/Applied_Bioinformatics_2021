#Correlation analysis
#Script to plot top [top_num] negatively correlated GO term from one sample (sample_name) and their correlation value in the other samples.

# Choose [top_num](1-133) and [sample_name](sa_sex, BA_sex, BA_veg, sa_veg) 
top_num <- 2
sample_name <- "sa_sex"


library("stringi")
library("ggplot2")
library("tidyverse")
# Read in tables
pearson_corr_avgs <- read.table("../../data_tables/GO_terms/1_data_for_corr_avg_and_std_GO_terms_filtered_4_avgs.tsv", header = TRUE) 
pearson_corr_avgs_with_pro <- read.table("../../data_tables/GO_terms/1_data_for_corr_avg_and_std_GO_terms_filtered_4_avgs_pro.tsv", header = TRUE) 
pearson_corr_reps <- read.table("../../data_tables/GO_terms/1_data_for_corr_avg_and_std_GO_terms_filtered_4_reps.tsv", header = TRUE) 
pearson_corr_reps_with_pro <- read.table("../../data_tables/GO_terms/1_data_for_corr_avg_and_std_GO_terms_filtered_4_reps_pro.tsv", header = TRUE) 


# top_pearson_corr_avgs <- head(pearson_corr_avgs[order(pearson_corr_avgs$Value),], top_num)
# top_pearson_corr_avgs_with_pro <- head(pearson_corr_avgs_with_pro[order(pearson_corr_avgs_with_pro$Value),], top_num)
# top_pearson_corr_rep <- head(pearson_corr_rep[order(pearson_corr_rep$Value),], top_num)
# top_pearson_corr_rep_with_pro <- head(pearson_corr_rep_with_pro[order(pearson_corr_rep_with_pro$Value),], top_num)

corrs <- list(pearson_corr_avgs, pearson_corr_avgs_with_pro, pearson_corr_reps, pearson_corr_reps_with_pro)
avg_df <- as.data.frame(matrix(nrow=0, ncol = 7))
top_df <- list()
for (corr in corrs){
  # Dividing data into subsets
  pearson_corr_sa_veg <- corr %>% filter(Sample == "sa_veg")
  pearson_corr_BA_sex <- corr %>% filter(Sample == "BA_sex")
  pearson_corr_BA_veg <- corr %>% filter(Sample == "BA_veg")
  pearson_corr_sa_sex <- corr %>% filter(Sample == "sa_sex")
  
  rownames(pearson_corr_BA_veg) <- pearson_corr_BA_veg$GO_term
  rownames(pearson_corr_sa_veg) <- pearson_corr_sa_veg$GO_term
  rownames(pearson_corr_sa_sex) <- pearson_corr_sa_sex$GO_term
  rownames(pearson_corr_BA_sex) <- pearson_corr_BA_sex$GO_term
  
  pearson_corr_sex <- rbind(pearson_corr_BA_sex, pearson_corr_sa_sex)
  pearson_corr_veg <- rbind(pearson_corr_BA_veg, pearson_corr_sa_veg)
  
  
  top_pearson_corr_sa_veg <- head(pearson_corr_sa_veg[order(pearson_corr_sa_veg$Avg_corr),], top_num)
  top_pearson_corr_BA_veg <- head(pearson_corr_BA_veg[order(pearson_corr_BA_veg$Avg_corr),], top_num)
  top_pearson_corr_sa_sex <- head(pearson_corr_sa_sex[order(pearson_corr_sa_sex$Avg_corr),], top_num)
  top_pearson_corr_BA_sex <- head(pearson_corr_BA_sex[order(pearson_corr_BA_sex$Avg_corr),], top_num)
  
  top_pearson_corr_sex <- head(pearson_corr_sex[order(pearson_corr_sex$Avg_corr),], top_num)
  top_pearson_corr_veg <- head(pearson_corr_veg[order(pearson_corr_veg$Avg_corr),], top_num)
  
  
  #Choose one of these for top 100 of that sample 
  
  # BA_sex
  if(sample_name == "BA_sex"){
    top_veg <- data.frame(GO_term=character(), BA_veg = double(), sa_veg = double(),sa_sex = double())
    for (top_corr in top_pearson_corr_BA_sex$GO_term){
      corr_veg <- data.frame(GO_term = top_corr, BA_veg = pearson_corr_BA_veg[top_corr, "Avg_corr"], sa_veg = pearson_corr_sa_veg[top_corr, "Avg_corr"], sa_sex = pearson_corr_sa_sex[top_corr, "Avg_corr"])
      top_veg <- rbind(top_veg, corr_veg)
    }
    a <- list(top_pearson_corr_BA_sex$Avg_corr, top_veg$sa_sex, top_veg$BA_veg, top_veg$sa_veg, top_veg$GO_term ) #BA_sex
  }
  

  # sa_sex
  if (sample_name == "sa_sex"){
    top_veg <- data.frame(GO_term=character(), BA_veg = double(), sa_veg = double(),BA_sex = double())
    for (top_corr in top_pearson_corr_sa_sex$GO_term){
      corr_veg <- data.frame(GO_term = top_corr, BA_veg = pearson_corr_BA_veg[top_corr, "Avg_corr"], sa_veg = pearson_corr_sa_veg[top_corr, "Avg_corr"], BA_sex = pearson_corr_BA_sex[top_corr, "Avg_corr"])
      top_veg <- rbind(top_veg, corr_veg)
    }
    a <- list(top_veg$BA_sex, top_pearson_corr_sa_sex$Avg_corr, top_veg$BA_veg, top_veg$sa_veg, top_veg$GO_term ) #sa_sex
  }
  


  # BA_veg
  if (sample_name == "BA_veg"){
    top_veg <- data.frame(GO_term=character(), sa_sex = double(), sa_veg = double(),BA_sex = double())
    for (top_corr in top_pearson_corr_BA_veg$GO_term){
      corr_veg <- data.frame(GO_term = top_corr, sa_sex = pearson_corr_sa_sex[top_corr, "Avg_corr"], sa_veg = pearson_corr_sa_veg[top_corr, "Avg_corr"], BA_sex = pearson_corr_BA_sex[top_corr, "Avg_corr"])
      top_veg <- rbind(top_veg, corr_veg)
    }
  
    a <- list(top_veg$BA_sex, top_veg$sa_sex, top_pearson_corr_BA_veg$Avg_corr, top_veg$sa_veg, top_veg$GO_term ) #BA_veg
  }
  
  # sa_veg
  if (sample_name == "sa_veg"){
    top_veg <- data.frame(GO_term=character(), sa_sex = double(), BA_veg = double(),BA_sex = double())
    for (top_corr in top_pearson_corr_sa_veg$GO_term){
      corr_veg <- data.frame(GO_term = top_corr, sa_sex = pearson_corr_sa_sex[top_corr, "Avg_corr"], BA_veg = pearson_corr_BA_veg[top_corr, "Avg_corr"], BA_sex = pearson_corr_BA_sex[top_corr, "Avg_corr"])
      top_veg <- rbind(top_veg, corr_veg)
    }
    a <- list(top_veg$BA_sex, top_veg$sa_sex, top_veg$BA_veg, top_pearson_corr_sa_veg$Avg_corr, top_veg$GO_term) #sa_veg
  }
  top_df <- list(top_df, a)
}

BA_sex = c()
for (i in 1:(top_num)){
  BA_sex <- append(BA_sex, "BA_sex")
}
sa_sex = c()
for (i in 1:(top_num)){
  sa_sex <- append(sa_sex, "sa_sex")
}
BA_veg = c()
for (i in 1:(top_num)){
  BA_veg <- append(BA_veg, "BA_veg")
}
sa_veg = c()
for (i in 1:(top_num)){
  sa_veg <- append(sa_veg, "sa_veg")
}


top_data_avgs <- as.data.frame(top_df[[1]][[1]][[1]][[2]])
colnames(top_data_avgs) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg", "GO_term")
BA_sex_df <- cbind(BA_sex, top_data_avgs$BA_sex, top_data_avgs$GO_term)
sa_sex_df <- cbind(sa_sex, top_data_avgs$sa_sex, top_data_avgs$GO_term)
BA_veg_df <- cbind(BA_veg, top_data_avgs$BA_veg, top_data_avgs$GO_term)
sa_veg_df <- cbind(sa_veg, top_data_avgs$sa_veg, top_data_avgs$GO_term)
plot_data_avgs <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_avgs) = c("sample", "correlation", "GO_term")
plot_data_avgs$sample <- as.factor(plot_data_avgs$sample)
plot_data_avgs$correlation <- as.double(plot_data_avgs$correlation)

p_avgs <- ggplot(data=plot_data_avgs, aes(x=sample, y=correlation, fill=sample, label=GO_term)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged")
if (top_num < 4){
  p_avgs <- p_avgs + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_text(size=3) +
    geom_label() 
}

top_data_avgs_pro <- as.data.frame(top_df[[1]][[1]][[2]])
colnames(top_data_avgs_pro) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg", "GO_term")
BA_sex_df <- cbind(BA_sex, top_data_avgs_pro$BA_sex, top_data_avgs_pro$GO_term)
sa_sex_df <- cbind(sa_sex, top_data_avgs_pro$sa_sex, top_data_avgs_pro$GO_term)
BA_veg_df <- cbind(BA_veg, top_data_avgs_pro$BA_veg, top_data_avgs_pro$GO_term)
sa_veg_df <- cbind(sa_veg, top_data_avgs_pro$sa_veg, top_data_avgs_pro$GO_term)
plot_data_avgs_pro <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_avgs_pro) = c("sample", "correlation", "GO_term")
plot_data_avgs_pro$sample <- as.factor(plot_data_avgs_pro$sample)
plot_data_avgs_pro$correlation <- as.double(plot_data_avgs_pro$correlation)

p_avgs_pro <- ggplot(data=plot_data_avgs_pro, aes(x=sample, y=correlation, fill=sample, label=GO_term)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, fill="white")+ 
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged with promotor region")
if (top_num < 4){
  p_avgs_pro <- p_avgs_pro + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_text(size=3) +
    geom_label() 
}

top_data_reps <- as.data.frame(top_df[[1]][[2]])
colnames(top_data_reps) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg", "GO_term")

BA_sex_df <- cbind(BA_sex, top_data_reps$BA_sex, top_data_reps$GO_term)
sa_sex_df <- cbind(sa_sex, top_data_reps$sa_sex, top_data_reps$GO_term)
BA_veg_df <- cbind(BA_veg, top_data_reps$BA_veg, top_data_reps$GO_term)
sa_veg_df <- cbind(sa_veg, top_data_reps$sa_veg, top_data_reps$GO_term)
plot_data_reps <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_reps) = c("sample", "correlation", "GO_term")
plot_data_reps$sample <- as.factor(plot_data_reps$sample)
plot_data_reps$correlation <- as.double(plot_data_reps$correlation)

p_reps <- ggplot(data=plot_data_reps, aes(x=sample, y=correlation, fill=sample, label=GO_term)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, fill="white")+ 
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg"))+
  labs(title="All gene expression replicates")
if (top_num < 4){
  p_reps <- p_reps + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_text(size=3) +
    geom_label() 
}

top_data_reps_pro <- as.data.frame(top_df[[2]])
colnames(top_data_reps_pro) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg", "GO_term")

BA_sex_df <- cbind(BA_sex, top_data_reps_pro$BA_sex, top_data_reps_pro$GO_term)
sa_sex_df <- cbind(sa_sex, top_data_reps_pro$sa_sex, top_data_reps_pro$GO_term)
BA_veg_df <- cbind(BA_veg, top_data_reps_pro$BA_veg, top_data_reps_pro$GO_term)
sa_veg_df <- cbind(sa_veg, top_data_reps_pro$sa_veg, top_data_reps_pro$GO_term)
plot_data_reps_pro <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_reps_pro) = c("sample", "correlation", "GO_term")
plot_data_reps_pro$sample <- as.factor(plot_data_reps_pro$sample)
plot_data_reps_pro$correlation <- as.double(plot_data_reps_pro$correlation)

p_reps_pro <- ggplot(data=plot_data_reps_pro, aes(x=sample, y=correlation, fill=sample, label=GO_term)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, fill="white")+ 
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="All gene expression replicates with promotor region")
if (top_num < 4){
  p_reps_pro <- p_reps_pro + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    geom_text(size=3) +
    geom_label() 
}

library(cowplot)

p_row <- plot_grid(p_avgs, p_avgs_pro, p_reps, p_reps_pro, labels = "AUTO")

title <- ggdraw() + 
  draw_label(
    paste0("Average correlation values for top ",top_num," correlated GO terms for ",sample_name),
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


