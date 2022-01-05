
library("stringi")

pearson_corr_avgs <- read.table("1_pear_corr_GO_terms_avgs.tsv", header = TRUE) 
GO_df <- read.delim("1_genes_per_GO_term.tsv")

# 
# ulst <- lapply(GO_df, unique)
# k <- lengths(ulst)-1
# dis <- table(k)
# dis_df <- as.data.frame(dis)
# gene_GO_dis <- ggplot(dis_df, aes(x=k, y=Freq)) +
#   geom_bar(stat="identity", fill="steelblue")+
#   geom_text(aes(label=Freq), vjust=-0.3, size=3.5)+
#   labs(title="Distribution of number of genes in each GO:term", 
#        x="Number of genes in a GO:term", y = "Frequency of the GO:term")+
#   theme_minimal()
#pearson_corr <- pearson_corr_avgs[order(pearson_corr_avgs$Value),]

# len = 0
# i = 1
# c = 0
# while (pearson_corr$Value[i] < 0){
#   a <- pearson_corr[i,1]
#   b <- gsub("[.]", " ", a)
#   len <- length(unique(GO_df[,a]))
#   #print(i)
#   i <- i + 1
#   if (len > 10){
#     print(b)
#     print(len)
#     print(pearson_corr$Value[i])
#     c = c +1 
#     print(c)
#   }
#  
#  }
GO_filt <- as.data.frame(matrix(nrow=165, ncol=0))
for (i in 1:(ncol(GO_df))) {
  if (length(unique(GO_df[,i])) > 10) {
    GO_filt <- cbind(GO_filt, GO_df[,i])
    colnames(GO_filt)[ncol(GO_filt)] = colnames(GO_df)[i]
  }
}

write.table(GO_filt, "1_filtered_GO_terms_per_gene_10_or_more_genes_per_GO.tsv", quote=FALSE, sep='\t', row.names = FALSE)

# 
# #Correlation analysis
# library("stringi")
# 
# # Read in tables
pearson_corr_avgs <- read.table('1_pear_corr_GO_terms_filtered_10_avgs.tsv', header = TRUE) 
pearson_corr_avgs_with_pro <- read.table("1_pear_corr_GO_terms_avgs_with_promotor_region.tsv", header = TRUE)
pearson_corr_rep <- read.table("1_pear_corr_GO_terms_all_replicates.tsv", header = TRUE)
pearson_corr_rep_with_pro <- read.table("1_pear_corr_GO_terms_all_replicates_with_promotor_region.tsv", header = TRUE)


plot_avgs <- ggplot(data=pearson_corr_avgs, aes(x=Sampletype, y=Value, fill=Sampletype)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("sample_BA_sex", "sample_sa_sex", "sample_BA_veg", "sample_sa_veg")) +
  labs(title="Gene expression values averaged", x="Sample", y="Correlation")

plot_avgs_pro <- ggplot(data=pearson_corr_avgs_with_pro, aes(x=Sampletype, y=Value, fill=Sampletype)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("sample_BA_sex", "sample_sa_sex", "sample_BA_veg", "sample_sa_veg")) +
  labs(title="Gene expression values averaged with promotor region", x="Sample", y="Correlation")

plot_reps <- ggplot(data=pearson_corr_rep, aes(x=Sampletype, y=Value, fill=Sampletype)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("sample_BA_sex", "sample_sa_sex", "sample_BA_veg", "sample_sa_veg")) +
  labs(title="All gene expression replicates", x="Sample", y="Correlation")

plot_reps_pro <- ggplot(data=pearson_corr_rep_with_pro, aes(x=Sampletype, y=Value, fill=Sampletype)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("sample_BA_sex", "sample_sa_sex", "sample_BA_veg", "sample_sa_veg")) +
  labs(title="All gene expression replicates with promotor region", x="Sample", y="Correlation")

library(cowplot)

plot_row <- plot_grid(plot_avgs, plot_avgs_pro, plot_reps, plot_reps_pro, labels = "AUTO")

tit <- ggdraw() +
  draw_label(
    "Correlation values between methylation and gene expression for GO terms for each sample",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 17)
  )
p_all <- plot_grid(
  tit, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

# Take top number values from each dataframe
top_num <- 100
# top_pearson_corr_avgs <- head(pearson_corr_avgs[order(pearson_corr_avgs$Value),], top_num)
# top_pearson_corr_avgs_with_pro <- head(pearson_corr_avgs_with_pro[order(pearson_corr_avgs_with_pro$Value),], top_num)
# top_pearson_corr_rep <- head(pearson_corr_rep[order(pearson_corr_rep$Value),], top_num)
# top_pearson_corr_rep_with_pro <- head(pearson_corr_rep_with_pro[order(pearson_corr_rep_with_pro$Value),], top_num)

corrs <- list(pearson_corr_avgs, pearson_corr_avgs_with_pro, pearson_corr_rep, pearson_corr_rep_with_pro)
avg_df <- as.data.frame(matrix(nrow=0, ncol = 7))
top_df <- list()
for (corr in corrs){
  # Dividing data into subsets
  pearson_corr_sa_veg <- corr %>% filter(Sampletype == "sample_sa_veg")
  pearson_corr_BA_sex <- corr %>% filter(Sampletype == "sample_BA_sex")
  pearson_corr_BA_veg <- corr %>% filter(Sampletype == "sample_BA_veg")
  pearson_corr_sa_sex <- corr %>% filter(Sampletype == "sample_sa_sex")

  rownames(pearson_corr_BA_veg) <- pearson_corr_BA_veg$GO_term
  rownames(pearson_corr_sa_veg) <- pearson_corr_sa_veg$GO_term
  rownames(pearson_corr_sa_sex) <- pearson_corr_sa_sex$GO_term
  rownames(pearson_corr_BA_sex) <- pearson_corr_BA_sex$GO_term

  pearson_corr_sex <- rbind(pearson_corr_BA_sex, pearson_corr_sa_sex)
  pearson_corr_veg <- rbind(pearson_corr_BA_veg, pearson_corr_sa_veg)


  top_pearson_corr_sa_veg <- head(pearson_corr_sa_veg[order(pearson_corr_sa_veg$Value),], top_num)
  top_pearson_corr_BA_veg <- head(pearson_corr_BA_veg[order(pearson_corr_BA_veg$Value),], top_num)
  top_pearson_corr_sa_sex <- head(pearson_corr_sa_sex[order(pearson_corr_sa_sex$Value),], top_num)
  top_pearson_corr_BA_sex <- head(pearson_corr_BA_sex[order(pearson_corr_BA_sex$Value),], top_num)

  top_pearson_corr_sex <- head(pearson_corr_sex[order(pearson_corr_sex$Value),], top_num)
  top_pearson_corr_veg <- head(pearson_corr_veg[order(pearson_corr_veg$Value),], top_num)

  #Mean correlation values
  avg_tot <- mean(corr$Value, na.rm=TRUE)
  avg_sex <- mean(pearson_corr_sex$Value, na.rm=TRUE)
  avg_veg <- mean(pearson_corr_veg$Value, na.rm=TRUE)
  avg_BA_sex <- mean(pearson_corr_BA_sex$Value, na.rm=TRUE)
  avg_sa_sex <- mean(pearson_corr_sa_sex$Value, na.rm=TRUE)
  avg_BA_veg <- mean(pearson_corr_BA_veg$Value, na.rm=TRUE)
  avg_sa_veg <- mean(pearson_corr_sa_veg$Value, na.rm=TRUE)
  df <- data.frame(tot=avg_tot, sex=avg_sex, veg=avg_veg, BA_sex=avg_BA_sex, sa_sex=avg_sa_sex, BA_veg=avg_BA_veg, sa_veg=avg_sa_veg)
  avg_df <- rbind(avg_df, df)
  print(avg_df)

  # Plotting
  x <- list(top_pearson_corr_BA_sex$Value, top_pearson_corr_sa_sex$Value, top_pearson_corr_BA_veg$Value, top_pearson_corr_sa_veg$Value)
  top_df <- list(top_df, x)
}


sample_BA_sex = c()
for (i in 1:(top_num)){
  sample_BA_sex <- append(sample_BA_sex, "BA_sex")
}
sample_sa_sex = c()
for (i in 1:(top_num)){
  sample_sa_sex <- append(sample_sa_sex, "sa_sex")
}
sample_BA_veg = c()
for (i in 1:(top_num)){
  sample_BA_veg <- append(sample_BA_veg, "BA_veg")
}
sample_sa_veg = c()
for (i in 1:(top_num)){
  sample_sa_veg <- append(sample_sa_veg, "sa_veg")
}


top_data_avgs <- as.data.frame(top_df[[1]][[1]][[1]][[2]])
colnames(top_data_avgs) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg")
BA_sex_df <- cbind(sample_BA_sex, top_data_avgs$BA_sex)
sa_sex_df <- cbind(sample_sa_sex, top_data_avgs$sa_sex)
BA_veg_df <- cbind(sample_BA_veg, top_data_avgs$BA_veg)
sa_veg_df <- cbind(sample_sa_veg, top_data_avgs$sa_veg)
plot_data_avgs <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_avgs) = c("sample", "correlation")
plot_data_avgs$sample <- as.factor(plot_data_avgs$sample)
plot_data_avgs$correlation <- as.double(plot_data_avgs$correlation)

p_avgs <- ggplot(data=plot_data_avgs, aes(x=sample, y=correlation, fill=sample)) +
  ylim(c(-1, 1)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged")


top_data_avgs_pro <- as.data.frame(top_df[[1]][[1]][[2]])
colnames(top_data_avgs_pro) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg")
BA_sex_df <- cbind(sample_BA_sex, top_data_avgs_pro$BA_sex)
sa_sex_df <- cbind(sample_sa_sex, top_data_avgs_pro$sa_sex)
BA_veg_df <- cbind(sample_BA_veg, top_data_avgs_pro$BA_veg)
sa_veg_df <- cbind(sample_sa_veg, top_data_avgs_pro$sa_veg)
plot_data_avgs_pro <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_avgs_pro) = c("sample", "correlation")
plot_data_avgs_pro$sample <- as.factor(plot_data_avgs_pro$sample)
plot_data_avgs_pro$correlation <- as.double(plot_data_avgs_pro$correlation)

p_avgs_pro <- ggplot(data=plot_data_avgs_pro, aes(x=sample, y=correlation, fill=sample)) +
  ylim(c(-1, -0.45)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="Gene expression values averaged with promotor region")


top_data_reps <- as.data.frame(top_df[[1]][[2]])
colnames(top_data_reps) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg")

BA_sex_df <- cbind(sample_BA_sex, top_data_reps$BA_sex)
sa_sex_df <- cbind(sample_sa_sex, top_data_reps$sa_sex)
BA_veg_df <- cbind(sample_BA_veg, top_data_reps$BA_veg)
sa_veg_df <- cbind(sample_sa_veg, top_data_reps$sa_veg)
plot_data_reps <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_reps) = c("sample", "correlation")
plot_data_reps$sample <- as.factor(plot_data_reps$sample)
plot_data_reps$correlation <- as.double(plot_data_reps$correlation)

p_reps <- ggplot(data=plot_data_reps, aes(x=sample, y=correlation, fill=sample)) +
  ylim(c(-1, -0.45)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg"))+
  labs(title="All gene expression replicates")

top_data_reps_pro <- as.data.frame(top_df[[2]])
colnames(top_data_reps_pro) = c("BA_sex", "sa_sex", "BA_veg", "sa_veg")

BA_sex_df <- cbind(sample_BA_sex, top_data_reps_pro$BA_sex)
sa_sex_df <- cbind(sample_sa_sex, top_data_reps_pro$sa_sex)
BA_veg_df <- cbind(sample_BA_veg, top_data_reps_pro$BA_veg)
sa_veg_df <- cbind(sample_sa_veg, top_data_reps_pro$sa_veg)
plot_data_reps_pro <- as.data.frame(rbind(BA_sex_df, sa_sex_df, BA_veg_df, sa_veg_df))
colnames(plot_data_reps_pro) = c("sample", "correlation")
plot_data_reps_pro$sample <- as.factor(plot_data_reps_pro$sample)
plot_data_reps_pro$correlation <- as.double(plot_data_reps_pro$correlation)

p_reps_pro <- ggplot(data=plot_data_reps_pro, aes(x=sample, y=correlation, fill=sample)) +
  ylim(c(-1, -0.45)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")+
  scale_x_discrete(limits=c("BA_sex", "sa_sex", "BA_veg", "sa_veg")) +
  labs(title="All gene expression replicates with promotor region")

library(cowplot)

p_row <- plot_grid(p_avgs, p_avgs_pro, p_reps, p_reps_pro, labels = "AUTO")

title <- ggdraw() +
  draw_label(
    "Correlation values between methylation and gene expression for top 100 correlated GO terms for each sample",
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

#par(mfrow=c(2,2))
# stripchart(top_df[[1]][[1]][[1]][[2]],
#            main="Pearson correlation value: avgs",
#            xlab="BA_sex   sa_sex    BA_veg    sa_veg",
#            ylab="Correlation value",
#            method="jitter",
#            col=c("green", "blue","red", "orange"),
#            pch=1,
#            vertical=TRUE,
#            ylim = c(-1, 1)
# )
# stripchart(top_df[[1]][[1]][[2]],
#            main="Pearson correlation value: avgs with pro",
#            xlab="BA_sex   sa_sex    BA_veg    sa_veg",
#            ylab="Correlation value",
#            method="jitter",
#            col=c("green", "blue","red", "orange"),
#            pch=1,
#            vertical=TRUE,
#            ylim = c(-1, 1)
# )
# stripchart(top_df[[1]][[2]],
#            main="Pearson correlation value: reps",
#            xlab="BA_sex   sa_sex    BA_veg    sa_veg",
#            ylab="Correlation value",
#            method="jitter",
#            col=c("green", "blue","red", "orange"),
#            pch=1,
#            vertical=TRUE,
#            ylim = c(-1, 1)
# )
# stripchart(top_df[[2]],
#            main="Pearson correlation value: reps with pro",
#            xlab="BA_sex   sa_sex    BA_veg    sa_veg",
#            ylab="Correlation value",
#            method="jitter",
#            col=c("green", "blue","red", "orange"),
#            pch=1,
#            vertical=TRUE,
#            ylim = c(-1, 1)
#)





