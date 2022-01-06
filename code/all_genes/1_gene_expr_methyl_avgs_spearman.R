# Using avg gene expression values, 
# keep one of line 5 or 8 commented depending on want with or without promotor region  

# Methylation and gene expression table
gene_expr <- read.table("Applied_Bioinformatics_2021/data_tables/1_rlog_normalised_avgs.tsv", header = TRUE) 
met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation.tsv", header = TRUE)

# Methylation with promotor region 
#met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)

# Spearman correlation data table 
spearman_corr <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Geneid", "Allel", "Tissue", "Value", "k_value"))

# Set Geneid as index
rownames(met) <- met$Geneid
rownames(gene_expr) <- gene_expr$Geneid

# each is a new plot 
sample_BA_veg <- grep("_BA_veg", colnames(gene_expr), value = TRUE)
sample_BA_sex <- grep("_BA_sex", colnames(gene_expr), value = TRUE)
sample_sa_veg <- grep("_sa_veg", colnames(gene_expr), value = TRUE)
sample_sa_sex <- grep("_sa_sex", colnames(gene_expr), value = TRUE)

samples <- list(sample_BA_veg, sample_BA_sex, sample_sa_veg, sample_sa_sex)
samples_allel <- list("BA", "BA", "sa", "sa")
samples_tissue <- list("veg", "sex", "veg", "sex")

genes = rownames(met)
# genes <- c("Ntsc1042", "Ntsc1043")  # For test of code

for (gene in genes) {
  num = 1
  for (sample in samples){
    per_gene <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Sample", "Met", "Expr")) 
    
    # Loop over all elements in sample 
    for(i in sample) {
      new <- data.frame(i, met[gene, i], gene_expr[gene, i])
      names(new) <- c("Sample", "Met", "Expr")
      per_gene <- rbind(per_gene, new)
    }
    #rownames(per_gene) <- per_gene$Sample
    
    corr_value <- cor(per_gene$Met, per_gene$Expr, method = "spearman") # Can also use cor.test()
    if (any(is.na(per_gene$Met))){
      k_value = NA
    }
    else {
      linearMod <- lm(Expr ~ Met, data=per_gene)
      k_value <- unname(linearMod$coefficients[2])
    }
    corr_new <- data.frame(gene, samples_allel[num], samples_tissue[num], corr_value, k_value)
    names(corr_new) <- c("Geneid", "Allel", "Tissue", "Value", "k_value")
    spearman_corr <- rbind(spearman_corr, corr_new)
    num = num +1
  }
}
write.table(spearman_corr, file = "1_spearman_corr_avgs.tsv")
#write.table(spearman_corr, file = "1_spearman_corr_avgs_with_pro.tsv")

# Not separated on Ba sa ------------------------------------------------------------------

# Methylation and gene expression table
gene_expr <- read.table("Applied_Bioinformatics_2021/1_rlog_normalised_avgs.tsv", header = TRUE) 
met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation.tsv", header = TRUE)

# Methylation with promotor region 
#met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)

# Spearman correlation data table 
spearman_corr <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Geneid", "Tissue", "Value", "k_value"))

# Set Geneid as index
rownames(met) <- met$Geneid
rownames(gene_expr) <- gene_expr$Geneid

# each is a new plot 
sample_veg <- grep("_veg", colnames(gene_expr), value = TRUE)
sample_sex <- grep("_sex", colnames(gene_expr), value = TRUE)

samples <- list(sample_veg, sample_sex)
samples_tissue <- list("veg", "sex")


genes = rownames(met)
# genes <- c("Ntsc1042", "Ntsc1043")  # For test 

for (gene in genes) {
  num = 1
  for (sample in samples){
    per_gene <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Sample", "Met", "Expr")) 
    
    # Loop over all elements in sample 
    for(i in sample) {
      new <- data.frame(i, met[gene, i], gene_expr[gene, i])
      names(new) <- c("Sample", "Met", "Expr")
      per_gene <- rbind(per_gene, new)
    }
    #rownames(per_gene) <- per_gene$Sample
    
    corr_value <- cor(per_gene$Met, per_gene$Expr, method = "spearman") # Can also use cor.test()
    if (any(is.na(per_gene$Met))){
      k_value = NA
    }
    else {
      linearMod <- lm(Expr ~ Met, data=per_gene)
      k_value <- unname(linearMod$coefficients[2])
    }
    corr_new <- data.frame(gene, samples_tissue[num], corr_value, k_value)
    names(corr_new) <- c("Geneid", "Tissue", "Value", "k_value")
    spearman_corr <- rbind(spearman_corr, corr_new)
    num = num +1
  }
}
write.table(spearman_corr, file = "1_spearman_corr_avgs_non_allel_sep.tsv")
#write.table(spearman_corr, file = "1_spearman_corr_avgs_with_pro_non_allel_sep.tsv")
