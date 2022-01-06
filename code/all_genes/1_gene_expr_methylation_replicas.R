# Keeping replicates separate 

# Methylation and gene expression table
gene_expr <- read.table("Applied_Bioinformatics_2021/1_rlog_normalised_allsamples.tsv", header = TRUE) 
# met <- read.table("Applied_Bioinformatics_2021/1_orthologous_genes_methylation.tsv", header = TRUE)

# Methylation with promotor region 
met <- read.table("Applied_Bioinformatics_2021/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)

# Pearson correlation data table 
pearson_corr <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Geneid", "Sampletype", "Value"))

# Set Geneid as index
rownames(met) <- met$Geneid
rownames(gene_expr) <- gene_expr$Geneid

# each is a new plot 
sample_BA_veg <- grep("_BA_P[1-3]_veg", colnames(gene_expr), value = TRUE)
sample_BA_sex <- grep("_BA_P[1-3]_sex", colnames(gene_expr), value = TRUE)
sample_sa_veg <- grep("_sa_P[1-3]_veg", colnames(gene_expr), value = TRUE)
sample_sa_sex <- grep("_sa_P[1-3]_sex", colnames(gene_expr), value = TRUE)

samples <- list(sample_BA_veg, sample_BA_sex, sample_sa_veg, sample_sa_sex)
samples_name <- list("sample_BA_veg", "sample_BA_sex", "sample_sa_veg", "sample_sa_sex")

# genes = rownames(met)
genes <- c("Ntsc1042", "Ntsc1043")  # For test 

for (gene in genes) {
  num = 1
  for (sample in samples){
    per_gene <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Sample", "Met", "Expr")) 
    
    # Loop over all elements in sample 
    for(i in sample) {
      # Removes P1-3 part if sample name to find methylation match
      met_match <- gsub("_P[1-3]", "", i )
      new <- data.frame(i, met[gene, met_match], gene_expr[gene, i])
      names(new) <- c("Sample", "Met", "Expr")
      per_gene <- rbind(per_gene, new)
    }
    rownames(per_gene) <- per_gene$Sample
    # Plot # if time use ggplot2 instead since they look better
    #plot(per_gene$Met, per_gene$Expr, ann=FALSE)
    #title(xlab="Methylation")
    #title(ylab="Gene expression")                                
    
    corr_value <- cor(per_gene$Met, per_gene$Expr) # Can also use cor.test() 
    corr_new <- data.frame(gene, samples_name[num], corr_value)
    names(corr_new) <- c("Geneid", "Sampletype", "Value")
    pearson_corr <- rbind(pearson_corr, corr_new)
    num = num +1
  }
}
stripchart(pearson_corr$Value,
           main="Pearson correlation number",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col="blue",
           pch=1, 
           vertical=TRUE
)
