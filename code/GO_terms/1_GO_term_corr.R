

#met <- read.table("1_orthologous_genes_methylation_with_promotor_region.tsv", TRUE)
met <- read.table("1_orthologous_genes_methylation.tsv", TRUE)
gene_expr <-read.table("../Applied_Bioinformatics_2021/data_tables/1_rlog_normalised_allsamples.tsv", TRUE)
#GO_terms <- read.delim("1_genes_per_GO_term.tsv")
GO_terms  <- read.delim("1_filtered_GO_terms_per_gene_4_or_more_genes_per_GO.tsv")

# Pearson correlation data table
pearson_corr <-
  setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Geneid", "Sampletype", "Value"))

# Set Geneid as index
rownames(met) <- met$Geneid
rownames(gene_expr) <- gene_expr$Geneid

# each is a new plot
sample_BA_veg <-
  grep("_BA_P[1-3]_veg", colnames(gene_expr), value = TRUE)
sample_BA_sex <-
  grep("_BA_P[1-3]_sex", colnames(gene_expr), value = TRUE)
sample_sa_veg <-
  grep("_sa_P[1-3]_veg", colnames(gene_expr), value = TRUE)
sample_sa_sex <-
  grep("_sa_P[1-3]_sex", colnames(gene_expr), value = TRUE)

samples <-
  list(sample_BA_veg, sample_BA_sex, sample_sa_veg, sample_sa_sex)
samples_name <-
  list("sample_BA_veg",
       "sample_BA_sex",
       "sample_sa_veg",
       "sample_sa_sex")

genes = rownames(met)
# genes <- c("Ntsc1042", "Ntsc1043")  # For test
for (i in 1:ncol(GO_terms)) {
  num = 1
  for (sample in samples) {
    per_gene <-
      setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Sample", "Met", "Expr"))
    
    # Loop over all elements in sample
    for (k in sample) {
      for (j in 1:nrow(GO_terms)) {
        if (is.na(GO_terms[j, i])) {
          break
        }
        # Removes P1-3 part if sample name to find methylation match
        met_match <- gsub("_P[1-3]", "", k)
        new <-
          data.frame(k, met[GO_terms[j, i], met_match], gene_expr[GO_terms[j, i], k])
        names(new) <- c("Sample", "Met", "Expr")
        per_gene <- rbind(per_gene, new)
      }
    }
    
    corr_value <-
      cor(per_gene$Met, per_gene$Expr) # Can also use cor.test()
    corr_new <-
      data.frame(colnames(GO_terms)[i], samples_name[num], corr_value)
    names(corr_new) <- c("GO_term", "Sampletype", "Value")
    pearson_corr <- rbind(pearson_corr, corr_new)
    num = num + 1
    
  }
}
stripchart(
  pearson_corr$Value,
  main = "Pearson correlation number",
  xlab = "",
  ylab = "Correlation value",
  method = "jitter",
  col = "blue",
  pch = 1,
  vertical = TRUE
)

#Write to file
write.table(pearson_corr, file='1_pear_corr_GO_terms_filtered_all_replicates.tsv', quote=FALSE, sep='\t', row.names = FALSE)
