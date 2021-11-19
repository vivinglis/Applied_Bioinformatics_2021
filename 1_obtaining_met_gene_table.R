library("dplyr")

#Script to obtain a table of the methylation level of all orthologous genes
#It computes the methylation level for a gene as number of methylated C-positions 
#divided by total number of C-positions in the gene. 


position_to_gene_vec <- function(met_df, tt){
  met_df <- met_df[order(met_df$Location),] 
  tt <- tt[order(tt$StartWithPR),]
  zeros = c()
  for (i in 1:(length(tt$orthoID))){
    zeros <- append(zeros, 0)
  }
  
  met_genes <- data.frame(zeros)
  loc_index = 0
  for (i in 1:(length(tt$orthoID))) {
    c_count = 0
    met_count = 0
    for (j in 1:(length(met_df$Location))){
      if(loc_index+j > length(met_df$Location)){
        break
      }
      if(met_df$Location[loc_index + j] > tt$stop[i]) {
        loc_index = loc_index + j
        break
      }
      if(met_df$Location[loc_index + j] > tt$StartWithPR[i] & met_df$Location[loc_index + j] < tt$stop[i]){
        c_count = c_count + 1

        if (met_df$methylated[loc_index + j] == 1){
          met_count = met_count +1
        }
      }
    }
    met_genes$zeros[i] = (met_count/c_count)
  }
  return(met_genes)
}

translation_table <- read.delim("1_translation_table.tsv", TRUE)
tt_filtered <- translation_table %>% filter(orthoID != is.na(orthoID))

samples = c("crassa_BA_sex",	"crassa_BA_veg",	"crassa_sa_sex",	"crassa_sa_veg",	"L1_BA_sex",	"L1_BA_veg",	"L1_sa_sex",	"L1_sa_veg", "L10_BA_sex",	"L10_BA_veg",	"L10_sa_sex",	"L10_sa_veg",	"L6_BA_sex",	"L6_BA_veg",	"L6_sa_sex",	"L6_sa_veg",	"sito_BA_sex",	"sito_BA_veg",	"sito_sa_sex",	"sito_sa_veg")
chroms <- c("I", "II", "III", "IV", "V", "VI", "VII")

for(chrom in chroms){
  tt<- tt_filtered %>% filter(Chromosome == chrom)
  df <- data.frame(tt$orthoID)
  for(sam in samples){
    met_sample <- met %>% filter(sample == sam)
    met_sample_chr <- met_sample %>% filter(Chr == chrom)
    df_sample <- position_to_gene_vec(met_sample_chr, tt)
    df <- cbind(df, df_sample)
    colnames(df)[ncol(df)] <- paste0(sam)
  }
    if(chrom == "I"){
      res <- df
    } else {
      res <- rbind(res, df)
    }
}
colnames(res)[1] <- "Geneid"

#Write to file
write.table(res, file='1_orthologous_genes_methylation_with_promotor_region.tsv', quote=FALSE, sep='\t', row.names = FALSE)
