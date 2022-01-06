
GO_trans <- read.delim("1_GO_translation_table.tsv", TRUE)
GO_trans <- GO_trans[order(GO_trans$Geneid),]

GO_terms <- unique(GO_trans$GO.term)

`%!in%` <- Negate(`%in%`)
for (GO_term in GO_terms){
  GO_col <- data.frame(matrix(NA_character_, ncol = 1, nrow = length(GO_trans$Geneid)))
  count = 1
  for (i in 1:length(GO_trans$Geneid)){
    if(GO_term == GO_trans$GO.term[i] & GO_trans$Geneid[i] %!in% GO_col$matrix.NA_character_..ncol...1..nrow...length.GO_trans.Geneid..){
      GO_col$matrix.NA_character_..ncol...1..nrow...length.GO_trans.Geneid..[count] = GO_trans$Geneid[i]
      count = count + 1
    }
  }
    if(GO_term == GO_terms[1]){
      GO_df <- GO_col
      colnames(GO_df)[ncol(GO_df)] <- paste0(GO_term)
    } else {
      GO_df <- cbind(GO_df, GO_col)
      colnames(GO_df)[ncol(GO_df)] <- paste0(GO_term)
    }
}

GO_df <- GO_df %>% distinct()

#Write to file
#write.table(GO_df, file='1_genes_per_GO_term.tsv', quote=FALSE, sep='\t', row.names = FALSE)
