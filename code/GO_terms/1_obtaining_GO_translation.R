library("dplyr")

translation_table <- read.delim("1_translation_table.tsv", TRUE)
GO_terms <- read.delim("../annotation/GO/go_for_nc12.tsv", TRUE)

GO_terms <- GO_terms[order(GO_terms$Locus),] 
tt <- translation_table[order(translation_table$geneID),]

GO_terms$Locus <- gsub(".7","",GO_terms$Locus)

GO_translation <- cbind(tt$orthoID[FALSE],GO_terms[FALSE,])
colnames(GO_translation)[1] <- "Geneid"

loc_index = 0
for (i in 1:(length(tt$geneID))) {
  for (j in 1:(length(GO_terms$Locus))){
    # if(loc_index+j > length(GO_terms$Locus)){
    #   print("INDEX TOO BIG")
    #   break
    # }
    # if(GO_terms$Locus[loc_index + j] > tt$geneID[i]) {
    #   loc_index = loc_index + j
    #   print("LOC_INDEX: ")
    #   print(loc_index)
    #   break
    # }
    if(GO_terms$Locus[loc_index + j] == tt$geneID[i]){
      
      GO_row <- cbind(tt$orthoID[i],GO_terms[loc_index + j ,])
      GO_translation[nrow(GO_translation) + 1,] = GO_row
      #print(GO_translation)
    }
  }
}

GO_translation_filtered <- GO_translation %>% filter(Geneid != is.na(Geneid))

#Write to file
write.table(GO_translation_filtered, file='1_GO_translation_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)