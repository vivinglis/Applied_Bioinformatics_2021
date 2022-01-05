install.packages("dplyr")
install.packages("stringr")
library("stringr")
library("dplyr")

#Reading gff-file (Path needs to be correct)
gff <- read.delim("../annotation/genomeAnnotations/Neurospora_crassa.NC12.gff", TRUE, comment.char = "#")

#Filtering the gff-file
gff_rel_col <- gff[, c(1, 3, 4, 5, 9)]
gff_rel_rows <- gff_rel_col %>% filter(V3 == "gene")
gff_rel_rows$geneID = str_extract(gff_rel_rows$V9,"(NCU[0-9][0-9][0-9][0-9][0-9])")
gff_rel_rows$orthoID = str_extract(gff_rel_rows$V9,"(Ntsc[0-9]*)")
translation_table <- gff_rel_rows[, c(1, 3, 4, 6, 7)]

#Adding new column for startposition with promotor region
translation_table$StartWithPR = translation_table$V4-100

#Renaming column names
translation_table <- rename(translation_table, Chromosome = V1)
translation_table <- rename(translation_table, start = V4)
translation_table <- rename(translation_table, stop = V5)

#Write to file
write.table(translation_table, file='1_translation_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)
