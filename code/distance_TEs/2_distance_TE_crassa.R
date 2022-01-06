# Distance calculations 
rm(list = ls())

# Crassa

# Packages
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

# Import datasets 
TEs <- read.delim("~/Applied_bioinformatics/NeurosporaMethylation/annotation/TE/crassaBA_2489_TEs", header=FALSE)
translation_table <- read.delim("~/Applied_bioinformatics/Applied_Bioinformatics_2021/data_tables/1_translation_table.tsv")
methylation_table <- read.delim("~/Applied_bioinformatics/Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation.tsv")

rownames(methylation_table) <- methylation_table$Geneid

# Add correct headers 
names(TEs) <- c("chromosome", "start", "stop", "perc_div", "matching_repeat", "repeat_class/family", "id")

gene_table <- filter(translation_table, is.na(orthoID) == FALSE) 
rownames(gene_table) <- gene_table$orthoID
# gene is orthoID in the rest of script  
genes <- gene_table$orthoID
distance_matrix <- setNames(data.frame(matrix(ncol = 3, nrow = length(genes))), c("distance", "met_veg", "met_sex"))
distance_matrix <- cbind(genes, distance_matrix)
names(distance_matrix)[1] <- "geneID"    #called othoID in translational table
rownames(distance_matrix) <- distance_matrix$geneID

# Chromosomes in common TEs and organism
chr_te <- unique(TEs[c("chromosome")])  # Check the other chr names of other organism 
chr_genes <- unique(gene_table[c("Chromosome")])
chromosomes <- intersect(chr_te[,"chromosome"], chr_genes[,"Chromosome"])

# Long TEs 
# TEs <- subset(TEs, stop - start > 400)
# Short TEs 
TEs <- subset(TEs, stop - start <= 400)

for (chr in chromosomes) {
  print(chr)
  TEs_chr <- subset(TEs, chromosome == chr) 
  genes_chr <- subset(gene_table, Chromosome == chr) 
  
  #Sort by start position for TEs 
  TEs_chr <- TEs_chr[order(TEs_chr$start),]
  
  # Take each gene and check against TE positions
  genes <- genes_chr$orthoID
  
  for (gene in genes){
    # print(gene)
  
    gene_start <- gene_table[gene, "start"]
    gene_end <- gene_table[gene, "stop"]
    distances <- list()
  
    # check against all TEs until cond satisfied
    for (row in 1:nrow(TEs_chr)) {
      #print(row)
      if (length(distances) == 2) {
        break
      }
      # are before the first TE 
      else if (row == 1 && gene_end < TEs_chr[row, "start"]){
        min_distance <- TEs_chr[row, "start"] - gene_end
        break
      }
      #gene after last TE
      else if(row == nrow(TEs_chr) && gene_start > TEs_chr[row, "stop"]){
        min_distance <- gene_start - TEs_chr[row, "stop"]
        break
      }
      # are between two TEs 
      else if (TEs_chr[row, "stop"] < gene_start && TEs_chr[row+1, "start"] > gene_end){
        distance1 <- gene_start - TEs_chr[row, "stop"]
        distance2 <- TEs_chr[row+1, "start"] - gene_end
        distances <- list(distance1, distance2)
        min_distance <- Reduce(min, distances)
        break # should be here ? 
      }
      # gene in TE
      else if (TEs_chr[row, "start"] < gene_start && TEs_chr[row, "stop"] > gene_end){
        min_distance <- 0 
        break
      }
      # first half of gene in TE
      else if (gene_start > TEs_chr[row, "start"] && gene_start < TEs_chr[row, "stop"] &&
               gene_end > TEs_chr[row, "stop"]) {
        min_distance <- 0 
        break
      }
      # last half of gene in TE
      else if (gene_end > TEs_chr[row, "start"] && gene_end < TEs_chr[row, "stop"] &&
               gene_start < TEs_chr[row, "start"]){
        min_distance <- 0 
        break
      }
    } # TE loop end
    # inside gene loop but outside TE loop 
    distance_matrix[gene, "distance"] <- min_distance
    sampletype_sex <- "crassa_BA_sex"   # Fix better
    sampletype_veg <- "crassa_BA_veg"
    distance_matrix[gene, "met_sex"] <- methylation_table[gene, sampletype_sex]
    distance_matrix[gene, "met_veg"] <- methylation_table[gene, sampletype_veg]
  } #gene loop end

} # chr loop end

# long TEs
#distance_matrix_long <- distance_matrix

#Short TEs
distance_matrix_short <- distance_matrix

# Sex tissue 
long <- ggplot(distance_matrix_long, aes(distance, met_sex)) + 
              geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 300000, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 300000, aes(label = ..rr.label..))+
              labs(title="Compare to long TEs", y = "Methylation level",
                   x = "Distance to closest TE")

close_distance_long <- subset(distance_matrix_long, distance <=  1000)
long_close <- ggplot(close_distance_long, aes(distance, met_sex )) + 
              geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 750, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 750, aes(label = ..rr.label..))+
              labs(title="Compare to long TEs", y = "Methylation level",
                   x = "Distance to closest TE")

short <- ggplot(distance_matrix_short, aes(distance, met_sex )) + 
         geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 62500, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 62500, aes(label = ..rr.label..))+
         labs(title="Compare to short TEs", y = "Methylation level",
              x = "Distance to closest TE")

close_distance_short <- subset(distance_matrix_short, distance <=  1000)
short_close <- ggplot(close_distance_short, aes(distance, met_sex )) + 
               geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 750, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 750, aes(label = ..rr.label..))+
               labs(title="Compare to short TEs", y = "Methylation level",
                     x = "Distance to closest TE")
s <- 15
plot_row <- plot_grid(long + theme(text = element_text(size = s)), 
                      long_close + theme(text = element_text(size = s)),
                      short + theme(text = element_text(size = s)), 
                      short_close + theme(text = element_text(size = s)), 
                      labels = c("A", "B", "C", "D"))

title <- ggdraw() + 
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)



# Veg tissue 
long <- ggplot(distance_matrix_long, aes(distance, met_veg)) + 
  geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 300000, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 300000, aes(label = ..rr.label..))+
  labs(title="Compare to long TEs", y = "Methylation level",
       x = "Distance to closest TE") 
  
close_distance_long <- subset(distance_matrix_long, distance <=  1000)
long_close <- ggplot(close_distance_long, aes(distance, met_veg )) + 
  geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 750, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 750, aes(label = ..rr.label..))+
  labs(title="Compare to long TEs", y = "Methylation level",
       x = "Distance to closest TE")

short <- ggplot(distance_matrix_short, aes(distance, met_veg )) + 
  geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 62500, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 62500, aes(label = ..rr.label..))+
  labs(title="Compare to short TEs", y = "Methylation level",
       x = "Distance to closest TE")

close_distance_short <- subset(distance_matrix_short, distance <=  1000)
short_close <- ggplot(close_distance_short, aes(distance, met_veg )) + 
  geom_point()+ geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1.00, label.x = 750, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, label.x = 750, aes(label = ..rr.label..))+
  labs(title="Compare to short TEs", y = "Methylation level",
       x = "Distance to closest TE")

plot_row <- plot_grid(long + theme(text = element_text(size = s)), 
                      close_distance_long + theme(text = element_text(size = s)),
                      short + theme(text = element_text(size = s)), 
                      close_distance_short + theme(text = element_text(size = s)), 
                      labels = c("A", "B", "C", "D"))

title <- ggdraw() + 
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)



