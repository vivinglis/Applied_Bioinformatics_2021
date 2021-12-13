# When trying to push to git and it still asks for credentials
# Use this code in the shell under the git tab in RStudio
# git config remote.origin.url git@github.com:vivinglis/Applied_Bioinformatics_2021.git

getwd()
setwd("../")
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(gtools)

# install.packages('gtools')
# BiocManager::install("apeglm")

############################## TESTS ################################
#max10veg <- read.table("Applied_Bioinformatics_2021/data_tables/1_top10_exp_genes_veg.txt", header = TRUE)
#max10sex <- read.table("Applied_Bioinformatics_2021/data_tables/1_top10_exp_genes_sex.txt", header = TRUE)
#maxexpvals <- gene_expr %>% summarise_if(is.numeric, max)
#testing <- gene_expr %>% group_by()
#gotermstabl <- read.delim("Applied_Bioinformatics_2021/1_GO_translation_table.tsv", header = TRUE)
#transtabl <- read.delim("Applied_Bioinformatics_2021/1_translation_table.tsv", header = TRUE)

####################### DATA PREP ##################################
met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)



########### DESeq2 - Differential Expression analysis ######################
getwd()

# Load RNA count matrix
rna_count_matrix <- read.delim("Data/NeurosporaMethylation-main/data/RNA/RNAcountMatrix.tsv", header = TRUE) 

# Made the coldata file in excel because it was easier.
coldataRNA <- read.csv("Applied_Bioinformatics_2021/coldata_deseq2.csv", header = TRUE, row.names = 1)
rnacount <- rna_count_matrix[,2:61]

# Create a dds object and run DESeq on it
dds <- DESeqDataSetFromMatrix(countData = rnacount, colData = coldataRNA, design = ~ lineage + tissue)
dds <- DESeq(dds)
res <- results(dds, contrast=c("tissue","sex","veg"))
results(dds)
res
resultsNames(dds)

# Shrink for visualisation and ranking of genes according to manual
resLFC <- lfcShrink(dds, coef="tissue_veg_vs_sex", type="apeglm")
resLFC



# Plots that I don't understand
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


############## Table of Ranked Genes ####################

# Create table of genes and rank them (using shrunk LFC)
# So shrink makes veg B and sex A but should be vice versa
# Need to flip the signs (so multiply by -1)

head(resLFC[[2]])
ranked_genes <- cbind(rna_count_matrix[1], "log2foldchange_exp" = -1*(resLFC[[2]]))


########### Find fold change in methylation level #######

# List of genes enriched with conserved methylation

enriched_met <- read.delim("Data/NeurosporaMethylation-main/annotation/genomeAnnotations/Neurospora_crassa.methylated.conserved.genes.tsv", header = TRUE) 
en_met_genes <- c(enriched_met[9])

met$ismetcon <- met$Geneid %in% en_met_genes$OrthoID

# Divide into sex and veg tissue with the info about which are orthoIDs
sex_met <- met[seq_len(ncol(met)-1) %% 2 == 0]
sex_met <- cbind(met[1], sex_met, met[22])
sex_met <- sex_met[order(sex_met$Geneid),]
veg_met <- met[seq_len(ncol(met)) %% 2 == 0 + 1]
veg_met <- cbind(veg_met, met[22])
veg_met <- veg_met[order(veg_met$Geneid),]

# Subset all the rows with values of TRUE in ismetcon column

enrich_sexmet <- subset(sex_met, ismetcon == TRUE)
enrich_vegmet <- subset(veg_met, ismetcon == TRUE)

# Get an average methylation level for each gene
enrich_sexmet$row_mean_sex <- rowMeans(enrich_sexmet[2:11], na.rm=TRUE)
enrich_vegmet$row_mean_veg <- rowMeans(enrich_sexmet[2:11], na.rm=TRUE)


# Ok maybe easier if it's one dataset
met_fc <- cbind(enrich_sexmet[1], enrich_vegmet[13], enrich_sexmet[13])
met_fc <- na.omit(test)
colnames(met_fc) <- c('Geneid', 'veg_met_avg', 'sex_met_avg')

# fold change function, put second number on top....
# then take log2 of this fold change to have same values?

met_fc$foldchange <- log2(foldchange(met_fc$sex_met_avg, met_fc$veg_met_avg))

################# Trying this again ####################
met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)

sex_met <- met[seq_len(ncol(met)) %% 2 == 0]
sex_met <- cbind(met[1], sex_met)
sex_met <- sex_met[order(sex_met$Geneid),]
veg_met <- met[seq_len(ncol(met)) %% 2 == 0 + 1]
veg_met <- veg_met[order(veg_met$Geneid),]


sex_met$row_mean_sex <- rowMeans(sex_met[2:11], na.rm=TRUE)
veg_met$row_mean_veg <- rowMeans(veg_met[2:11], na.rm=TRUE)


# Ok maybe easier if it's one dataset
met_fc <- cbind(met[1], sex_met[12], veg_met[12])
colnames(met_fc) <- c('Geneid', 'veg_met_avg', 'sex_met_avg')

# fold change function, put second number on top....
# then take log2 of this fold change to have same values?

met_fc$foldchange <- foldchange(met_fc$sex_met_avg, met_fc$veg_met_avg)
met_fc$logfc <- foldchange2logratio(met_fc$foldchange, base = 2)

met_fc <- na.omit(met_fc)
met_fc <- met_fc[is.finite(rowSums(met_fc[2:5])),]


densfc <- density(met_fc$logfc)
## calculate the range of the graph
xlim <- range(densfc$x)
ylim <- range(0,densfc$y)
#pick the colours
fcCol <- rgb(0,1,1,0.2)
## plot the carrots and set up most of the plot parameters
plot(densfc, xlim = xlim, ylim = ylim, xlab = 'Expression',
     main = 'Distribution of expression', 
     panel.first = grid())


# Plots are terrible

##### Subset of genes above and below 2 stdevs #######

top2sd <- 2*sd(met_fc$logfc) + mean(met_fc$logfc)
bot2sd <- 2*sd(met_fc$logfc) - mean(met_fc$logfc)

fcmethylated_genes <- subset(met_fc, met_fc$logfc > top2sd &
                               met_fc$logfc < bot2sd)

# Get subset of ranked genes that are fcmethylated and are not
ranked_genes$methfc <- ranked_genes$Geneid %in% fcmethylated_genes$Geneid


ranked_genes_fcmeth <- subset(ranked_genes, methfc == TRUE)
ranked_genes_notmeth <- subset(ranked_genes, methfc == FALSE)


fcmeth <- ranked_genes_fcmeth[[2]]
notmeth <- ranked_genes_notmeth[[2]]

typeof(fcmeth)
is.numeric(fcmeth)




ks.test(notmeth, fcmeth)

# No significant difference! Drawn from the same distribution.

