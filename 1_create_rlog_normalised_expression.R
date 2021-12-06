getwd()
setwd("C:/Users/vivia/Documents/Master's Bioinformatics YEAR 2/Applied Bioinformatics/Project/Data")

# Packageds and libraries that I had to load / download
library(ape)
install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")
library(tidyverse)

n_crassa_annot <- read.gff("NeurosporaMethylation-main/annotation/genomeAnnotations/Neurospora_crassa.NC12.gff")

# Load Hugo's translation table
transtab <- read.delim("../Applied_Bioinformatics_2021/1_translation_table.tsv", header = TRUE)

# Load RNA count matrix
rna_count_matrix <- read.delim("NeurosporaMethylation-main/data/RNA/RNAcountMatrix.tsv", header = TRUE) 
type(rna_count_matrix)

# Made the coldata file in excel because it was easier.
coldataRNA <- read.csv("../coldata.csv", header = TRUE, row.names = 1)
mod_rnacount <- rna_count_matrix[,2:61]

# Check that all colnames are in the coldataRNA file
all(colnames(mod_rnacount) %in% rownames(coldataRNA))


# DESeq2 normalisation 
dds <- DESeqDataSetFromMatrix(countData = mod_rnacount, colData = coldataRNA, design = ~ lineage + tissue)
rld <- rlog(dds)
max(rld)

# Create a table of normalised counts for all samples
rld_mod <- assay(rld)
genes <- rna_count_matrix[1]
rld_mod <- cbind(genes,rld_mod)
write.table(rld_mod, file="../Applied_Bioinformatics_2021/1_rlog_normalised_allsamples.tsv", sep="\t", quote=F, col.names=NA)

max(rld_mod)


rld_mod_avg %>% summarise_if(is.numeric, max)


# Same table but with averages (not the individual samples)
type(rld_mod)
rld_mod <- as_tibble(rld_mod)

rld_mod_avg <- rld_mod %>% 
  transmute(Geneid, crassa_BA_sex = (crassa_BA_P1_sex+crassa_BA_P2_sex+crassa_BA_P3_sex) / 3,
            crassa_BA_veg = (crassa_BA_P1_veg+crassa_BA_P2_veg+crassa_BA_P3_veg) / 3,
            crassa_sa_sex = (crassa_sa_P1_sex+crassa_sa_P2_sex+crassa_sa_P3_sex) / 3,
            crassa_sa_veg = (crassa_sa_P1_veg+crassa_sa_P2_veg+crassa_sa_P3_veg) / 3,
            L1_BA_sex = (L1_BA_P1_sex+L1_BA_P2_sex+L1_BA_P3_sex) / 3,
            L1_BA_veg = (L1_BA_P1_veg+L1_BA_P2_veg+L1_BA_P3_veg) / 3,
            L1_sa_sex = (L1_sa_P1_sex+L1_BA_P2_sex+L1_sa_P3_sex) / 3,
            L1_sa_veg = (L1_sa_P1_veg+L1_BA_P2_veg+L1_sa_P3_veg) / 3,
            L10_BA_sex = (L10_BA_P1_sex+L10_BA_P2_sex+L10_BA_P3_sex) / 3,
            L10_BA_veg = (L10_BA_P1_veg+L10_BA_P2_veg+L10_BA_P3_veg) / 3,
            L10_sa_sex = (L10_sa_P1_sex+L10_sa_P2_sex+L10_sa_P3_sex) / 3,
            L10_sa_veg = (L10_sa_P1_veg+L10_sa_P2_veg+L10_sa_P3_veg) / 3,
            L6_BA_sex = (L6_BA_P1_sex+L6_BA_P2_sex+L6_BA_P3_sex) / 3,
            L6_BA_veg = (L6_BA_P1_veg+L6_BA_P2_veg+L6_BA_P3_veg) / 3,
            L6_sa_sex = (L6_sa_P1_sex+L1_sa_P2_sex+L1_sa_P3_sex) / 3,
            L6_sa_veg = (L6_sa_P1_veg+L1_sa_P2_veg+L1_sa_P3_veg) / 3,
            sito_BA_sex = (sito_BA_P1_sex+sito_BA_P2_sex+sito_BA_P3_sex) / 3,
            sito_BA_veg = (sito_BA_P1_veg+sito_BA_P2_veg+sito_BA_P3_veg) / 3,
            sito_sa_sex = (sito_sa_P1_sex+sito_sa_P2_sex+sito_sa_P3_sex) / 3,
            sito_sa_veg = (sito_sa_P1_veg+sito_sa_P2_veg+sito_sa_P3_veg) / 3)

write.table(rld_mod_avg, file="../Applied_Bioinformatics_2021/1_rlog_normalised_avgs.tsv", sep="\t", quote=F, col.names=NA)


# When trying to push to git and it still asks for credentials
# Use this code in the shell under the git tab in RStudio
# git config remote.origin.url git@github.com:vivinglis/Applied_Bioinformatics_2021.git