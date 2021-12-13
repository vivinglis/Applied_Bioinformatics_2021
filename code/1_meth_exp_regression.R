# When trying to push to git and it still asks for credentials
# Use this code in the shell under the git tab in RStudio
# git config remote.origin.url git@github.com:vivinglis/Applied_Bioinformatics_2021.git

getwd()
setwd("../")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(lattice)

######################## Data prep ###################################
# Load expression and methylation data tables

gene_expr <- read.table("Applied_Bioinformatics_2021/data_tables/1_rlog_normalised_avgs.tsv", header = TRUE) 
met <- read.table("Applied_Bioinformatics_2021/data_tables/1_orthologous_genes_methylation_with_promotor_region.tsv", header = TRUE)

# Expression data divided into tissues (already ordered by gene ID)

sex_exp <- gene_expr[seq_len(ncol(gene_expr)) %% 2 == 0]
sex_exp <- cbind(gene_expr[1],sex_exp)
veg_exp <- gene_expr[seq_len(ncol(gene_expr)) %% 2 == 0 + 1]

# Methylation data divided into tissues and ordered on gene ID

sex_met <- met[seq_len(ncol(met)) %% 2 == 0]
sex_met <- cbind(met[1], sex_met)
sex_met <- sex_met[order(sex_met$Geneid),]
veg_met <- met[seq_len(ncol(met)) %% 2 == 0 + 1]
veg_met <- veg_met[order(veg_met$Geneid),]


########### Testing making the dataset ################

#test_data <- sex_met[1]
#test_data <- cbind(test_data,tissue = rep("sex",7109),lineage= rep("crassa_BA",7109))
#test_data <- cbind(test_data,methylation= sex_met[2])
#names(test_data)[names(test_data) == "crassa_BA_sex"] <- "methylation"
#test_data <- test_data[order(test_data$Geneid),]

#sex_exp <- sex_exp[order(sex_exp$Geneid),]
#test_data <- cbind(test_data, sex_exp[2])
#names(test_data)[names(test_data) == "crassa_BA_sex"] <- "expression"

# Make a loop and do it for each lineage
# Also try making one big table and have expression ~ methylation + Geneid


######################################### SEX TISSUE ######################################

# Lineages in sex

lins <- c("crassa_BA_sex", "crassa_sa_sex", "L1_BA_sex", "L1_sa_sex", 
          "L10_BA_sex", "L10_sa_sex", "L6_BA_sex", "L6_sa_sex",
          "sito_BA_sex", "sito_sa_sex")


# Create datasets from each lineage in sex_met

for(i in 2:ncol(sex_met)){
  temp <- data.frame(sex_met[,i])
  colnames(temp) <- "methylation"
  assign(colnames(sex_met)[i], temp)
  rm(temp)
}

# Make into a list and rename list elements according to lineage, add gene IDs
sex_list <- list(crassa_BA_sex, crassa_sa_sex, L1_BA_sex, L1_sa_sex, 
                 L10_BA_sex, L10_sa_sex, L6_BA_sex, L6_sa_sex,
                 sito_BA_sex, sito_sa_sex)
names(sex_list) <- lins
sex_list <- lapply(sex_list, add_column, sex_met[1], .before = 1)


# Add expression data and change column name to expression
sex_list_indx <- 1

for(sexitems in sex_list){
  sex_list[[sex_list_indx]] <- cbind(sex_list[[sex_list_indx]], sex_exp[sex_list_indx+1])
  colnames(sex_list[[sex_list_indx]])[3] <- "expression"
  sex_list_indx <- sex_list_indx + 1
}

# Regression time!
reg_list <- list()
for(elements in sex_list){
  fit <- lm(expression ~ methylation, data=elements)
  print(summary(fit))
  reg_list <- cbind(reg_list, summary(fit))
}

data_reg <- as.data.frame(reg_list)
colnames(data_reg) <- lins


# And now to plot stuff!
plot_list = list()
pos <- 1
listpos <- 1
for(stuff in sex_list){
  plotstuff <- ggplot(stuff, aes(methylation, expression))+
    geom_point(na.rm = TRUE)+
    labs(title = names(sex_list)[listpos])+
    geom_smooth(method=lm, na.rm = TRUE)+
    stat_regline_equation(label.y = 16, label.x = 0.5, aes(label = ..rr.label..), na.rm = TRUE)
  plot_list[[pos]] <- plotstuff
  pos <- pos + 1
  listpos <- listpos + 1
  
}
title = textGrob("SEX: Expression ~ Methylation",gp=gpar(fontsize=15,font=3))
grid.arrange(grobs = plot_list, ncol = 5, top = title )

tit <- expression(paste("Methylation and expression in sex tissue of ", 
                        italic("Neurospora crassa"), " (BA)"))


# Plot of just N. crassa for popular science presentation
popsciplot <- ggplot(sex_list[[1]], aes(methylation, expression))+
  geom_point(na.rm = TRUE)+
  labs(title = tit)+
  geom_smooth(method=lm, na.rm = TRUE)+
  stat_regline_equation(label.y = 16, label.x = 0.75, aes(label = ..rr.label..), na.rm = TRUE)

popsciplot

######################################### VEG TISSUE ######################################

veglins <- c("crassa_BA_veg", "crassa_sa_veg", "L1_BA_veg", "L1_sa_veg", 
          "L10_BA_veg", "L10_sa_veg", "L6_BA_veg", "L6_sa_veg",
          "sito_BA_veg", "sito_sa_veg")


# Create datasets from each lineage in veg_met

for(j in 2:ncol(veg_met)){
  temp <- data.frame(veg_met[,j])
  colnames(temp) <- "methylation"
  assign(colnames(veg_met)[j], temp)
  rm(temp)
}

# Make into a list and rename list elements according to lineage, add gene IDs
veg_list <- list(crassa_BA_veg, crassa_sa_veg, L1_BA_veg, L1_sa_veg, 
                 L10_BA_veg, L10_sa_veg, L6_BA_veg, L6_sa_veg,
                 sito_BA_veg, sito_sa_veg)
names(veg_list) <- veglins
veg_list <- lapply(veg_list, add_column, veg_met[1], .before = 1)


# Add expression data and change column name to expression
veg_list_indx <- 1

for(vegitems in veg_list){
  veg_list[[veg_list_indx]] <- cbind(veg_list[[veg_list_indx]], veg_exp[veg_list_indx+1])
  colnames(veg_list[[veg_list_indx]])[3] <- "expression"
  veg_list_indx <- veg_list_indx + 1
}


# Regression time!
reg_list_veg <- list()
for(vegelements in veg_list){
  fit_veg <- lm(expression ~ methylation, data=vegelements)
  print(summary(fit_veg))
  reg_list_veg <- cbind(reg_list_veg, summary(fit_veg))
}

data_reg_veg <- as.data.frame(reg_list_veg)
colnames(data_reg_veg) <- veglins


# And now to plot stuff!
plot_list_veg = list()
vpos <- 1
vlistpos <- 1
for(vstuff in veg_list){
  vplotstuff <- ggplot(vstuff, aes(methylation, expression))+
    geom_point(na.rm = TRUE)+
    labs(title = names(veg_list)[vlistpos])+
    geom_smooth(method=lm, na.rm = TRUE)+
    stat_regline_equation(label.y = 16, label.x = 0.5, aes(label = ..rr.label..), na.rm = TRUE)
  plot_list_veg[[vpos]] <- vplotstuff
  vpos <- vpos + 1
  vlistpos <- vlistpos + 1
  
}

grid.arrange(grobs = plot_list_veg, ncol = 5, top = textGrob("VEG: Expression ~ Methylation",gp=gpar(fontsize=15,font=3)))


########## Exporting data tables #######

#write.table(sex_met, file="Applied_Bioinformatics_2021/data_tables/1_sex_methylation", sep="\t", quote=F, col.names=NA)
#write.table(sex_exp, file="Applied_Bioinformatics_2021/data_tables/1_sex_expression", sep="\t", quote=F, col.names=NA)
#write.table(veg_met, file="Applied_Bioinformatics_2021/data_tables/1_veg_methylation", sep="\t", quote=F, col.names=NA)
#write.table(veg_exp, file="Applied_Bioinformatics_2021/data_tables/1_veg_expression", sep="\t", quote=F, col.names=NA)


############# 2 Standard deviations ############
# Functions to calculate the values for 2 standard deviations above mean (top)
# and for 2 standard deviations below the mean (bot) for methylation levels

sd2_top <- function(dataset){
  two_top_sds <- mean(unlist(dataset), na.rm = TRUE) + 2*sd(unlist(dataset), na.rm = TRUE)
  return(two_top_sds)
}

sd2_bot <- function(dataset){
  two_bot_sds <- mean(unlist(dataset), na.rm = TRUE) - 2*sd(unlist(dataset), na.rm = TRUE)
  return(two_bot_sds)
}

sex_listsd <- list(crassa_BA_sex, crassa_sa_sex, L1_BA_sex, L1_sa_sex, 
                  L10_BA_sex, L10_sa_sex, L6_BA_sex, L6_sa_sex,
                  sito_BA_sex, sito_sa_sex)
names(sex_listsd) <- lins

veg_listsd <- list(crassa_BA_veg, crassa_sa_veg, L1_BA_veg, L1_sa_veg, 
                  L10_BA_veg, L10_sa_veg, L6_BA_veg, L6_sa_veg,
                  sito_BA_veg, sito_sa_veg)
names(veg_listsd) <- veglins



# Get a list of the 2stdevs above and below the mean for both tissues (methylation)
# Sex lists for top 2 sds and bottom 2 sds
sdlists_bot <- list()
index <- 1
for(lineages in sex_listsd){
  sdlists_bot[[index]] <- sd2_bot(lineages)
  index <- index + 1
}

sdlists_top <- list()
index <- 1
for(lineages in sex_listsd){
  sdlists_top[[index]] <- sd2_top(lineages)
  index <- index + 1
}

names(sdlists_bot) <- lins
names(sdlists_top) <- lins


# Veg lists for top 2sds and bottom 2 sds
sdlistv_bot <- list()
index <- 1
for(lineagesv in veg_listsd){
  sdlistv_bot[[index]] <- sd2_bot(lineagesv)
  index <- index + 1
}

sdlistv_top <- list()
index <- 1
for(lineagesv in veg_listsd){
  sdlistv_top[[index]] <- sd2_top(lineagesv)
  index <- index + 1
}

names(sdlistv_bot) <- veglins
names(sdlistv_top) <- veglins



################################ ABOVE 2 STDEVS ###############################
# For each item in e.g. the sex_list go through it and only plot the numbers
# that are above the threshold value as specified in the top list

# SEX
linsnames <- lins
tissuelist <- sex_list
top_list <- sdlists_top
bot_list <- sdlists_bot
# VEG
linsnames <- veglins
tissuelist <- veg_list
top_list <- sdlistv_top
bot_list <- sdlistv_bot

# Create subsets of the data
index_list <- 1:10
subset_list_top <- list()
index <- 1
idx <- 1
for(items in index_list){
  subset_list_top[[idx]] <- subset(tissuelist[[items]], methylation >= top_list[[index]])
  idx <- idx + 1
  index <- index + 1
}

names(subset_list_top) <- linsnames


# Plot this stuff
listpos1 <- 1
pos1 <- 1
plot_list_tops <- list()
for(valst in subset_list_top){
  topplots <- ggplot(valst, aes(methylation, expression))+
    geom_point(na.rm = TRUE)+
    labs(title = names(subset_list_top)[listpos1])+
    geom_smooth(method=lm, na.rm = TRUE)+
    stat_regline_equation(label.y = 16, label.x = 0.75, aes(label = ..rr.label..), na.rm = TRUE)
  plot_list_tops[[pos1]] <- topplots
  pos1 <- pos1 + 1
  listpos1 <- listpos1 + 1
  
}

grid.arrange(grobs = plot_list_tops, ncol = 5, top = textGrob("Above 2 stdevs: Expression ~ Methylation",gp=gpar(fontsize=15,font=3)))


###################### BELOW 2 STDEVS ########################################


# Create subsets of the data
subset_list_bot <- list()
index <- 1
idx <- 1
for(items in index_list){
  subset_list_bot[[idx]] <- subset(tissuelist[[items]], methylation <= bot_list[[index]])
  idx <- idx + 1
  index <- index + 1
}

names(subset_list_bot) <- linsnames


# Plot this stuff
listpos1 <- 1
pos1 <- 1
plot_list_bots <- list()
for(botnons in subset_list_bot){
  botplots <- ggplot(botnons, aes(methylation, expression))+
    geom_point(na.rm = TRUE)+
    labs(title = names(subset_list_bot)[listpos1])+
    geom_smooth(method=lm, na.rm = TRUE)+
    stat_regline_equation(label.y = 15, label.x = 0, aes(label = ..rr.label..), na.rm = TRUE)+
    xlim(0,0.17)
  plot_list_bots[[pos1]] <- botplots
  pos1 <- pos1 + 1
  listpos1 <- listpos1 + 1
  
}

grid.arrange(grobs = plot_list_bots, ncol = 5, top = textGrob("Below 2 stdevs: Expression ~ Methylation",gp=gpar(fontsize=15,font=3)))



################### -2 STDEVS <-> +2 STDEVS ################

# Create subset of the data
mean_met <- list()
topidx <- 1
botidx <- 1
idx <- 1
for(meanmetitems in index_list){
  mean_met[[idx]] <- subset(tissuelist[[meanmetitems]], 
                            methylation >= bot_list[[botidx]] & methylation <= top_list[[topidx]])
  idx <- idx + 1
  botidx <- botidx + 1
  topidx <- topidx + 1
}

names(mean_met) <- linsnames


# Plot this stuff
listpos1 <- 1
pos1 <- 1
plot_list_mean <- list()
for(meanitems in mean_met){
  meanplots <- ggplot(meanitems, aes(methylation, expression))+
    geom_point(na.rm = TRUE)+
    labs(title = names(mean_met)[listpos1])+
    geom_smooth(method=lm, na.rm = TRUE)+
    stat_regline_equation(label.y = 17, label.x = 0.15, aes(label = ..rr.label..), na.rm = TRUE)+
    xlim(0.15, 0.30)
  plot_list_mean[[pos1]] <- meanplots
  pos1 <- pos1 + 1
  listpos1 <- listpos1 + 1
  
}

grid.arrange(grobs = plot_list_mean, ncol = 5, top = textGrob("Mean Region: Expression ~ Methylation",gp=gpar(fontsize=15,font=3)))


###### Save DATA TABLES #####
#getwd()
#write.table
#write.table(mean_met, file="Applied_Bioinformatics_2021/data_tables/viv_meanmet.txt", sep=";")

#capture.output(mean_met, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/met_mean.txt")
#capture.output(subset_list_top, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/met_high.txt")
#capture.output(subset_list_bot, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/met_low.txt")

#options(max.print=1000000)

#topgenes <- lapply(subset_list_top, `[[`, 1)
#botgenes <- lapply(subset_list_bot, `[[`, 1)
#midgenes <- lapply(mean_met, `[[`, 1)

#capture.output(topgenes, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/topgenes.txt")
#capture.output(botgenes, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/botgenes.txt")
#capture.output(midgenes, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/midgenes.txt")

#options(max.print = 1000)

#capture.output(topgenes, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/topgenes_v.txt")
#capture.output(botgenes, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/botgenes_v.txt")
#capture.output(midgenes, file = "Applied_Bioinformatics_2021/data_tables/viv_till_olivia/midgenes_v.txt")

######################### Cumulative Distribution Graphs #################

# Popsci plot
a <- subset_list_top[[1]][[3]]
b <- subset_list_bot[[1]][[3]]
c <- mean_met[[1]][[3]]

tit_ecdf <- expression(paste("Expression ECDF for Methylation Levels in ", italic("Neurospora crassa")))


plot(ecdf(a), verticals=TRUE, do.points=FALSE, col='royalblue4',
     main=tit_ecdf, xlab="Expression", ylab="Density")
plot(ecdf(b), verticals=TRUE, do.points=FALSE, add=TRUE, col='red')
plot(ecdf(c), verticals=TRUE, do.points=FALSE, add=TRUE, col='sienna1')
legend('topleft',c('Mean methylation','High methylation', 'Low methylation'),
       fill = c('sienna1', 'royalblue4', 'red'), bty = 'n',
       border = NA)


# Actual plots
# create a data set containing the range you wish to use
min(gene_expr[2:11], na.rm = TRUE)
max(gene_expr[2:11], na.rm = TRUE)
axisdf <- data.frame(xax = c(-3,18))
# create a list of calls to `stat_function` with the colours you wish to use

# Try again with ggplot2
subset_idx <- 1:10
listpos1 <- 1
pos1 <- 1
cplot_list <- list()
for(subsets in subset_idx){
  a <- subset_list_top[[subsets]][[3]]
  b <- subset_list_bot[[subsets]][[3]]
  c <- mean_met[[subsets]][[3]]
  ecdflist <- Map(f  = stat_function, colour = c('royalblue4', 'red', 'sienna1'),
                  fun = list(ecdf(a), ecdf(b), ecdf(c)), geom = 'step')
  cplot <- ggplot(data = axisdf, aes(x = xax)) + 
    ecdflist +
    labs(title = names(mean_met)[listpos1])+
    xlab("Expression")+
    ylab("Density")+
    theme(legend.position = "none")
  cplot_list[[pos1]] <- cplot
  pos1 <- pos1 + 1
  listpos1 <- listpos1 + 1
}

grid.arrange(grobs = cplot_list, ncol = 5, 
             top = textGrob("ECDF: Expression for Varying Methylation", gp=gpar(fontsize=15,font=3)))

# For interpretation
# Red = low methylation
# Blue = high methylation
# Yellow = average methylation


##### Kolmogorov-Smirnov Tests #####

# install.packages('dgof')
library(dgof)

testx <- mean_met[[5]][[3]]
testy <- subset_list_top[[5]][[3]]

head(testx)
head(testy)
ks.test(testx, testy, alternative = "less")
#   cannot compute correct p-values with ties?

# Basically, yes, even the closest one is significantly lower than mean.
# Wanted to see that highly methylated was significantly lower expressed than averagely methylated.
# That's why I chose alternative = 'less'.


# A nice loop to calculate the same thing for all lineages
mean_met
subset_list_top
subset_list_bot

# Mid-high (less)
# Mid-low (greater), mid-low (twosided), mid-low(less)
# High-low (greater)
# Run this loop for veg and sex and all the options listed above
# index_list is a list of numbers from 1:10


ksvals <- data.frame()
idx <- 1
for(linnum in index_list){
  x <- mean_met[[idx]][[3]]
  y <- subset_list_top[[idx]][[3]]
  pval <- ks.test(x, y, alternative = 'less')[[2]]
  ksvals <- rbind(ksvals, c(veglins[idx], pval))
  idx <- idx + 1

}
colnames(ksvals) <- c('lineages', 'p-value')


#write.table(ksvals, file="Applied_Bioinformatics_2021/data_tables/ksvals_veg_mh.txt", sep=";")



###### Geometric Density ######

# Stolen code from the internet
# https://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r

library(grid)
library(gridGraphics)


densplotlist <- list()
idx <- 1
for(lineanum in index_list){
  testx <- mean_met[[lineanum]][[3]]
  testy <- subset_list_top[[lineanum]][[3]]
  testz <- subset_list_bot[[lineanum]][[3]]
  
  ## calculate the density - don't plot yet
  densCarrot <- density(testx, na.rm = TRUE)
  densCuke <- density(testy, na.rm = TRUE)
  densCauli <- density(testz, na.rm = TRUE)
  ## calculate the range of the graph
  xlim <- range(-6, 16)
  ylim <- range(0,0.25)
  #pick the colours
  carrotCol <- rgb(1,0,0,0.2)
  cukeCol <- rgb(0,0,1,0.2)
  cauliCol <- rgb(0,1,1,0.2)
  ## plot the carrots and set up most of the plot parameters
  plot(densCarrot, xlim = xlim, ylim = ylim, xlab = 'Expression',
       main = names(mean_met)[idx]) #panel.first = grid())
  #put our density plots in
  polygon(densCarrot, density = -1, col = carrotCol)
  polygon(densCuke, density = -1, col = cukeCol)
  polygon(densCauli, density = -1, col = cauliCol)
  ## add a legend in the corner
  #legend('topleft',c('Mean methylation','High methylation', 'Low methylation'),
  #       fill = c(carrotCol, cukeCol, cauliCol), bty = 'n',
  #       border = NA)
  
  #grid.echo()
  #a <- grid.grab()
  #densplotlist[[idx]] <- a
  #idx <- idx + 1
}

typeof(densplotlist)

par(mfrow=c(2,5))

grid.arrange(grobs = densplotlist[1:5], ncol = 3, nrow = 2)

