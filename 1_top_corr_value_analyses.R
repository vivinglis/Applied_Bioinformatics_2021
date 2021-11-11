# Top hundred genes with correlation value close to -1. 

# Packages
# install.packages("stringi")
library("stringi")

# Read in tables
pearson_corr_avgs <- read.table("1_2_pearson_corr_avgs.tsv", header = TRUE) 
pearson_corr_avgs_with_pro <- read.table("1_2_pearson_corr_avgs_with_pro.tsv", header = TRUE) 
pearson_corr_rep <- read.table("1_2_pearson_corr_rep.tsv", header = TRUE) 
pearson_corr_rep_with_pro <- read.table("1_2_pearson_corr_rep_with_pro.tsv", header = TRUE) 

# Take top number values from each dataframe 
top_num <- 100
top_pearson_corr_avgs <- head(pearson_corr_avgs[order(-pearson_corr_avgs$Value, decreasing =TRUE),], top_num)
top_pearson_corr_avgs_with_pro <- head(pearson_corr_avgs_with_pro[order(-pearson_corr_avgs_with_pro$Value, decreasing =TRUE),], top_num)
top_pearson_corr_rep <- head(pearson_corr_rep[order(-pearson_corr_rep$Value, decreasing =TRUE),], top_num)
top_pearson_corr_rep_with_pro <- head(pearson_corr_rep_with_pro[order(-pearson_corr_rep_with_pro$Value, decreasing =TRUE),], top_num)

# Gener samma 
a <- intersect(top_pearson_corr_avgs$Geneid, top_pearson_corr_avgs_with_pro$Geneid)
b <- intersect(top_pearson_corr_rep$Geneid, top_pearson_corr_rep_with_pro$Geneid)
c <- intersect(a,b)

# replicat and average utan promotor 
d <- intersect(top_pearson_corr_avgs$Geneid, top_pearson_corr_rep$Geneid)

# replicat and average med promotor 
e <- intersect(top_pearson_corr_avgs_with_pro$Geneid, top_pearson_corr_rep_with_pro$Geneid)

# Plot sex mot veg alla gene

avgs <- list(subset(pearson_corr_avgs, Tissue == "veg")$Value,
             subset(pearson_corr_avgs, Tissue == "sex")$Value)
avgs_p <- list(subset(pearson_corr_avgs_with_pro, Tissue == "veg")$Value,
               subset(pearson_corr_avgs_with_pro, Tissue == "sex")$Value)
rep <- list(subset(pearson_corr_rep, Tissue == "veg")$Value,
            subset(pearson_corr_rep, Tissue == "sex")$Value)
rep_p <- list(subset(pearson_corr_rep_with_pro, Tissue == "veg")$Value,
            subset(pearson_corr_rep_with_pro, Tissue == "sex")$Value)
  
  
par(mfrow=c(2,2))
stripchart(avgs,
           main="Pearson correlation value: avgs",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE
)
stripchart(avgs_p,
           main="Pearson correlation value: avgs with promotor",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE
)
stripchart(rep,
           main="Pearson correlation value: replicas",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE
)
stripchart(rep_p,
           main="Pearson correlation value: replicas with promotor",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE
)


# Sex mot veg for top 100
avgs <- list(subset(top_pearson_corr_avgs, Tissue == "veg")$Value,
             subset(top_pearson_corr_avgs, Tissue == "sex")$Value)
avgs_p <- list(subset(top_pearson_corr_avgs_with_pro, Tissue == "veg")$Value,
               subset(top_pearson_corr_avgs_with_pro, Tissue == "sex")$Value)
rep <- list(subset(top_pearson_corr_rep, Tissue == "veg")$Value,
            subset(top_pearson_corr_rep, Tissue == "sex")$Value)
rep_p <- list(subset(top_pearson_corr_rep_with_pro, Tissue == "veg")$Value,
              subset(top_pearson_corr_rep_with_pro, Tissue == "sex")$Value)

par(mfrow=c(2,2))
stripchart(avgs,
           main="Pearson correlation value: avgs",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE,
           ylim = c(-1, -0.9)
)
stripchart(avgs_p,
           main="Pearson correlation value: avgs with promotor",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE, 
           ylim = c(-1, -0.9)
)
stripchart(rep,
           main="Pearson correlation value: replicas",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE,
           ylim = c(-1, -0.9)
)
stripchart(rep_p,
           main="Pearson correlation value: replicas with promotor",
           xlab="",
           ylab="Correlation value",
           method="jitter",
           col=c("green", "blue"),
           pch=1, 
           vertical=TRUE,
           ylim = c(-1, -0.9)
)


