## 2022.10.14 KIM JS
## trial for clustering & heatmap by using pheatmap

############################################################################
######################## Z-standardization & heapmap
library(pheatmap)
library(viridis) ##for a color palette
library(factoextra) ##for elbow plots finding-out optimal clusters
library(cluster) ## for Clusgap() method

#### Excel-prepared -log10 P-values 
d <- read.csv("sample.csv", header=T, row.names=1)
head(d, c(5,5))
#       Day.0_UP  Day.1.DR_UP Day.3.DR_UP Day.4.DR_UP Day.1.WTR_UP
#10 0.9868054860 0.4172589920 0.006833294 0.006492316  1.000000000
#20 0.6293700630 0.3582790790 0.357276872 0.151149394  1.000000000
#30 1.0000000000 0.7260749530 0.177509484 0.028276453  0.867369286
#40 0.0000176386 0.0000189053 0.165929924 0.987246154  0.000296418
#51 1.0000000000 0.4434475220 0.242329370 0.087546826  1.000000000
d1 <- d # as native data
d2 <- -log(d, 10) #log transformation

#Z-standardization
d1z <- t(apply(d1, 1, scale)) #need transpose
colnames(d1z) <- colnames(d1) #z-scaling losing colnames
d1z <- na.omit(d1z) #127 rows -> 106 rows

d2z <- t(apply(d2, 1, scale))
colnames(d2z) <- colnames(d2)
d2z <- na.omit(d2z) #127 rows -> 106 rows

## visualization
#different colorscheme like magma(100), viridis(100) also should be considered
palette <- colorRampPalette(c("blue", "gray88","gray88","gray88", "red"))(100) # for custom palette
pheatmap(d1z, color = inferno(100), cluster_cols = F, clustering_method = "ward.D2", 
         clustering_distance_rows = "correlation", cutree_rows = 8)
#the plot saved as PDF as a portrait (12 x 6 in) for KEGG
#the plot saved as PDF as a portrait (10 x 6 in) for GO-MF



###### Find-out Optimal clusters
##https://www.r-bloggers.com/2019/01/10-tips-for-choosing-the-optimal-number-of-clusters/
#elbow plot 
fviz_nbclust(d1z, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")
#silhouette plot
fviz_nbclust(e1z, kmeans, method = "silhouette", k.max = 24) 
+ theme_minimal() + ggtitle("The Silhouette Plot")
#gap stats
gap_stat <- clusGap(d1z, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")
