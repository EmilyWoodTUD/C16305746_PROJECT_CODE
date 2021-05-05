#Heatmaps using ward clustering and correlation distancing
#Emily Wood
#22nd of April 2021


library(ComplexHeatmap)
library(gplots)
library(pheatmap)
library(heatmaply)


#load in the data and make annotation dat for the norm_vec value
dat <- read.csv('total_dataset_com copy.csv', row.names=1)

#data for annotation and new column stating whether norm_vec was pos or neg
df1 <- read.csv('total_dataset_com.csv', row.names=1)
df1$score <- ifelse(df1$norm_vec >=0, "Positive", "Negative")
df1 <- df1[,c(1,20)]


#annotation for the norm_vec pos/neg
annotationc = as.data.frame(df1[,"score"], row.names=rownames(df1))
STV <- df1$score
norm_vec1 = as.data.frame(STV, row.names=rownames(df1))


# Pairwise correlation between samples (columns)
cols.cor <- cor(dat, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(dat), use = "pairwise.complete.obs", method = "pearson")


## Row- and column-wise clustering using correlation 
hclust.col <- hclust(as.dist(1-cols.cor)) 
hclust.row <- hclust(as.dist(1-rows.cor))


#writing the file and generating the heatmap using correlation distancing and ward clustering
png(file="correlationdist.png", width = 900, height = 28000, units = "px", bg = "white") 
pheatmap(dat,scale="row",
         trace = "none", density.info = "none",
         Colv = as.dendrogram(hclust.col),
         Rowv = as.dendrogram(hclust.row)
         ,show_rownames= T,show_colnames=T,
         density.info="histogram",   # Plot histogram of data and colour key
         trace="none",               # Turn of trace lines from heat map
         cexRow=0.5,cexCol=0.75, keysize=0.75,    # Amend row and column label fonts
         main = "Gene Contributions Correaltion distance", annotation_row = norm_vec1)


# Close file
dev.off()

#compressed heatmap with no gene names
png(file="compressed.correlation.png", width = 400, height = 800, bg = "white") 

pheatmap(dat,scale="row",
         trace = "none", density.info = "none",
         Colv = as.dendrogram(hclust.col),
         Rowv = as.dendrogram(hclust.row)
         ,show_rownames= F,show_colnames=T,
         density.info="histogram",   # Plot histogram of data and colour key
         trace="none",               # Turn of trace lines from heat map
         cexRow=0.5,cexCol=0.75, keysize=0.75,    # Amend row and column label fonts
         main = "Gene Contributions Correaltion distance", annotation_row = norm_vec1)


# Close file
dev.off()

