#Heatmaps - Ward clust, Euclidean distancing - Rocaglate untreated groups 
#Compressed and stretched heatmaps 
#April 20th 2021
#Emily Wood



library(d3heatmap)
library(gplots)
library(pheatmap)
library(heatmaply)

#read the data
dat <- read.csv('total_dataset_com.csv', row.names=1)
dat <- dat[,8:19]
dat <- dat[, c(1,5,9,10)]

#data for annotation and new column stating whether norm_vec was pos or neg
df1 <- read.csv('total_dataset_com.csv', row.names=1)
df1$score <- ifelse(df1$norm_vec >=0, "Positive", "Negative")
df1 <- df1[,c(1,20)]


#scale the data
#dat.tdy <- dat
#dat.n <- scale(t(dat.tdy))
#dat.tn <- t(dat.n)



#annotation for the norm_vec pos/neg
annotationc = as.data.frame(df1[,"score"], row.names=rownames(df1))
STV <- df1$score
norm_vec1 = as.data.frame(STV, row.names=rownames(df1))


#apply euclidean data
d1 <- dist(dat,method = "euclidean", diag = FALSE, upper = FALSE)
round(d1,3)
d2 <- dist(dat,method = "euclidean", diag = FALSE, upper = TRUE)


# Clustering distance between experiments using Ward linkage
c1 <- hclust(d1, method = "ward.D2", members = NULL)
# Clustering distance between proteins using Ward linkage
c2 <- hclust(d2, method = "ward.D2", members = NULL)


# Check clustering by plotting dendrograms
par(mfrow=c(2,1),cex=0.5) # Make 2 rows, 1 col plot frame and shrink labels


png(file="RocaglateTreat_Hmap.png", width = 800, height = 28000, units = "px", bg="white") 
par(cex.main=0.75) # Shrink title fonts on plot
pheatmap(dat,scale="row",                   # Tidy, normalised data
         Colv=as.dendrogram(c1),     # Experiments clusters in cols
         Rowv=as.dendrogram(c2),     # Protein clusters in rows
         density.info="histogram",   # Plot histogram of data and colour key
         trace="none",               # Turn of trace lines from heat map
         cexRow=0.5,cexCol=0.75, keysize=0.75,    # Amend row and column label fonts
         main = "Log2FC Rocaglate Drug perturbations Heatmap", annotation_row = norm_vec1, angle_col = 90)

# Close file
dev.off()


#compressed euclidean hmap, no gene names
png(file="RocaglateTreatedcompressed.png", width = 400, height = 800, bg = "white") 
par(cex.main=0.75) # Shrink title fonts on plot
pheatmap(dat,scale="row" ,                    # Tidy, normalised data
         Colv=as.dendrogram(c1),     # Experiments clusters in cols
         Rowv=as.dendrogram(c2),     # Protein clusters in rows
         density.info="histogram",   # Plot histogram of data and colour key
         trace="none",  show_rownames= F,show_colnames=T,             # Turn of trace lines from heat map
         cexRow=0.5,cexCol=0.75, keysize=0.75,    # Amend row and column label fonts
         main = "Log2FC Rocaglate Drug perturbations Heatmap" ,annotation_row = norm_vec1, angle_col = 90)

# Close file
dev.off()
