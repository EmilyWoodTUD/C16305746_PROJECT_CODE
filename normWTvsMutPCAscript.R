#This is the PCA plot for wildtupe 1h versus Mutant B6 with respect to controls - all genes

#packages for PCA plot
library(gplots) 
library(preprocessCore)
library(factoextra)

#load in the data
data <- read.delim(file = "RNA_seq_TNF_full_new_Nov2018.txt", header = TRUE, row.names = 1)
data <- data[-1]

#convert to numerical data
y <- data.matrix(data)

#vectors containing metadata 
#mutant dataset - B6
mutant <- y[, c(1, 2, 3, 9, 10, 11)]
#Wild-type dataset - 1h
WT <- y[, c(5,6,7,13,14,15)]
#metadata of the control sets for 1h - 
controlWT <- rowMeans(y[,c(4,12)])
#metadata of the control sets for B6 - 
controlmut <- rowMeans(y[,c(8,16)])

#calculating fold change between wildtype and controls, +1 to omit any N/A and to allow for plotting on PCA
WTfoldchange <- log2(WT+1) - log2(controlWT+1)
Mutfoldchange <- log2(mutant+1) - log2(controlmut+1)

#combine foldchanges
foldchange <-cbind(WTfoldchange, Mutfoldchange)

#stripping out any N/A values, converting to data matrix and data frame
foldchange1 <- na.omit(foldchange)
yyy <- data.matrix(foldchange1)
fcdat <- data.frame(yyy)

#form PCA plot 
res.pca <- prcomp(t(fcdata[-1:-2]))
get_eig(res.pca)
fviz_eig(res.pca)
fviz_pca_ind(res.pca)

# apply quantile normalization and plot PCA
norm_data = normalize.quantiles(as.matrix(fcdat[-1:-2]),copy=TRUE)
res_norm.pca <- prcomp(t(norm_data))
fviz_eig(res_norm.pca, choice = "eigenvalue", addlabels = TRUE)
fviz_eig(res_norm.pca , geom = "line")
fviz_pca_ind(res_norm.pca)
