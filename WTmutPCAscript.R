#This is the PCA plot of the RNA_seq_TNF_full_new_Nov2018 data without being normalised to controls
#packages for PCA
library(gplots) 
library(ggfortify)
library(preprocessCore)
library(factoextra)


#load in the data
data <- read.delim(file = "RNA_seq_TNF_full_new_Nov2018.txt", header = TRUE, row.names = 1)
TPMdata <- data[-1]


#stripping out any N/A values
TPMdata <- na.omit(TPMdata)

#convert data to matrix
data.matrix(TPMdata)

res.pca <- prcomp(t(TPMdata[-1:-2]))
get_eig(res.pca)
fviz_eig(res.pca)
fviz_pca_ind(res.pca)

# apply quantile normalization and plot PCA
norm_data = normalize.quantiles(as.matrix(TPMdata[-1:-2]),copy=TRUE)
res_norm.pca <- prcomp(t(norm_data))
fviz_eig(res_norm.pca)
fviz_pca_ind(res_norm.pca)
