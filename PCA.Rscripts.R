'''
#PLOTTING PCA on R studio normalised for controls
#Feburary 2021
#Emily Wood
'''

#packages for heatmap
library(gplots) 
library(pheatmap)
library(heatmaply)


#load in the data
f<- read.delim(file = "RNA_seq_TNF_full_new_Nov2018.txt", header = TRUE, row.names = 1)

#selecting genes of interest and renaming row to gene names
data <- f[c("ENSMUSG00000027187","ENSMUSG00000024661", "ENSMUSG00000050708", "ENSMUSG00000031400", "ENSMUSG00000075706", "ENSMUSG00000058135", "ENSMUSG00000015839", "ENSMUSG00000038615","ENSMUSG00000050114","ENSMUSG00000070034","ENSMUSG00000070031","ENSMUSG00000032802","ENSMUSG00000005354","ENSMUSG00000020250") ,]
TPMdata <- data.frame(data, row.names = "external_gene_name")


#stripping out any N/A values
TPMdata <- na.omit(TPMdata)

#convert data to matrix
data.matrix(TPMdata)

res.pca <- prcomp(t(TPMdata[-1:-2]))
get_eig(res.pca)
fviz_eig(res.pca)
fviz_pca_ind(res.pca)

# apply quantile normalization and plot PCA
library(preprocesscore)
norm_data = normalize.quantiles(as.matrix(TPMdata[-1:-2]),copy=TRUE)
res_norm.pca <- prcomp(t(norm_data))
fviz_eig(res_norm.pca)
fviz_pca_ind(res_norm.pca)


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
