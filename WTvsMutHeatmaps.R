#PLOTTING Heatmaps on R studio for Wild types
#Feburary 2021
#Emily Wood
#packages for heatmap
library(gplots) 
library(pheatmap)
library(heatmaply)


#load in the data
f<- read.delim(file = "RNA_seq_TNF_full_new_Nov2018.txt", header = TRUE, row.names = 1)

#selecting genes of interest and renaming row to gene names
data <- f[c("ENSMUSG00000027187","ENSMUSG00000024661", "ENSMUSG00000050708", "ENSMUSG00000031400", "ENSMUSG00000075706", "ENSMUSG00000058135", "ENSMUSG00000015839", "ENSMUSG00000038615","ENSMUSG00000050114","ENSMUSG00000070034","ENSMUSG00000070031","ENSMUSG00000032802","ENSMUSG00000005354","ENSMUSG00000020250") ,]
data1 <- data.frame(data, row.names = "external_gene_name")


#convert to numerical data
y <- data.matrix(data1)

#convert to log scale
log_y <- log10(y+1)
sc_L <- scale(log_y, scale=TRUE)

#generating heatmap, scaling to row (each sample)
pheatmap(sc_L, main = "Gene expression profiles" , scale = "row")

