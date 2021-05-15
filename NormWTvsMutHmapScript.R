'''
#Author: Emily Wood
#Date: 14/03/21
#This is the heatmap for wildtupe 1h versus Mutant B6 with respect to controls - subset of genes.
'''



#packages for heatmap
library(pheatmap)

#load in the data
f <- read.delim(file = "RNA_seq_TNF_full_new_Nov2018.txt", header = TRUE, row.names = 1)

#selecting genes of interest and renaming row to gene names
data <- f[c("ENSMUSG00000027187","ENSMUSG00000024661", "ENSMUSG00000050708", "ENSMUSG00000031400", "ENSMUSG00000075706", "ENSMUSG00000058135", "ENSMUSG00000015839", "ENSMUSG00000038615","ENSMUSG00000050114","ENSMUSG00000070034","ENSMUSG00000070031","ENSMUSG00000032802","ENSMUSG00000005354","ENSMUSG00000020250") ,]
data <- data.frame(data, row.names = "external_gene_name")

#convert to numerical data
y <- data.matrix(data)

#vectors containing metadata 
#mutant dataset - B6
mutant <- y[c(1:14), c(1, 2, 3, 9, 10, 11)]
#Wild-type dataset - 1h
WT <- y[c(1:14), c(5,6,7,13,14,15)]
#metadata of the control sets for 1h - 
controlWT <- rowMeans(y[,c(4,12)])
#metadata of the control sets for B6 - 
controlmut <- rowMeans(y[,c(8,16)])

#calculating fold change between wildtype and controls, +1 to omit any N/A and to allow for plotting on heatmap
WTfoldchange <- log2(WT+1) - log2(controlWT+1)
Mutfoldchange <- log2(mutant+1) - log2(controlmut+1)

#combine foldchanges
foldchange <-cbind(WTfoldchange, Mutfoldchange)

#scaling columns of numeric matrix
data3 <- scale(foldchange)

#generating heatmap, scaling to row (each gene)
pheatmap(data3, main = "Fold change B6(mutant) vs 1h(WT)", scale = "row", show_colnames = "WTfoldchange")

