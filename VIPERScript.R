#VIPER ANALYSIS 
#17th April 2021
#Emily Wood

library("viper")
library("mixtools")
library("bcellViper")


#read in the network , convert to a value for arcane2regulon function 
nwfile <- readLines( file("NETWORK_Macrophage_onlyTFs_50.tsv"), n=50000)
afile <- textConnection(nwfile)


#eset 
eset <- read.delim("Macrophage_expression_matrix.tsv")
eset <- as.matrix(eset)


#generating the regulon 
regul <- aracne2regulon(afile = afile, eset = eset, gene = FALSE, format = "3col", verbose = T)
print(regul)


#reading in the expression set used to generate the network file
pd <- read.csv("Log2FC.Samples.csv")
rownames(pd) <- make.names(pd[,"X"], unique = TRUE)
pd <- pd[,2:36]
pd <- as.matrix(pd)


#generating signatures
signature <- rowTtest(pd[,30:32], pd[,1:2])
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) *
                + sign(signature$statistic))[, 1]


#generate the null models
dnull <- ttestNull(pd[,30:32], pd[,1:2], per=1000, repos = TRUE, verbose = T)


#run the viper analysis with viper function
mrs <- msviper(signature, regul, dnull, pleiotropy = FALSE,
               minsize = 25, adaptive.size = FALSE, ges.filter = TRUE,
               synergy = 0, level = 10, pleiotropyArgs = list(regulators = 0.05,
                                                              shadow = 0.05, targets = 10, penalty = 20, method = "adaptive"),
               cores = 1, verbose = TRUE)
summary(mrs)


#plotting data on MSplot
plot(mrs, mrs = 10, color = c("cornflowerblue",
                              "salmon"), pval = NULL, bins = 500, cex = 0, density = 0,
     smooth = 0, sep = 0.2, hybrid = TRUE, include = c("expression",
                                                       "activity"), gama = 2)



#write to files
mrs[["regulon"]]
write.csv(as.data.frame(mrs[["regulon"]]), file="VIPERregulon.csv")
write.csv(as.data.frame(mrs[["es"]][["size"]]), file="VIPERsize.csv")
write.csv(as.data.frame(mrs[["es"]][["nes"]]), file="VIPERnes.csv")
write.csv(as.data.frame(mrs[["es"]][["p.value"]]), file="VIPERpval.csv")
