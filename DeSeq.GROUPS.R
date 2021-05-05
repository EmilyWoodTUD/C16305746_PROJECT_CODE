#DESEQ2 Analysis for all drug groups (3 replicated each) is drug perturbation data 
#Performed to obtain the log2fc
#15th of March 2021
#Emily Wood


#DESEQ ANALYSIS FOR FOLD CHANGES FOR EACH GROUP WITH RESPECT TO CONTROLS S1-S3(UNTREAT)

#read the count data set gene names to ID, need to make gene names unique as it cannot use duplicate names
cts <- read.delim("readcount_genename.xls")
rownames(cts) <- make.names(cts[,"gene_name"], unique = TRUE)
cts <- as.matrix(cts[,2:40])
cts <- cts[,-11]

#read in the column/sample info 
coldata <- read.csv("Sample details RNA seq experiment.csv", row.names = "Sample.name.in.report")
coldata <- coldata[,c("Sample.ID", "Group.name", "Sample.Detail")]
coldata <- coldata[-11,]

#check that col data and count data are in the same order -required for DeSeq
coldata$Group.name <- factor(coldata$Group.name)
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))


#begin differential expression from the count matrix and column data 
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group.name)
dds

#additional feature data,  added to the DESeqDataSet by adding to the metadata columns
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

#pre-filter low count genes before running the DESeq2 functions
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Tell the DESeq2 functions which level you want to compare against S1-S3 UNTREAT
dds$Group.name <- relevel(dds$Group.name, ref = "UNTREAT")

#differential expression analysis steps results, which extracts a results table with log2 fold changes, p values and adjusted p values
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)



#write each results into a table
res1_T_vs_UNTREAT <- results(dds, contrast = list("Group.name_T_vs_UNTREAT"))
res2_M_UNTREAT <- results(dds, contrast = list("Group.name_M_vs_UNTREAT"))
res3_TM_vs_UNTREAT <- results(dds, contrast = list("Group.name_TM_vs_UNTREAT"))
res4_TMR_vs_UNTREAT<- results(dds, contrast = list("Group.name_TMR_vs_UNTREAT"))
res5_TMI_vs_UNTREAT <- results(dds, contrast = list("Group.name_TMI_vs_UNTREAT"))
res6_TMRI_vs_UNTREAT <- results(dds, contrast = list("Group.name_TMRI_vs_UNTREAT"))

res7_TM_b_vs_UNTREAT<- results(dds, contrast = list("Group.name_TM__vs_UNTREAT"))
res8_TMR_b_vsUNTREAT <- results(dds, contrast = list("Group.name_TMR__vs_UNTREAT"))
res9_TMS_vs_UNTREAT <- results(dds, contrast = list("Group.name_TMS_vs_UNTREAT"))
res10_TMH_vs_UNTREAT <- results(dds, contrast = list("Group.name_TMH_vs_UNTREAT"))
res11_TMRS_vs_UNTREAT <- results(dds, contrast = list("Group.name_TMRS_vs_UNTREAT"))
res12_TMRH_vs_UNTREAT <- results(dds, contrast = list("Group.name_TMRH_vs_UNTREAT"))


#write all results into a csv files
write.csv(as.data.frame(res1_T_vs_UNTREAT), 
          file="res1_T_vs_UNTREAT.csv")

write.csv(as.data.frame(res2_M_UNTREAT), 
          file="res2_M_UNTREAT.csv")

write.csv(as.data.frame(res3_TM_vs_UNTREAT), 
          file="res3_TM_vs_UNTREAT.csv")

write.csv(as.data.frame(res4_TMR_vs_UNTREAT), 
          file="res4_TMR_vs_UNTREAT.csv")

write.csv(as.data.frame(res5_TMI_vs_UNTREAT), 
          file="res5_TMI_vs_UNTREAT.csv")

write.csv(as.data.frame(res6_TMRI_vs_UNTREAT ), 
          file="res6_TMRI_vs_UNTREAT .csv")

write.csv(as.data.frame(res7_TM_b_vs_UNTREAT), 
          file="res7_TM_b_vs_UNTREAT.csv")

write.csv(as.data.frame(res8_TMR_b_vsUNTREAT), 
          file="res8_TMR_b_vsUNTREAT.csv")

write.csv(as.data.frame(res9_TMS_vs_UNTREAT), 
          file="res9_TMS_vs_UNTREAT.csv")

write.csv(as.data.frame(res10_TMH_vs_UNTREAT), 
          file="res10_TMH_vs_UNTREAT.csv")

write.csv(as.data.frame(res11_TMRS_vs_UNTREAT), 
          file="res11_TMRS_vs_UNTREAT.csv")

write.csv(as.data.frame(res12_TMRH_vs_UNTREAT), 
          file="res12_TMRH_vs_UNTREAT.csv")