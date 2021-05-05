#DESEQ2 Analysis for all samples replicates is drug perturbation data 
#Performed to obtain the log2fc
#15th of March 2021
#Emily Wood

#DESEQ ANALYSIS FOR FOLD CHANGES FOR EACH GROUP WITH RESPECT TO CONTROLS S1-S3(UNTREAT)

#read the count data set gene names to ID, need to make gene names unique as it cannot use duplicate names
cts <- read.delim("readcount_genename.xls")
rownames(cts) <- make.names(cts[,"gene_name"], unique = TRUE)
cts <- as.matrix(cts[,2:40])


#read in the column/sample info 
coldata <- read.csv("Sample details RNA seq experiment.csv", row.names = "Sample.name.in.report")
coldata <- coldata[,c("Sample.ID", "Group.name", "Sample.Detail")]
coldata$samples <- c("untreat", "untreat", "untreat", "T_1", "T_2", "T_3", "M_1", "M_2", "M_3", "TM_1", "TM_2", "TM_3", "TMR_1","TMR_2","TMR_3","TMI_1","TMI_2","TMI_3", "TMRI_1","TMRI_2","TMRI_3","TM_b_1","TM_b_2","TM_b_3", "TMR_b_1","TMR_b_2","TMR_b_3", "TMS_1", "TMS_2", "TMS_3", "TMH_1", "TMH_2", "TMH_3","TMRS_1", "TMRS_2", "TMRS_3", "TMRH_1", "TMRH_2", "TMRH_3" )

#check that col data and count data are in the same order -required for DeSeq
coldata$samples <- factor(coldata$samples)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))


#begin differential expression from the count matrix and column data 
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ samples)
dds

#additional feature data,  added to the DESeqDataSet by adding to the metadata columns
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

#pre-filter low count genes before running the DESeq2 functions
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


#Tell the DESeq2 functions which level you want to compare against S1-S3 UNTREAT
dds$samples <- relevel(dds$samples, ref = "untreat")

#differential expression analysis steps results, which extracts a results table with log2 fold changes, p values and adjusted p values
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)


#write each results into a table
samples_M_1_vs_untreat <- results(dds, contrast = list("samples_M_1_vs_untreat"))
samples_M_2_vs_untreat <- results(dds, contrast = list("samples_M_2_vs_untreat"))
samples_M_3_vs_untreat <- results(dds, contrast = list("samples_M_3_vs_untreat"))
samples_T_1_vs_untreat<- results(dds, contrast = list("samples_T_1_vs_untreat"))
samples_T_2_vs_untreat <- results(dds, contrast = list("samples_T_2_vs_untreat"))
samples_T_3_vs_untreat <- results(dds, contrast = list("samples_T_3_vs_untreat"))
samples_TM_1_vs_untreat <- results(dds, contrast = list("samples_TM_1_vs_untreat"))
samples_TM_2_vs_untreat <- results(dds, contrast = list("samples_TM_2_vs_untreat"))
samples_TM_3_vs_untreat <- results(dds, contrast = list("samples_TM_3_vs_untreat"))
samples_TM_b_1_vs_untreat <- results(dds, contrast = list("samples_TM_b_1_vs_untreat"))
samples_TM_b_2_vs_untreat <- results(dds, contrast = list("samples_TM_b_2_vs_untreat"))
samples_TM_b_3_vs_untreat <- results(dds, contrast = list("samples_TM_b_3_vs_untreat"))
samples_TMH_1_vs_untreat <- results(dds, contrast = list("samples_TMH_1_vs_untreat"))
samples_TMH_2_vs_untreat <- results(dds, contrast = list("samples_TMH_2_vs_untreat"))
samples_TMH_3_vs_untreat <- results(dds, contrast = list("samples_TMH_3_vs_untreat"))
samples_TMI_1_vs_untreat <- results(dds, contrast = list("samples_TMI_1_vs_untreat"))
samples_TMI_2_vs_untreat <- results(dds, contrast = list("samples_TMI_2_vs_untreat"))
samples_TMI_3_vs_untreat <- results(dds, contrast = list("samples_TMI_3_vs_untreat"))
samples_TMR_1_vs_untreat <- results(dds, contrast = list("samples_TMR_1_vs_untreat"))
samples_TMR_2_vs_untreat <- results(dds, contrast = list("samples_TMR_2_vs_untreat"))
samples_TMR_3_vs_untreat <- results(dds, contrast = list("samples_TMR_3_vs_untreat"))
samples_TMR_b_1_vs_untreat <- results(dds, contrast = list("samples_TMR_b_1_vs_untreat"))
samples_TMR_b_2_vs_untreat <- results(dds, contrast = list("samples_TMR_b_2_vs_untreat"))
samples_TMR_b_3_vs_untreat <- results(dds, contrast = list("samples_TMR_b_3_vs_untreat"))
samples_TMRH_1_vs_untreat <- results(dds, contrast = list("samples_TMRH_1_vs_untreat"))
samples_TMRH_2_vs_untreat <- results(dds, contrast = list("samples_TMRH_2_vs_untreat"))
samples_TMRH_3_vs_untreat <- results(dds, contrast = list("samples_TMRH_3_vs_untreat"))
samples_TMRI_1_vs_untreat <- results(dds, contrast = list("samples_TMRI_1_vs_untreat"))
samples_TMRI_2_vs_untreat <- results(dds, contrast = list("samples_TMRI_2_vs_untreat"))
samples_TMRI_3_vs_untreat <- results(dds, contrast = list("samples_TMRI_3_vs_untreat"))
samples_TMRS_1_vs_untreat <- results(dds, contrast = list("samples_TMRS_1_vs_untreat"))
samples_TMRS_2_vs_untreat <- results(dds, contrast = list("samples_TMRS_2_vs_untreat"))
samples_TMRS_3_vs_untreat <- results(dds, contrast = list("samples_TMRS_3_vs_untreat"))
samples_TMS_1_vs_untreat <- results(dds, contrast = list("samples_TMS_1_vs_untreat"))
samples_TMS_2_vs_untreat <- results(dds, contrast = list("samples_TMS_2_vs_untreat"))
samples_TMS_3_vs_untreat <- results(dds, contrast = list("samples_TMS_3_vs_untreat"))



#write all results into a csv files
write.csv(as.data.frame(samples_M_1_vs_untreat), 
          file="samples_M_1_vs_untreat.csv")

write.csv(as.data.frame(samples_M_2_vs_untreat), 
          file="samples_M_2_vs_untreat.csv")

write.csv(as.data.frame(samples_M_3_vs_untreat), 
          file="samples_M_3_vs_untreat.csv")

write.csv(as.data.frame(samples_T_1_vs_untreat), 
          file="samples_T_1_vs_untreat.csv")

write.csv(as.data.frame(samples_T_2_vs_untreat), 
          file="samples_T_2_vs_untreat.csv")

write.csv(as.data.frame(samples_T_3_vs_untreat), 
          file="samples_T_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TM_1_vs_untreat), 
          file="samples_TM_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TM_2_vs_untreat), 
          file="samples_TM_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TM_3_vs_untreat), 
          file="samples_TM_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TM_b_1_vs_untreat), 
          file="samples_TM_b_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TM_b_2_vs_untreat), 
          file="samples_TM_b_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TM_b_3_vs_untreat), 
          file="samples_TM_b_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMH_1_vs_untreat), 
          file="samples_TMH_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMH_2_vs_untreat), 
          file="samples_TMH_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMH_3_vs_untreat), 
          file="samples_TMH_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMI_1_vs_untreat), 
          file="samples_TMI_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMI_2_vs_untreat), 
          file="samples_TMI_2_vs_untreat")

write.csv(as.data.frame(samples_TMI_3_vs_untreat), 
          file="samples_TMI_3_vs_untreat")

write.csv(as.data.frame(samples_TMR_1_vs_untreat), 
          file="samples_TMR_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMR_2_vs_untreat), 
          file="samples_TMR_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMR_3_vs_untreat), 
          file="samples_TMR_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMR_b_1_vs_untreat), 
          file="samples_TMR_b_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMR_b_2_vs_untreat), 
          file="samples_TMR_b_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMR_b_3_vs_untreat), 
          file="samples_TMR_b_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRH_1_vs_untreat), 
          file="samples_TMRH_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRH_2_vs_untreat), 
          file="samples_TMRH_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRH_3_vs_untreat), 
          file="samples_TMRH_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRI_1_vs_untreat), 
          file="samples_TMRI_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRI_2_vs_untreat), 
          file="samples_TMRI_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRI_3_vs_untreat), 
          file="samples_TMRI_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRS_1_vs_untreat), 
          file="samples_TMRS_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRS_2_vs_untreat), 
          file="samples_TMRS_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMRS_3_vs_untreat), 
          file="samples_TMRS_3_vs_untreat.csv")

write.csv(as.data.frame(samples_TMS_1_vs_untreat), 
          file="samples_TMS_1_vs_untreat.csv")

write.csv(as.data.frame(samples_TMS_2_vs_untreat), 
          file="samples_TMS_2_vs_untreat.csv")

write.csv(as.data.frame(samples_TMS_3_vs_untreat), 
          file="samples_TMS_3_vs_untreat.csv")
