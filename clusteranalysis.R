''''
##phenotypic contribtiion of each clusters
#Author:EMilt Wood
#5th of May 2021
''''

#read in the files
file <- read.csv("total_dataset_com.csv", header = TRUE)
file <- file[,c(1,2,9,10,11,12,13,14,15,16,17,18,19,20)]

#read in the file with the clusters list
subset <- read.csv("Clusters.csv")
  
 #multiple the FC for each gene in the cluster by the STV 
file$T <- file$norm_vec * file$T_log2FC
file$TM <- file$norm_vec * file$TM_log2FC
file$M <- file$norm_vec * file$M_log2FC
file$TMR <- file$norm_vec * file$TMR_log2FC
file$TMI <- file$norm_vec * file$TMI_log2FC
file$TMRI <- file$norm_vec * file$TMRI_log2FC
file$TM_b <- file$norm_vec * file$TM_b_log2FC
file$TMR_b <- file$norm_vec * file$TMR_b_log2FC
file$TMS <- file$norm_vec * file$TMS_log2FC
file$TMH <- file$norm_vec * file$TMH_log2FC
file$TMRS <- file$norm_vec * file$TMRS_log2FC
file$TMRH <- file$norm_vec * file$TMRH_log2FC

#new file with all the caclulated contributions
file <- file[,c(1,2,15,16,17,18,19,20,21,22,23,24,25,26)]


#subset the genes from each cluster into new file
C1<- file[file$gene %in% subset$Cluster.1,]
C2<- file[file$gene %in% subset$Cluster.2,]
C3<- file[file$gene %in% subset$Cluster.3,]
C4<- file[file$gene %in% subset$Cluster.4,]
C5<- file[file$gene %in% subset$Cluster.5,]
C6<- file[file$gene %in% subset$Cluster.6,]
C7<- file[file$gene %in% subset$Cluster.7,]
C8<- file[file$gene %in% subset$Cluster.8,]
C9<- file[file$gene %in% subset$Cluster.9,]
C10<- file[file$gene %in% subset$Cluster.10,]



#Sum each product over the cluster
sum(C10$T)
sum(C10$M)
sum(C10$TM)
sum(C10$TMR)
sum(C10$TMI)
sum(C10$TMRI)
sum(C10$TM_b)
sum(C10$TMR_b)
sum(C10$TMS)
sum(C10$TMH)
sum(C10$TMRS)
sum(C10$TMRH)

