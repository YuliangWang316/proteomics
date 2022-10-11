if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(tidyverse)
library(DESeq2)
#import data
setwd("E:/TEST3/")
mycounts<-read.table("proteinGroups2.txt",header = TRUE,row.names = 1,sep ="\t" )
condition<-factor(c(rep("AF",10),rep("SR",10)),levels = c("SR","AF"))
colData<-data.frame(row.names = colnames(mycounts),condition)
mycounts_new<-as.data.frame(apply(mycounts,MARGIN = 2,as.integer))
rownames(mycounts_new)<-rownames(mycounts)
mycounts_new<-na.omit(mycounts_new)

dds <- DESeqDataSetFromMatrix(mycounts_new, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="All_results_protein2.csv")

