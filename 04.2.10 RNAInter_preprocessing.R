rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
RNAinter <- read.delim('10.RNAInter/RNAInter.txt')
table(RNAinter$Species1)
RNAinter <- RNAinter[(RNAinter$Species1=='Homo sapiens')&
                       (RNAinter$Species2=='Homo sapiens')&
                       (RNAinter$strong!='N/A'),]
colnames(RNAinter)
RNAinter <- RNAinter[,c("Interactor1.Symbol","Interactor2.Symbol")]
colnames(RNAinter) <- c('Gene1','Gene2')
sum(str_detect(RNAinter$Gene1,'hsa-mir'))
sum(str_detect(RNAinter$Gene2,'hsa-mir'))
RNAinter$Gene1 <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',RNAinter$Gene1)
RNAinter$Gene2 <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',RNAinter$Gene2)
RNAinter <- RNAinter[!duplicated(paste0(RNAinter$Gene1,RNAinter$Gene2)),]
sum(duplicated(paste0(RNAinter$Gene1,RNAinter$Gene2)))
sum(!duplicated(c(RNAinter$Gene1,RNAinter$Gene2)))
write.table(RNAinter,'12.ALL/10.RNAInter.txt',row.names = F,quote = F,sep='\t')
