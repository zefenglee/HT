rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
miRNet1 <- read.csv('04.miRNet/miRNet-mir-gene-hsa-mirtarbase.csv')
sum(str_detect(miRNet1$mir_id,'hsa-mir'))   
miRNet1$mir_id[str_detect(miRNet1$mir_id,'hsa-mir')] = sub(pattern = 'hsa-mir',replacement = 'hsa-miR',miRNet1$mir_id[str_detect(miRNet1$mir_id,'hsa-mir')])
miRNet1 <- miRNet1[,c(2,4)]
colnames(miRNet1) <- c('Gene1','Gene2')          
sum(duplicated(paste0(miRNet1$Gene1,'_',miRNet1$Gene2)))
write.table(miRNet1,'12.ALL/04.miRNet1.txt',row.names = F,quote = F)

miRNet2 <- read.csv('04.miRNet/miRNet-mir-gene-hsa-tarbase.csv')
sum(str_detect(miRNet2$mir_id,'hsa-mir'))   
miRNet2$mir_id[str_detect(miRNet2$mir_id,'hsa-mir')] = sub(pattern = 'hsa-mir',replacement = 'hsa-miR',miRNet2$mir_id[str_detect(miRNet2$mir_id,'hsa-mir')])
miRNet2 <- miRNet2[,c(2,4)]
colnames(miRNet2) <- c('Gene1','Gene2')          
sum(duplicated(paste0(miRNet2$Gene1,'_',miRNet2$Gene2)))
miRNet2 <- miRNet2[!duplicated(paste0(miRNet2$Gene1,'_',miRNet2$Gene2)),]
write.table(miRNet1,'12.ALL/04.miRNet1.txt',row.names = F,quote = F)

miRNet = rbind(miRNet1,miRNet2)
sum(str_detect(miRNet$Gene1,'hsa-mir')) 
sum(str_detect(miRNet$Gene2,'hsa-mir')) 
sum(duplicated(paste0(miRNet$Gene1,'_',miRNet$Gene2)))
miRNet <- miRNet[!duplicated(paste0(miRNet$Gene1,'_',miRNet$Gene2)),]
write.table(miRNet,'12.ALL/04.miRNet.txt',row.names = F,quote = F)
