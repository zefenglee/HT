rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
NPinter1 <- read.table('05.NPinter4/circRNA_interaction.txt',sep = '\t',header = T)
NPinter1 <- NPinter1[(NPinter1$organism=='Homo sapiens')&(NPinter1$ncID!='-'),c("ncID","tarName" )]
NPinter1$ncID <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',NPinter1$ncID )
NPinter1$tarName <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',NPinter1$tarName )
colnames(NPinter1) <- c('Gene1','Gene2')

NPinter2 <- read.xlsx('05.NPinter4/lncRNA_interaction.xlsx',sheet = 'lncRNA_interaction')
NPinter2 <- NPinter2[(NPinter2$organism=='Homo sapiens')&(NPinter2$experiment!='-'),c("ncName","tarName" )]
sum(str_detect(NPinter2$tarName,'hsa-mir'))
NPinter2$tarName <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',NPinter2$tarName)
sum(str_detect(NPinter2$ncName,'hsa-mir'))
NPinter2$ncName <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',NPinter2$ncName)
NPinter2 <- NPinter2[!duplicated(paste0(NPinter2$ncName,NPinter2$tarName)),]
colnames(NPinter2) <- c('Gene1','Gene2')

NPinter3 <- read.xlsx('05.NPinter4/miRNA_interaction.xlsx',sheet = 'miRNA_interaction')
NPinter3 <- NPinter3[(NPinter3$organism=='Homo sapiens')&(NPinter3$experiment!='-'),c("ncName","tarName" )]
sum(str_detect(NPinter3$tarName,'hsa-mir'))
NPinter3$tarName <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',NPinter3$tarName)
sum(str_detect(NPinter3$ncName,'hsa-mir'))
NPinter3$ncName <- sub(pattern = 'hsa-mir',replacement = 'hsa-miR',NPinter3$ncName)
NPinter3 <- NPinter3[!duplicated(paste0(NPinter3$ncName,NPinter3$tarName)),]
colnames(NPinter3) <- c('Gene1','Gene2')

NPinter <- rbind(NPinter1,NPinter2,NPinter3)

sum(str_detect(NPinter$Gene1,'hsa-mir'))
sum(str_detect(NPinter$Gene2,'hsa-mir'))

NPinter <- NPinter[!duplicated(paste0(NPinter$Gene1,NPinter$Gene2)),]

sum(!duplicated(c(NPinter$Gene1,NPinter$Gene2)))    
write.table(NPinter,'12.ALL/05.NPinter.txt',row.names = F,quote = F,sep='\t')
