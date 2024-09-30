rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
PPI <- read.csv('11.文献/PPI.csv')[,c(1,2)]
colnames(PPI) <- c('Gene1','Gene2')
sum(duplicated(paste0(PPI$Gene1,PPI$Gene2)))

PPI_NC <- read.csv('11.文献/ncPPI_PPI.csv')[,c(1,2)]
colnames(PPI_NC) <- c('Gene1','Gene2')
sum(duplicated(paste0(PPI_NC$Gene1,PPI_NC$Gene2)))

Interact <- rbind(PPI,PPI_NC)
sum(duplicated(paste0(Interact$Gene1,Interact$Gene2)))
Interact <- Interact[!duplicated(paste0(Interact$Gene1,Interact$Gene2)),]

sum(!duplicated(c(Interact$Gene1,Interact$Gene2)))
dim(Interact)
write.table(Interact,'12.ALL/11.Article.txt',row.names = F,quote = F)
