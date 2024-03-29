rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
miRtarbase <- read.xlsx('03.miRTarBase/03.miRTarBase_SE_WR.xlsx',sheet = 'miRTarBase_SE_WR')
miRtarbase <- miRtarbase[miRtarbase$`Species.(miRNA)`=='Homo sapiens',]
table(miRtarbase$`Species.(miRNA)`)
sum(str_detect(miRtarbase$miRNA,'hsa-mir'))   
colnames(miRtarbase)
miRtarbase <- miRtarbase[,c("miRNA","Target.Gene")]
colnames(miRtarbase) <- c('Gene1','Gene2')
write.table(miRtarbase,'12.ALL/03.miRTarBase.txt',row.names = F,quote = F)
sum(!duplicated(c(miRtarbase$Gene1,miRtarbase$Gene2)))
