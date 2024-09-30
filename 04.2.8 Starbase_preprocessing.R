rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
Starbase1 <- read.delim("F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/08.Starbase/Human/circRNA_miRNA_interaction.txt")
colnames(Starbase1)
Starbase1 <- Starbase1[,c("miRNAname", "geneName")]
colnames(Starbase1) <- c('Gene1','Gene2')

Starbase2 <- read.delim("F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/08.Starbase/Human/lncRNA_miRNA_interaction.txt")
colnames(Starbase2)
Starbase2 <- Starbase2[,c("miRNAname", "geneName")]
colnames(Starbase2) <- c('Gene1','Gene2')

Starbase3 <- read.delim("F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/08.Starbase/Human/mRNA_miRNA_interaction.txt")
colnames(Starbase3)
Starbase3 <- Starbase3[,c("miRNAname", "geneName")]
colnames(Starbase3) <- c('Gene1','Gene2')

Starbase <- rbind(Starbase1,Starbase2,Starbase3)

sum(duplicated(paste0(Starbase$Gene1,Starbase$Gene2)))
Starbase <- Starbase[!duplicated(paste0(Starbase$Gene1,Starbase$Gene2)),]
sum(!duplicated(c(Starbase$Gene1,Starbase$Gene2)))
write.table(Starbase,'12.ALL/08.Starbase.txt',row.names = F,quote = F)
