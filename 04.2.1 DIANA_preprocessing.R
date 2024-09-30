library(openxlsx)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
DiANA <- read.table('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/01.DIANA/DIANA_miRNA_mRNA.txt',header = T)
head(DiANA)
DiANA <- DiANA[DiANA$interaction_score>=quantile(DiANA$interaction_score,0.75),]
id_trans <- bitr(DiANA$ensembl_gene_id
                 , fromType = 'ENSEMBL'
                 ,toType = 'SYMBOL'
                 ,OrgDb = 'org.Hs.eg.db')
colnames(DiANA)
colnames(id_trans) <- c('ensembl_gene_id','Gene2')
DiANA <- inner_join(DiANA,id_trans,by='ensembl_gene_id')
sum(is.na(DiANA$Gene2))
sum(is.null(DiANA$Gene2))
DiANA <- DiANA[,c("mirna","Gene2")]
colnames(DiANA) <- c("Gene1","Gene2")
write.table(DiANA,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/12.ALL/01.DIANA_miRNA_mRNA.txt',row.names = F,quote = FALSE)
rm(list=ls())
DiANA <- read.table('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/01.DIANA/tarbase_data.csv',sep='\t',header = T)
head(DiANA)
DiANA <- DiANA[(DiANA$species=='Homo sapiens')&(DiANA$direct_indirect=='DIRECT'),]
colnames(DiANA)
DiANA <- DiANA[,c("mirna","geneName")]
colnames(DiANA) <- c("Gene1","Gene2")
sum(str_detect(DiANA$Gene1,'hsa-mir'))

DiANA_a <- read.table('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/12.ALL/01.DIANA_miRNA_mRNA.txt',header = T)

DiANA <- rbind(DiANA_a,DiANA)

sum(DiANA$Gene1==DiANA$Gene2)

DiANA <- DiANA[!duplicated(paste0(DiANA$Gene1,DiANA$Gene2)),]

sum(!duplicated(c(DiANA$Gene1,DiANA$Gene2)))

write.table(DiANA,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/12.ALL/01.DIANA.txt',row.names = F,quote = FALSE)
