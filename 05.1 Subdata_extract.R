rm(list = ls())
library(openxlsx)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(GGally)
library(stringr)
biomarkers_list <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/Intersect_2p_1/features_selected_from_factor.csv')
biomarkers_list <- list('mRNA' = biomarkers_list$feature[biomarkers_list$view=='mRNA']
                        ,'miRNA' = biomarkers_list$feature[biomarkers_list$view=='miRNA']
                        ,'lncRNA' = biomarkers_list$feature[biomarkers_list$view=='lncRNA']
                        ,'circRNA' = biomarkers_list$feature[biomarkers_list$view=='circRNA'])
load_data <- function(){
  mRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'mrna_fpkm',rowNames =T))
  columns <- colnames(mRNA_expr)  
  miRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'mirna_fpkm',rowNames = T)[,columns])
  # rownames(miRNA_expr) <- str_replace_all(pattern = '-',replacement = '_',rownames(miRNA_expr))
  lncRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'lncrna_fpkm',rowNames = T)[,columns])
  # rownames(lncRNA_expr) <- str_replace_all(pattern = '-',replacement = '_',rownames(lncRNA_expr))
  circRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'circrna_fpkm_filtered',rowNames = T)[,columns])
  datas_list <- list('mRNA'=mRNA_expr,'miRNA'=miRNA_expr,'lncRNA'=lncRNA_expr,'circRNA'=circRNA_expr)
  return(datas_list)
}
datas_list <- load_data()
Group <- substr(colnames(datas_list$mRNA),1,1)
Group
fpkm_extract <- function(datas_list,biomarkers_list){
  mRNA_fpkm <- datas_list$mRNA[biomarkers_list$mRNA,]
  miRNA_fpkm <- datas_list$miRNA[biomarkers_list$miRNA,]
  lncRNA_fpkm <- datas_list$lncRNA[biomarkers_list$lncRNA,]
  circRNA_fpkm <- datas_list$circRNA[biomarkers_list$circRNA,]
  fpkm_list <- list('mRNA'=mRNA_fpkm,'miRNA'=miRNA_fpkm,'lncRNA'=lncRNA_fpkm,'circRNA'=circRNA_fpkm)
  return(fpkm_list)
}
fpkm_list <- fpkm_extract(datas_list,biomarkers_list)
fpkm <- as.data.frame(t(rbind(fpkm_list$mRNA,fpkm_list$miRNA,fpkm_list$lncRNA,fpkm_list$circRNA)))
fpkm$Group <- Group
network <- read.table('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/Intersect_2p/network.txt',header = T)
df <- data.frame('NuLL'=c(1:dim(fpkm)[1]),row.names = rownames(fpkm))
for(i in seq(dim(network)[1])){
  gene1 <- network[i,]$Gene1
  gene2 <- network[i,]$Gene2
  colname_multi <- paste0(gene1,'*',gene2)
  colname_divi <- paste0(gene1,'/',gene2)
  df[,colname_divi] <- (fpkm[,gene1])/(fpkm[,gene2])
  df[,colname_multi] <- (fpkm[,gene1])*(fpkm[,gene2])
  
}
df <- df[,-1]

write.table(fpkm,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/subdata/Intersect_2p/subdata_fpkm.csv',sep=',')
write.table(df,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/subdata/Intersect_2p/subdata_fpkm_derived.csv',sep=',')
