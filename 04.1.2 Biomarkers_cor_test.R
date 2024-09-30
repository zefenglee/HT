rm(list = ls())
library(openxlsx)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(GGally)
library(stringr)
# biomarkers_list <- read.delim('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/Intersect_2p/biomarker_included_in_Interactome.txt')
biomarkers_list <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/Intersect_2p_1/features_selected_from_factor.csv')
biomarkers_list <- list('mRNA' = biomarkers_list$feature[biomarkers_list$view=='mRNA']
                        ,'miRNA' = str_replace_all(pattern = '-',replacement = '_',biomarkers_list$feature[biomarkers_list$view=='miRNA'])
                        ,'lncRNA' = str_replace_all(pattern = '-',replacement = '_',biomarkers_list$feature[biomarkers_list$view=='lncRNA'])
                        ,'circRNA' = biomarkers_list$feature[biomarkers_list$view=='circRNA'])
load_data <- function(){
  mRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'mrna_fpkm',rowNames =T))
  columns <- colnames(mRNA_expr)  
  miRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'mirna_fpkm',rowNames = T)[,columns])
  rownames(miRNA_expr) <- str_replace_all(pattern = '-',replacement = '_',rownames(miRNA_expr))
  lncRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'lncrna_fpkm',rowNames = T)[,columns])
  rownames(lncRNA_expr) <- str_replace_all(pattern = '-',replacement = '_',rownames(lncRNA_expr))
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
DEA_m <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/mRNA/Intersection_2p_mRNA_1.csv')
DEA_mi <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/miRNA/Intersection_2p_miRNA_1.csv')
DEA_mi$Id <- str_replace_all(pattern = '-',replacement = '_',DEA_mi$Id)
DEA_lnc <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/lncRNA/Intersection_2p_lncRNA_1.csv')
DEA_lnc$Id <- str_replace_all(pattern = '-',replacement = '_',DEA_lnc$Id)
DEA_circ <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/circRNA/Intersection_2p_circRNA_1.csv')
DEAs <- rbind(DEA_m[DEA_m$Id %in% biomarkers_list$mRNA,]
      ,DEA_mi[DEA_mi$Id %in% biomarkers_list$miRNA,]
      ,DEA_lnc[DEA_lnc$Id %in% biomarkers_list$lncRNA,]
      ,DEA_circ[DEA_circ$Id %in% biomarkers_list$circRNA,])
rownames(DEAs) <- DEAs$Id

gene1 <- c()
gene2 <- c()
expr1 <- t(fpkm_list$circRNA)
biomarkers1 <-colnames(expr1)
expr2 <- t(fpkm_list$miRNA)
biomarkers2 <- colnames(expr2)
expr <- as.data.frame(cbind(expr1,expr2))

for(i in biomarkers1){
  for(j in biomarkers2){
    t <- cor.test(expr[,i],expr[,j])
    p <- t$p.value
    cor_coef <- round(as.numeric(t$estimate),2)
    if(p<0.05 & abs(cor_coef)>=0.2){
      gene1 <- c(gene1,i)
      gene2 <- c(gene2,j)
    }
  }
}
cor_expression_network <- data.frame('Gene1'=gene1
                                     ,'Gene2'=gene2)
cor_expression_network$Gene1_Change <- 'Up'
cor_expression_network$Gene2_Change <- 'Up'
for(i in seq(length(gene1))){
  cor_expression_network$Gene1_Change[cor_expression_network$Gene1==gene1[i]]=DEAs[gene1[i],'Change']
  cor_expression_network$Gene2_Change[cor_expression_network$Gene2==gene2[i]]=DEAs[gene2[i],'Change']
}
cor_expression_network$Gene2 <- str_replace_all(pattern = '_',replacement = '-',cor_expression_network$Gene2)
write.table(cor_expression_network,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/06.Co_expression/circ_mi/Co_expression_network_circ_mi.txt',,row.names = F,quote = F,sep='\t')

# expr$Group <- Group
# expr$Group <- as.factor(expr$Group)
df <- matrix(nrow = length(biomarkers1),ncol = length(biomarkers2),0,dimnames = list(biomarkers1,biomarkers2))

for(i in biomarkers1){
  for(j in biomarkers2){
    t <- cor.test(expr[,i],expr[,j])
    p <- t$p.value
    cor_coef <- round(as.numeric(t$estimate),2)
    if(p<0.05){
      cor_coef <- paste0(cor_coef,'*')
    }
    df[i,j] <- cor_coef
  }
}
write.table(df,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/06.Co_expression/circ_mi/circRNA_miRNA_cor_test.csv',sep = ',',row.names = T)
