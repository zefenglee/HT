library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(tidyverse)
library(openxlsx)
rm(list = ls())
gene_list <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/Intersect_2p_1/features_selected_from_factor.csv')
gene_list <- gene_list$Symbol[1:3]
gene_list
load('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/edgeR/biomarkers_targetgenes.RData')
Enrichment_analysis <- function(gene_list,keyType){
  ego_CC <- data.frame(enrichGO(gene=gene_list, OrgDb= 'org.Hs.eg.db', keyType=keyType, ont="CC", pvalueCutoff= 1,qvalueCutoff= 1))
  ego_MF <- data.frame(enrichGO(gene=gene_list, OrgDb= 'org.Hs.eg.db', keyType=keyType, ont="MF", pvalueCutoff= 1,qvalueCutoff= 1))
  ego_BP <- data.frame(enrichGO(gene=gene_list, OrgDb= 'org.Hs.eg.db', keyType=keyType, ont="BP", pvalueCutoff= 1,qvalueCutoff= 1))
  ego_CC$resource="CC"
  ego_MF$resource="MF"
  ego_BP$resource="BP"
  id_trans <- bitr(gene = gene_list,fromType=keyType, toType="ENTREZID", OrgDb='org.Hs.eg.db')
  id_trans_ <- pull(id_trans,ENTREZID)
  ekegg <- enrichKEGG(gene = id_trans_, organism = 'hsa', pvalueCutoff = 1, qvalueCutoff = 1)
  ekegg <- data.frame(ekegg)
  ekegg$resource='KEGG'
  a = c()
  for(i in ekegg$geneID){
    geneids <- str_split(i,'/')[[1]]
    geneids_<- bitr(gene = geneids,fromType="ENTREZID", toType=keyType, OrgDb='org.Hs.eg.db')
    geneids_ <- pull(geneids_,keyType)
    #print(geneids_)
    b <- paste0(geneids_,'/',collapse = '')
    a<-c(a,b)
  }
  ekegg$geneID=a
  pathways <-rbind(ego_BP,ego_CC,ego_MF,ekegg)
  # write.table(pathways,paste0(filepath,'_pathways.csv'),sep=',',row.names = F)
  return(pathways)
}
df <- Enrichment_analysis(gene_list = gene_list ,keyType = 'SYMBOL')
mRNA_pathways <- Enrichment_analysis(gene_list = mRNAs,keyType = 'SYMBOL',filepath = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/05.Pathways/edgeR/mRNA')
miRNA_targets_pathways <- Enrichment_analysis(gene_list = miRNA_targets_validated,keyType = 'SYMBOL',filepath = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/05.Pathways/edgeR/miRNA_targets')
lncRNA_targets_pathways <- Enrichment_analysis(gene_list = lncRNA_targets,keyType = 'SYMBOL',filepath = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/05.Pathways/edgeR/lncRNA_targets')
