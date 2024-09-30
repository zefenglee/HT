library(stringr)
library(limma)
library(openxlsx)
library(tidyverse)
library(edgeR)
rm(list = ls())
limma_DEA <- function(sheet,Log2FC_cutoff,file_path){
  print(paste0('data type: ',sheet))
  data <- read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = sheet,rowNames = T)
  print(paste0('data dim(unfiltered): ',dim(data)[1],'*',dim(data)[2]))
  group_list <- ifelse(str_detect(colnames(data),'T'),'T','C')
  group_list <-factor(group_list,levels = c('T','C'))
  group_list
  
  if((sheet == 'circrna_count')|(sheet == 'mirna_count')){
    data <- data[rowSums(data==0)<=0.5*length(colnames(data)),]
  }else{
    data <- data[rowSums(data==0)<=0.75*length(colnames(data)),]
    data <- data[rowMeans(data)>10,]
    }
  print(paste0('data dim(filtered): ',dim(data)[1],'*',dim(data)[2]))
  #差异分析
  design<-model.matrix(~0+group_list)
  colnames(design)=levels(group_list)
  rownames(design)=colnames(data)
  dge<-DGEList(counts = data)
  dge<-calcNormFactors(dge)
  v<-voom(dge,design )
  fit<-lmFit(v,design)
  contrast=paste(levels(group_list),collapse = '-')
  contrast
  cont.matrix<-makeContrasts(contrasts = contrast,levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  DEG<-topTable(fit2,coef = contrast,n=Inf)
  colnames(DEG)
  colnames(DEG)=c("Log2FC","AveExpr","t","P_Value","adj_P_Value" ,"B")
  if(sheet=='circrna_count'){
    anno <- read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'circrna_anno',rowNames = F)
    DEG$circRNA_ID=rownames(DEG)
    DEG <- inner_join(DEG,anno,by='circRNA_ID')
    DEG <- na.omit(DEG)
    rownames(DEG) <- DEG$circBase_ID
    DEG <- DEG[,c("Log2FC","AveExpr","t","P_Value","adj_P_Value" ,"B")]
    }
  #定义筛选标准，输出表格
  K_Up<-(DEG$P_Value<0.05)&(DEG$Log2FC>Log2FC_cutoff)
  K_Down<-(DEG$P_Value<0.05)&(DEG$Log2FC< -Log2FC_cutoff)
  DEG$Change='Normal'
  DEG$Change[K_Up]='Up'
  DEG$Change[K_Down]='Down'
  print(paste0('cutoff: ',Log2FC_cutoff))
  save_path = paste0(file_path,'/limma_',sheet,'_','l2fc_',Log2FC_cutoff,'.csv')
  write.table(DEG,save_path,sep=',')
  DEG_col <- data.frame('Id'=rownames(DEG)[DEG$Change!='Normal'],
                        'Change'=DEG$Change[DEG$Change!='Normal'])
  save_col <- paste0(file_path,'/limma_',sheet,'_','l2fc_',Log2FC_cutoff,'_cols.csv')
  write.table(DEG_col,save_col,sep=',')
  print(table(DEG$Change))
  stat <- data.frame(
    c(
      as.numeric(
        ifelse(
          is.na(table(DEG$Change)['Down']),0,table(DEG$Change)['Down'])
      )
      ,as.numeric(
        ifelse(
          is.na(table(DEG$Change)['Up']),0,table(DEG$Change)['Up'])
        )
        )
      )
  colnames(stat) <- paste0(sheet,'_','l2fc_',Log2FC_cutoff)
  res <- list(DEG,stat)
  return(res)
}

DEGs_stats <- function(cutoffs,file_path,sheet){
  stats <- data.frame('NULL'=c(0,0),row.names = c('Down','Up'))
  for(cutoff in cutoffs){
  res <- limma_DEA(sheet = sheet,Log2FC_cutoff = cutoff,file_path = file_path)
  stat <- res[2]
  stats <- cbind(stats,stat) 
  }
  stats <- stats[,-1]
  save_path <- paste0(file_path,'/limma_',sheet,'_','stats','.csv')
  write.table(stats,save_path,sep=',',row.names = T)
  return(stats)
}
cutoffs <- c(0,0.5,1,1.5,2,2.5,3)
file_path = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/limma/mRNA/'
sheet = 'mrna_count'
stats_mRNA <- DEGs_stats(cutoffs = cutoffs ,file_path = file_path,sheet = sheet)

file_path = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/limma/miRNA/'
sheet = 'mirna_count'
stats_mRNA <- DEGs_stats(cutoffs = cutoffs ,file_path = file_path,sheet = sheet)

file_path = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/limma/circRNA/'
sheet = 'circrna_count'
stats_mRNA <- DEGs_stats(cutoffs = cutoffs ,file_path = file_path,sheet = sheet)

file_path = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/limma/lncRNA/'
sheet = 'lncrna_count'
stats_mRNA <- DEGs_stats(cutoffs = cutoffs ,file_path = file_path,sheet = sheet)
