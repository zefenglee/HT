rm(list = ls())
library(tidyverse)
interact_DEGs <- function(sheet,Log2FC_cutoff){
  limma_DEG <- read.csv(
    paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/limma/',sheet,'/limma_',tolower(sheet),'_count_l2fc_',Log2FC_cutoff,'.csv')
  )
  limma_DEG$Id <- rownames(limma_DEG)
  limma_Up <- limma_DEG[limma_DEG$Change=='Up',c('Id','Change')]
  limma_Down <- limma_DEG[limma_DEG$Change=='Down',c('Id','Change')]
  edgeR_DEG <- read.csv(
    paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/edgeR/',sheet,'/edgeR_',tolower(sheet),'_count_l2fc_',Log2FC_cutoff,'.csv')
  )
  edgeR_DEG$Id <- rownames(edgeR_DEG)
  edgeR_Up <- edgeR_DEG[edgeR_DEG$Change=='Up',c('Id','Change')]
  edgeR_Down <- edgeR_DEG[edgeR_DEG$Change=='Down',c('Id','Change')]
  Up_inter <- inner_join(limma_Up,edgeR_Up,by='Id')
  Down_inter <- inner_join(limma_Down,edgeR_Down,by='Id')
  print(paste0('Up: ',dim(Up_inter)[1]))
  print(paste0('Down: ',dim(Down_inter)[1]))
  df <- rbind(Down_inter,Up_inter)
  df <- df[,c('Id','Change.x')]
  colnames(df) <- c('Id','Change')
  write.table(df,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/',sheet,'/Intersection_2p_',sheet,'_',Log2FC_cutoff,'.csv'),sep=',')
  stat <- data.frame(
    c(
      as.numeric(
        ifelse(
          is.na(table(df$Change)['Down']),0,table(df$Change)['Down'])
      )
      ,as.numeric(
        ifelse(
          is.na(table(df$Change)['Up']),0,table(df$Change)['Up'])
      )
    )
  )
  colnames(stat) <- paste0('Intersection_',sheet,'_',Log2FC_cutoff)
  return(list(df,stat))
}
intersect_stats <- function(cutoffs,sheet){
  stats <- data.frame('NULL'=c(0,0),row.names = c('Down','Up'))
  for(cutoff in cutoffs){
    res <- interact_DEGs(sheet = sheet,Log2FC_cutoff = cutoff)
    stat <- res[2]
    stats <- cbind(stats,stat)
  }
  stats <- stats[,-1]
  save_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/',sheet,'/Intersection_2p_',sheet,'_','stats','.csv')
  write.table(stats,save_path,sep=',',row.names = T)
  return(stats)
}
cutoffs <- c(0,0.5,1,1.5,2,2.5,3)
sheet <- 'mRNA'
stats <- intersect_stats(cutoffs = cutoffs,sheet = sheet )

sheet <- 'miRNA'
stats <- intersect_stats(cutoffs = cutoffs,sheet = sheet )

sheet <- 'lncRNA'
stats <- intersect_stats(cutoffs = cutoffs,sheet = sheet )

sheet <- 'circRNA'
stats <- intersect_stats(cutoffs = cutoffs,sheet = sheet )
