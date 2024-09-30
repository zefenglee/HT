rm(list = ls())
library(openxlsx)
library(MOFA2)
load_data <- function(){
  mRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'mrna_fpkm',rowNames =T))
  columns <- colnames(mRNA_expr)  
  miRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'mirna_fpkm',rowNames = T)[,columns])
  lncRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'lncrna_fpkm',rowNames = T)[,columns])
  circRNA_expr <- as.matrix(read.xlsx('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/RNA_seq.xlsx',sheet = 'circrna_fpkm_filtered',rowNames = T)[,columns])
  datas_list <- list('mRNA'=mRNA_expr,'miRNA'=miRNA_expr,'lncRNA'=lncRNA_expr,'circRNA'=circRNA_expr)
  return(datas_list)
}
# load_cols <- function(Log2FC_cutoff){
#   mRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect/mRNA/Intersection_mRNA_',Log2FC_cutoff,'_cols.csv')
#   mRNA_cols <- read.csv(mRNA_path)$Id
#   miRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect/miRNA/Intersection_miRNA_',Log2FC_cutoff,'_cols.csv')
#   miRNA_cols <- read.csv(miRNA_path)$Id
#   lncRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect/lncRNA/Intersection_lncRNA_',Log2FC_cutoff,'_cols.csv')
#   lncRNA_cols <- read.csv(lncRNA_path)$Id
#   circRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect/circRNA/Intersection_circRNA_',Log2FC_cutoff,'_cols.csv')
#   circRNA_cols <- read.csv(circRNA_path)$Id
#   return(list('mRNA'=mRNA_cols,'miRNA'=miRNA_cols,'lncRNA'=lncRNA_cols,'circRNA'=circRNA_cols))
# }
load_cols <- function(Log2FC_cutoff){
  mRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/mRNA/Intersection_2p_mRNA_',Log2FC_cutoff,'.csv')
  mRNA_cols <- read.csv(mRNA_path)$Id
  miRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/miRNA/Intersection_2p_miRNA_',Log2FC_cutoff,'.csv')
  miRNA_cols <- read.csv(miRNA_path)$Id
  lncRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/lncRNA/Intersection_2p_lncRNA_',Log2FC_cutoff,'.csv')
  lncRNA_cols <- read.csv(lncRNA_path)$Id
  circRNA_path <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/02.DEA/Intersect_2p/circRNA/Intersection_2p_circRNA_',Log2FC_cutoff,'.csv')
  circRNA_cols <- read.csv(circRNA_path)$Id
  return(list('mRNA'=mRNA_cols,'miRNA'=miRNA_cols,'lncRNA'=lncRNA_cols,'circRNA'=circRNA_cols))
}
load_subdata <- function(cols_list,data_list){
  data_list$mRNA <- scale(data_list$mRNA[cols_list$mRNA,])
  data_list$miRNA <- scale(data_list$miRNA[cols_list$miRNA,])
  data_list$lncRNA <- scale(data_list$lncRNA[cols_list$lncRNA,])
  data_list$circRNA <- scale(data_list$circRNA[cols_list$circRNA,])
  return(data_list)
}

datas_list <- load_data()
lapply(datas_list, dim)
Log2FC_cutoff <- 1
cols_list <- load_cols(Log2FC_cutoff = Log2FC_cutoff )
data_list <- load_subdata(cols_list,datas_list)
lapply(data_list, dim)
HT_data <- data_list
lapply(HT_data,dim)
# 2.metadata导入
#HT_metadata <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/样本临床信息_预处理.csv') 
#head(HT_metadata)
# 3.Create the MOFA object
MOFAobject <- create_mofa(HT_data)
print(MOFAobject)
#load('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Codes/202312_Hashimotos_thyroiditis/MOFA_exploration/CLL_data.RData')
# 2.可视化数据结构
plot_data_overview(MOFAobject)
# 3.定义选项
# 3.1 定义数据选项
# scale_groups：如果组有不同的范围/方差，那么将每个组缩放到单位方差是一个好的做法，默认值为  FALSE
# scale_views：如果视图有不同的范围/方差，那么将每个视图缩放到单位方差是一个好的做法，默认值为  FALSE
data_opts <- get_default_data_options(MOFAobject)
#data_opts$views <- c('mRNA','miRNA','lncRNA','circRNA')
data_opts
# 3.2 定义模型选项，这部分一般不改动，用默认参数
model_opts <- get_default_model_options(MOFAobject)
model_opts$spikeslab_weights <- T
model_opts$num_factors <- 15
model_opts
# 3.3 定义训练选项
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$freqELBO <- 1
train_opts$maxiter <- 10000
train_opts

# 4.训练和构建MOFA模型
# 4.1 准备mofa对象
MOFAobject <- prepare_mofa(object = MOFAobject
                           ,data_options = data_opts
                           ,model_options = model_opts
                           ,training_options = train_opts)
# 4.2 训练mofa模型
outfile <- file.path(paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/model_HT_Intersect_2p_',Log2FC_cutoff,'.hdf5'))
print(outfile)
MOFAobject.trained <- run_mofa(MOFAobject,outfile = outfile,save_data = T,
                               use_basilisk = F)
saveRDS(MOFAobject,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/MOFA2_HT_Intersect_2p_',Log2FC_cutoff,'.rds'))

