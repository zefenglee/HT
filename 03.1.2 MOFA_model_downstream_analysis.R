library(MOFA2)
library(ggplot2)
library(tidyverse)
rm(list = ls())
# 下游分析
# 1.加载模型
packages <- c('edgeR','Intersect_2p','Intersect_3p')
package <- packages[2]
Log2FC_cutoff <- 1
filepath <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                   ,package,'/',package,'_',Log2FC_cutoff,'/model_HT_',package,'_',Log2FC_cutoff,'.hdf5')
print(filepath)
model <- load_model(filepath)
#添加metadata
HT_metadata <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/样本临床信息_预处理.csv')
head(HT_metadata)
samples_metadata(model) <- HT_metadata
head(model@samples_metadata)
# 1.方差解释度
variance_explained <- data.frame(model@cache$variance_explained$r2_total[[1]])
colnames(variance_explained) = 'variance_explained'
plot_variance_explained(model, x='group', y='factor', plot_total = T)[[2]]
variance_explained
write.table(variance_explained,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                                      ,package,'/',package,'_',Log2FC_cutoff,'/variance_explained.csv'),sep=',')
# 2.factor之间的相关性分析，好的因子分解，每个因子之间的相关性是不一样的
p_corr <- plot_factor_cor(model)
write.table(data.frame(p_corr$corr),paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                                           ,package,'/',package,'_',Log2FC_cutoff,'/factors_corr.csv'),sep=',')

# 3.各因子各组学解释度
variance_explained_per_factor <- as.data.frame(model@cache$variance_explained$r2_per_factor[[1]])
write.table(variance_explained_per_factor,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                                           ,package,'/',package,'_',Log2FC_cutoff,'/variance_explained_per_factor.csv'),sep=',')
plot_variance_explained(model, x="view", y="factor",factors = 'all')
# 4.提取每个样本各因子值，比较HT和HC样本组的因子显著性，用于筛选组间存在显著性差异的因子
factors <- get_factors(model, as.data.frame = T)
factors <- as.data.frame(t(pivot_wider(factors,names_from = sample)))
colnames(factors) <- factors['factor',]
factors <- factors[c(-1,-2),]
factors$Group <- substr(rownames(factors),1,1)
factors['P',] <- 0
for(i in seq((dim(factors)[2]-1))){
  c <- as.numeric(factors[,i][factors$Group=='C'])
  t <- as.numeric(factors[,i][factors$Group=='T'])
  p <- wilcox.test(c,t)
  factors['P',i] <- p$p.value
  print(p$p.value)
}
write.table(factors,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                           ,package,'/',package,'_',Log2FC_cutoff,'/factor_values_per_sample.csv'),sep=',')
# Factor2组间差异显著
# 5.显著性因子的特征展示
# 5.1样本组间显著性因子值的分布图（小提琴图或箱线图）
plot_factor(model, 
            factors = 3, 
            color_by = "Group",
            add_violin = TRUE,
            dodge = TRUE
)
plot_factors(model, 
             factors = c(1,3),
             color_by = "Group"
)
# 6.筛选显著性因子的高权重特征（权重阈值待确定)
weights <- as.data.frame(pivot_wider(get_weights(model, as.data.frame = T),names_from = 'factor'))
head(weights, n=3)
write.table(weights,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                           ,package,'/',package,'_',Log2FC_cutoff,'/feature_weight_per_factor.csv'),sep=',')
# biomark_list <- c()
# for(i in seq(dim(weights)[2]-2)){
#   biomark_list = c(biomark_list,as.character(weights$view[abs(weights[,3:dim(weights)[2]][,i])>0.2]))
# }
# table(biomark_list)
factor1_up_cutoff <- mean(weights$Factor1)+2*sd(weights$Factor1)
factor1_down_cutoff <- mean(weights$Factor1)-2*sd(weights$Factor1)
biomark_list1 <- weights[(weights$Factor1>factor1_up_cutoff)|(weights$Factor1<factor1_down_cutoff),c('feature','view','Factor1')]
biomark_list1$Factor <- 'Facror1'
table(biomark_list1$view)

factor2_up_cutoff <- mean(weights$Factor2)+2*sd(weights$Factor2)
factor2_down_cutoff <- mean(weights$Factor2)-2*sd(weights$Factor2)
biomark_list2 <- weights[(weights$Factor2>factor2_up_cutoff)|(weights$Factor2<factor2_down_cutoff),c('feature','view','Factor2')]
biomark_list2$Factor <- 'Facror2'
table(biomark_list2$view)

factor3_up_cutoff <- mean(weights$Factor3)+2*sd(weights$Factor3)
factor3_down_cutoff <- mean(weights$Factor3)-2*sd(weights$Factor3)
biomark_list3 <- weights[weights$Factor3>factor3_up_cutoff | weights$Factor3<factor3_down_cutoff,c('feature','view','Factor3')]
biomark_list3$Factor <- 'Facror3'
table(biomark_list3$view)

factor5_up_cutoff <- mean(weights$Factor5)+2*sd(weights$Factor5)
factor5_down_cutoff <- mean(weights$Factor5)-2*sd(weights$Factor5)
biomark_list5 <- weights[(weights$Factor5>factor5_up_cutoff)|(weights$Factor5<factor5_down_cutoff),c('feature','view','Factor5')]
biomark_list5$Factor <- 'Facror5'
table(biomark_list5$view)

factor10_up_cutoff <- mean(weights$Factor10)+2*sd(weights$Factor10)
factor10_down_cutoff <- mean(weights$Factor10)-2*sd(weights$Factor10)
biomark_list10 <- weights[weights$Factor10>factor10_up_cutoff | weights$Factor10<factor10_down_cutoff,c('feature','view','Factor10')]
biomark_list10$Factor <- 'Facror10'
table(biomark_list10$view)

biomarkers <- rbind(biomark_list1[,c(1,2,4)],biomark_list3[,c(1,2,4)])
biomarkers <- biomarkers[!duplicated(biomarkers$feature),]
write.table(biomarkers
            ,paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/'
                    ,package,'/',package,'_',Log2FC_cutoff,'/features_selected_from_factor.csv'),sep=',')
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
bitr(biomarkers$feature[str_detect(biomarkers$feature,'ENSG')],fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')

