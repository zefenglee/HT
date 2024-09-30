rm(list = ls())
library(openxlsx)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(GGally)
library(stringr)
biomarkers_list <- read.delim('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/Intersect_2p_1/features_selected_from_factor.csv',sep=',')
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

#expr1 <- data.frame(t(fpkm_list$circRNA))
expr1 <- t(fpkm_list$miRNA)
#biomarkers1 <- biomarkers_list$circRNA 
biomarkers1 <-colnames(expr1)
#rownames(expr1) <- biomarkers1
#expr2 <- as.data.frame(fpkm_list$miRNA)
expr2 <- t(fpkm_list$mRNA)
biomarkers2 <- colnames(expr2)
#expr <- as.data.frame(t(rbind(expr1,expr2)))
expr <- as.data.frame(cbind(expr1,expr2))
expr$Group <- Group
expr$Group <- as.factor(expr$Group)
for(x_ in biomarkers1){
  for(y_ in biomarkers2){
    t <- cor.test(expr[,x_],expr[,y_])
    if(t$p.value<1){
      p <- ggscatter(expr
                     ,x=x_
                     ,y=y_
                     ,add = 'reg.line'
                     ,conf.int = T
                     ,cor.coef = T
                     ,cor.method = 'pearson'
                     # ,xlab = paste0(x_,'(circRNA)')
                     ,xlab=paste0(str_replace_all(pattern = '_',replacement = '-',x_),'(miRNA)')
                     ,ylab=paste0(str_replace_all(pattern = '_',replacement = '-',y_),'(mRNA)')
                     ,title = paste("Relationship between\n"
                                    ,str_replace_all(pattern = '_',replacement = '-',x_)
                                    ,'and',str_replace_all(pattern = '_',replacement = '-',y_)
                                    ,sep = ' ')
                     # ,title = paste("Relationship between\n"
                     #                ,x_
                     #                ,'and',str_replace_all(pattern = '_',replacement = '-',y_)
                     #                ,sep = ' ')
                     )+
        geom_point(aes(colour=Group))+
        scale_discrete_manual(values=c("#9394E7","#F0988c")
                              ,aesthetics = 'colour'
                              ,labels = c("HC", "HT")
        )+
        theme_pubr()+
        theme(#axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          #axis.title.y = element_blank(),
          #panel.grid=element_blank(),
          #panel.border = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          plot.title = element_text(size=10,hjust=0.5),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8)
        )
      savepath <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/06.Co_expression/mi_m/',x_,'_',y_,'.jpg')
      jpeg(file = savepath,width =1000, height = 1000, res=300)
      plot(p)
      dev.off()
    }
  }
}




expr1 <- data.frame(t(fpkm_list$circRNA))
biomarkers1 <- biomarkers_list$circRNA 
rownames(expr1) <- biomarkers1
expr2 <- as.data.frame(fpkm_list$miRNA)
biomarkers2 <- rownames(expr2)
expr <- as.data.frame(t(rbind(expr1,expr2)))
expr$Group <- Group
expr$Group <- as.factor(expr$Group)
for(x_ in biomarkers1){
  for(y_ in biomarkers2){
    t <- cor.test(expr[,x_],expr[,y_])
    if(t$p.value<0.05){
      p <- ggscatter(expr
                     ,x=x_
                     ,y=y_
                     ,add = 'reg.line'
                     ,conf.int = T
                     ,cor.coef = T
                     ,cor.method = 'pearson'
                     ,xlab = paste0(x_,'(circRNA)')
                     #,xlab=paste0(str_replace_all(pattern = '_',replacement = '-',x_),'(miRNA)')
                     ,ylab=paste0(str_replace_all(pattern = '_',replacement = '-',y_),'(miRNA)')
                     # ,title = paste("Relationship between\n"
                     #                ,str_replace_all(pattern = '_',replacement = '-',x_)
                     #                ,'and',str_replace_all(pattern = '_',replacement = '-',y_)
                     #                ,sep = ' ')
                     ,title = paste("Relationship between\n"
                                    ,x_
                                    ,'and',str_replace_all(pattern = '_',replacement = '-',y_)
                                    ,sep = ' ')
      )+
        geom_point(aes(colour=Group))+
        scale_discrete_manual(values=c("#9394E7","#F0988c")
                              ,aesthetics = 'colour'
                              ,labels = c("HC", "HT")
        )+
        theme_pubr()+
        theme(#axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          #axis.title.y = element_blank(),
          #panel.grid=element_blank(),
          #panel.border = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          plot.title = element_text(size=10,hjust=0.5),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8)
        )
      savepath <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/06.Co_expression/circ_mi//',x_,'_',y_,'.jpg')
      jpeg(file = savepath,width =1000, height = 1000, res=300)
      plot(p)
      dev.off()
    }
  }
}
# p <- ggscatter(expr
#           ,x=x_
#           ,y=y_
#           ,add = 'reg.line'
#           ,conf.int = T
#           ,cor.coef = T
#           ,cor.method = 'pearson'
#           ,xlab=str_replace_all(pattern = '_',replacement = '-',x_)
#           ,ylab=str_replace_all(pattern = '_',replacement = '-',y_)
#           ,title = paste("Relationship between\n"
#                          ,str_replace_all(pattern = '_',replacement = '-',x_)
#                          ,'and',str_replace_all(pattern = '_',replacement = '-',y_)
#                          ,sep = ' '))+
#   geom_point(aes(colour=Group))+
#   scale_discrete_manual(values=c("#9394E7","#F0988c")
#                         ,aesthetics = 'colour'
#                         ,labels = c("HC", "HT")
#                         )+
#   theme_pubr()+
#     theme(#axis.text.y = element_blank(),
#           #axis.ticks.y = element_blank(),
#           #axis.title.y = element_blank(),
#           #panel.grid=element_blank(),
#           #panel.border = element_blank(),
#           axis.title.x = element_text(size=10),
#           axis.title.y = element_text(size=10),
#           plot.title = element_text(size=12,hjust=0.5),
#           legend.title = element_text(size=10),
#           legend.text = element_text(size=10)
#           )
# p  
# savepath <- paste0('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/06.Co_expression/lnc_mi/',x_,'_',y_,'.jpg')
# jpeg(file = savepath,width =1000, height = 1000, res=300)
# plot(p)
# dev.off()
# cor_mRNA_miRNA <- cor(t(rbind(fpkm_list$mRNA,fpkm_list$miRNA)))[biomarkers_list$mRNA,biomarkers_list$miRNA]
# cor_lncRNA_miRNA <- cor(t(rbind(fpkm_list$lncRNA,fpkm_list$miRNA)))[biomarkers_list$lncRNA,biomarkers_list$miRNA]
# cor_lncRNA_mRNA <- cor(t(rbind(fpkm_list$lncRNA,fpkm_list$mRNA)))[biomarkers_list$lncRNA,biomarkers_list$mRNA]
# cor_circRNA_miRNA <- cor(t(rbind(fpkm_list$circRNA,fpkm_list$miRNA)))[biomarkers_list$circRNA,biomarkers_list$miRNA]

# ggpairs(expr[,1:(dim(expr)[2]-1)]
#         ,showStrips = T,ggplot2::aes(color=Group)) + 
#   theme(axis.text = element_text(colour = "black", size = 11),
#         strip.background = element_rect(fill = "#d63d2d"),
#         strip.text = element_text(colour = "white", size = 12,
#                                   face = "bold"))
# 
# 
# 
# 
# 
# 
# # x_ = biomarkers1[1]
# # print(x_)
# # y_ = biomarkers2[1]
# # print(y_)
# 
# p <- ggscatterstats(expr, 
#                 x ="ENSG00000273448",
#                 y ="hsa-miR-374a-5p",
#                type = "pearson",
#                centrality.para = "median",                              
#                margins = "both",                                         
#                xfill = 'red',
#                yfill = 'blue', 
#                marginal.type = "densigram",   #类型可以换成density,boxplot,violin,densigram
#                title = paste("Relationship between",x_,'and',y_,sep = ' '),
#                messages = FALSE,
#                
#                )+
#   geom_point(aes(colour=Group))+
#   # scale_discrete_manual(values=c("#9394E7","#F0988c")
#   #                       ,aesthetics = 'colour'
#   #                       ,labels = c("HC", "HT")
#   # 
#   #                       )+
#   theme_pubr()+
#     theme(#axis.text.y = element_blank(),
#           #axis.ticks.y = element_blank(),
#           #axis.title.y = element_blank(),
#           #panel.grid=element_blank(),
#           #panel.border = element_blank(),
#           axis.title.x = element_text(size=15),
#           axis.title.y = element_text(size=15),
#           plot.title = element_text(size=17,hjust=0.5),
#           legend.title = element_text(size=15),
#           legend.text = element_text(size=12)
#           )
# p
# # color_set1 <- c("#934B43",'#F1D77E','#9394E7','#D76364','#B1CE46','#5F97D2','#EF7A6D','#63E398','#9DC3E7')
# # color_set2 = c('#A1A9D0','#F0988c','#B883D4','#96CCCB','#CFEAF1','#C4A5DE','#F6CAE5','#9DC3E7','#9E9E9E')
# # color_set3 = c('#8ECFC9','#FFBE7A','#FA7F6F','#82B0D2','#BEB8DC','#E7DAD2','#999999')




