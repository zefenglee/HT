rm(list = ls())
gc()
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)
setwd('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/')
RAIN <- read.table('06.RAIN/9606.v1.experiments.tsv',sep = '\t')
summary(RAIN$V6)
RAIN <- RAIN[RAIN$V6>=quantile(RAIN$V6,0.75),c('V2','V3')]
colnames(RAIN) <- c('Gene1','Gene2')
sum(str_detect(RAIN$Gene1,'hsa-mir'))
sum(str_detect(RAIN$Gene2,'hsa-mir'))
id_trans <- bitr(c(RAIN$Gene1[str_detect(RAIN$Gene1,'ENSP')],RAIN$Gene2[str_detect(RAIN$Gene2,'ENSP')])
     ,fromType = 'ENSEMBLPROT',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
id_trans[id_trans$ENSEMBLPROT=='ENSP00000300086','SYMBOL']="TERF2IP"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000355511','SYMBOL']="CHML"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000420182','SYMBOL']="IQCJ-SCHIP1"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000246166','SYMBOL']="FNTB"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000424598','SYMBOL']="DEFB4B"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000291554','SYMBOL']="CRYAA"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000296600','SYMBOL']="ZNF474"
id_trans[id_trans$ENSEMBLPROT=='ENSP00000361894','SYMBOL']="GNMT"
id_trans <- id_trans[!duplicated(id_trans$ENSEMBLPROT),]
sum(duplicated(id_trans$ENSEMBLPROT))
for(i in seq(length(RAIN$Gene1))){
  gene1 <- RAIN$Gene1[i]
  if(str_detect(gene1,'ENSP')){
    bool_ <- id_trans$ENSEMBLPROT %in% gene1
    if(sum(bool_)!=0){
      RAIN$Gene1[i] <- id_trans$SYMBOL[bool_]
    }else{
      RAIN$Gene1[i] <- ''
    }
  }
  
}
for(i in seq(length(RAIN$Gene2))){
  gene2 <- RAIN$Gene2[i]
  if(str_detect(gene2,'ENSP')){
    bool_ <- id_trans$ENSEMBLPROT %in% gene2
    if(sum(bool_)!=0){
      RAIN$Gene2[i] <- id_trans$SYMBOL[bool_]
    }else{
      RAIN$Gene2[i] <- ''
    }
  }
  
}
RAIN <- RAIN[RAIN$Gene1!=''&RAIN$Gene2!='',]
sum(duplicated(paste0(RAIN$Gene1,RAIN$Gene2)))
sum(!duplicated(c(RAIN$Gene1,RAIN$Gene2)))
write.table(RAIN,'12.ALL/06.RAIN.txt',row.names = F,quote = F)
