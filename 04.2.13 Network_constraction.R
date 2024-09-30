library(clusterProfiler)
library(org.Hs.eg.db)
rm(list = ls())
Interactome <- read.delim('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/12.ALL/Interactome.txt')
biomarkers_df <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/edgeR/edgeR_1/features_selected_from_factor.csv')
table(biomarkers_df$view)
biomarkers <- biomarkers_df$feature
length(biomarkers)
sum(biomarkers %in% c(Interactome$Gene1,Interactome$Gene2))
ENSGs <- biomarkers[str_detect(biomarkers,'ENSG')]
bitr(ENSGs,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
network <- Interactome[(Interactome$Gene1 %in% biomarkers)&(Interactome$Gene2 %in% biomarkers), ]
sum(!duplicated(c(network$Gene1,network$Gene2)))
biomarkers_included <- c(network$Gene1,network$Gene2)[!duplicated(c(network$Gene1,network$Gene2))]
biomarkers_included_df <- biomarkers_df[biomarkers %in% biomarkers_included,]
table(biomarkers_included_df$view)
write.table(network,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/edgeR/network.txt',row.names = F,quote = F,sep='\t')
write.table(biomarkers_included_df,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/edgeR/biomarker_included_in_Interactome.txt',row.names = F,quote = F,sep='\t')

rm(list = ls())
Interactome <- read.delim('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Intersctomes/12.ALL/Interactome.txt')
biomarkers_df <- read.csv('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/03.MOFA/Intersect_2p/Intersect_2p_1/features_selected_from_factor.csv')
table(biomarkers_df$view)
biomarkers <- biomarkers_df$feature
length(biomarkers)
sum(biomarkers %in% c(Interactome$Gene1,Interactome$Gene2))
ENSGs <- biomarkers[str_detect(biomarkers,'ENSG')]
bitr(ENSGs,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
network <- Interactome[(Interactome$Gene1 %in% biomarkers)&(Interactome$Gene2 %in% biomarkers), ]
sum(!duplicated(c(network$Gene1,network$Gene2)))
biomarkers_included <- c(network$Gene1,network$Gene2)[!duplicated(c(network$Gene1,network$Gene2))]
biomarkers_included_df <- biomarkers_df[biomarkers %in% biomarkers_included,]
table(biomarkers_included_df$view)
write.table(network,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/Intersect_2p/network.txt',row.names = F,quote = F,sep='\t')
write.table(biomarkers_included_df,'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/04.Network/Intersect_2p/biomarker_included_in_Interactome.txt',row.names = F,quote = F,sep='\t')
