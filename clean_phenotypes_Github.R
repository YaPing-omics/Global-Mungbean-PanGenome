#phenotype
f1=list.files(pattern='assoc.txt')

pheno=read.delim('~/Documents/mungbean/Trait/GWAS_phenotype_num.txt', sep='\t')
temp=read.delim('~/Documents/mungbean/Trait/GWAS_phenotype_two.txt', sep='\t')
pheno=merge(pheno, temp, by='MMC')
temp=read.delim('~/Documents/mungbean/Trait/GWAS_phenotype_asia_new_worep.txt', sep='\t')
pheno=merge(pheno, temp, by.x='MMC',by.y='Genotype')
temp=read.delim('~/Documents/mungbean/Trait/GWAS_phenotype_aus_new.txt', sep='\t')
pheno=merge(pheno, temp, by.x='MMC',by.y='Genotype')

pheno=pheno[,c(1,grep('SW',colnames(pheno)))]

colnames(pheno)=gsub('_Minicore.evaluation','',colnames(pheno))
colnames(pheno)=gsub('_NA','',colnames(pheno))
colnames(pheno)=gsub('X100SW_','',colnames(pheno))
colnames(pheno)=gsub('X1000SW_','',colnames(pheno))
colnames(pheno)[3]='Bangladesh_Barisal_2017'
colnames(pheno)[4]='Bangladesh_Gazipur_2017'
colnames(pheno)[5]='Bangladesh_Ishwardi_2017'
colnames(pheno)[15]='Pakistan_2016'
pheno=pheno[,order(colnames(pheno))]
pheno=pheno[,c(20,1:19,21:30)]
colnames(pheno)

write.table(pheno,'~/Documents/mungbean/Trait/GWAS_phenotype_SW_all_renames.txt', sep='\t', quote=F, row.names = F)

library(corrplot)
cor1=cor(pheno[,-1], use='pairwise.complete.obs')
pdf('SW_cor.pdf', width=10, height=10)
corrplot(cor1, type='upper', addCoef.col = 'black', number.cex=0.5, tl.col='black', tl.srt=45, tl.cex=0.8)
dev.off()

library(factoextra)
library(FactoMineR)
pheno=read.table('~/Documents/mungbean/Trait/GWAS_phenotype_SW_all_renames.txt', sep='\t', header=T, stringsAsFactors = F)

pca1=PCA(pheno[,-1], graph=F)
pdf('SW_PCA.pdf', width=10, height=10)
fviz_pca_biplot(pca1, axes=c(1,2), geom.ind='point', pointshape=20, label='var', repel=T, col.var = 'dodgerblue')+theme_minimal()
dev.off()

pca1=PCA(t(pheno[,-1]), graph=F)
pdf('SW_PCA_var.pdf', width=10, height=10)
fviz_pca_ind(pca1, axes=c(1,2), geom.var=c('point','text'), pointshape=20, label='ind', repel=T, col.ind = 'dodgerblue')+labs(title ="PCA")+theme_minimal()
dev.off()

SW=which(colnames(pheno) %in% gsub('.assoc.txt','', gsub('MMC_','',f1)))
pheno1=pheno[,c(1,SW)]
pheno1=pheno1[,c(1,3:14,2,15:19)]
write.table(pheno1,'~/Documents/mungbean/Trait/GWAS_phenotype_SW_PVE02.txt', sep='\t', quote=F, row.names = F)

cor1=cor(pheno[,-1], use='pairwise.complete.obs')
pdf('SW_cor.pdf', width=8, height=8)
corrplot(cor1, type='upper', addCoef.col = 'black', number.cex=0.5, tl.col='black', tl.srt=45, tl.cex=0.8)
dev.off()


