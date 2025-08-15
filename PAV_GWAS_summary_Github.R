library(pheatmap)
library(qqman)

f1=list.files(pattern='log.txt')

pve=c()
for(i in seq(length(f1))){
  temp=read.delim(f1[i])
  pve[i]=temp[23,1]
}
pve=gsub('## pve estimate in the null model = ','',pve)
pve=as.numeric(pve)
key=which(pve>0.2)

f1=f1[key]
f1=gsub('.log.txt','.assoc.txt',f1)

f2=gsub('MMC_PAV_X','',f1)
f2=gsub('MMC_PAV_','',f2)
f2=gsub('100SW_','',f2)
f2=gsub('1000SW_','',f2)
f2=gsub('.assoc.txt','',f2)
f2=gsub('_Minicore.evaluation','',f2)
f2=gsub('NA_','',f2)
f2[3]='Bangladesh_Barisal_2017'
f2[4]='Bangladesh_Ishwardi_2017'
f2[1]='mean'

#combine GWAS result
d1=read.delim(f1[1], sep='\t', stringsAsFactors = F)

result=data.frame(rs=d1[,2])
result1=data.frame()
for(i in 1:length(f1)){
  d1=read.table(f1[i], sep='\t', header=T, stringsAsFactors = F)
  colnames(d1)[15]=f2[i]
  result=merge(result,d1[,c(2,15)], by='rs', all.x=T)
  
  d1=d1[which(d1[,15] < 0.001),]
  colnames(d1)[15]='p_score'
  d1$trait=rep(paste0('SW_', f2[i]))
  result1=rbind(result1,d1) 
}

length(unique(result1$rs))
write.table(result1, 'SW_PAV_sig_summary.txt', quote=F, row.names = F, sep='\t')

#modify LOD 
result.g=result
padj=-log10(0.05/nrow(result.g))
for(i in 2:ncol(result.g)){
  result.g[,i]=-log10(result.g[,i])
  result.g[which(result.g[,i] < 3),i]=1
  result.g[which(result.g[,i] >= 3 & result.g[,i] < padj),i]=2
  result.g[which(result.g[,i] >= padj),i]=3
}

qtl=c()
for(i in 2:ncol(result.g)){
  qtl[[i]]=which(result.g[,i]>1)
}
qtl=unique(unlist(qtl))
qtl=result.g[qtl,'rs']

con=c()
for(i in 1:nrow(result.g)) con[i]=length(which(na.omit(result.g[i,-1])>1))
#con.1=result.g[which(con>10),'rs']
con.1=which(con>10)

result.g[con.1,-1]=result.g[con.1,-1]+3

fig=result.g[which(result.g$rs %in% qtl),order(colnames(result.g))]
head(fig)
fig=fig[order(fig$rs),]
fig=fig[,c(16,18,1:15,17,19:20)]
colnames(fig)
row.names(fig)=fig$rs
fig$chr=unlist(lapply(strsplit(fig$rs,split='[.]'),'[[',1))
fig$pos=unlist(lapply(strsplit(fig$rs,split='[.]'),'[[',2))

fig1=fig[grep('Chr',fig$rs),]
fig1$chr=as.numeric(gsub('Chr','',fig1$chr))
fig1$pos=as.numeric(fig1$pos)
fig1=fig1[order(fig1$chr, fig1$pos),]
fig2=fig[-grep('Chr',fig$rs),]
fig3=data.frame(rbind(fig1,fig2))
fig3=fig3[rev(seq(nrow(fig3))),]
fig3$No=rev(seq(nrow(fig3)))
colnames(fig3)
fig3=fig3[,c(1,21:23,2:20)]
fig3=fig3[,c(1:12,14,13,16,15,17:23)]

colnames(fig3)[5]='Mean'
colnames(fig3)[6]='Bangladesh_Barisal_2017'
colnames(fig3)[7]='Bangladesh_Ishwardi_2017'
colnames(fig3)[8]='Bangladesh_Ishwardi_2018'
colnames(fig3)[9]='Bangladesh_Rangpur_2020'
colnames(fig3)[13]='India_Kanpur_2018'
colnames(fig3)[14]='India_Dharwad_2021'
colnames(fig3)[15]='Kenya_Katumani_2018'
colnames(fig3)[16]='Kenya_KampiMawe_2019'

anno=data.frame(Chr=fig3$chr)
rownames(anno)=fig3$rs
anno[-grep('Chr', row.names(anno)),'Chr']='0'

fig3.1=fig3[which(fig3$chr %in% seq(11)),]

# set colours
anno_colori=list(Chr=c('1'='grey20','2'='grey80','3'='grey20','4'='grey80','5'='grey20','6'='grey80','7'='grey20','8'='grey80','9'='grey20','10'='grey80','11'='grey20','0'='grey80'))

pdf('PAV_SW.pdf', width=12, height=4)
pheatmap(t(fig3.1[rev(seq(nrow(fig3.1))),-c(1:4)]),color=c('azure','skyblue','blue','mistyrose','pink','red'),legend_breaks=seq(6), legend_labels = c('specific = 0', 'specific > 3', 'specific > 5', 'consistent = 0', 'consistent > 3', 'consistent > 5'), cluster_rows = F, cluster_cols = F, fontsize = 10, fontsize_row = 10, fontsize_col = 2, main='', angle_col = '90', na_col='white', cellwidth = 2.5, cellheight= 8, show_colnames = T,annotation_col = anno, annotation_colors = anno_colori)
dev.off()

man1=result[which(result$rs %in% fig3.1$rs), c(1,2)]
man1=merge(fig3[,c(1:3)], man1, by='rs', )
man1$chr=as.numeric(man1$chr)
man1$pos1=as.numeric(unlist(lapply(strsplit(man1$rs,split='[.]'),'[[',2)))
man1=man1[order(man1$chr,man1$pos1),]
man1$pos=seq(nrow(man1))
head(man1)

man2=man1[which(man1$rs %in% result.g[which(con>10),'rs']),]
man3=man1[which(man1$rs %in% c('Chr1.2667','Chr1.3000','Chr2.2867','Chr7.790','Chr9.1538','Chr11.1261')),]

pdf('PAV_SW_mean.pdf', width=8, height=4)
manhattan(man1, chr = 'chr', bp = 'pos', p = 'SW_mean', snp = 'rs', col=c('dodgerblue','gold'), suggestiveline = F, genomewideline = F, ylim=c(0,5.2))
points(x=man2$pos, y=-log10(man2$SW_mean), pch=16, col='red')
points(x=man3$pos, y=-log10(man3$SW_mean), pch=4, col='red')
dev.off()

#PAV boxplot
library(ggplot2)
pheno=read.delim('~/Documents/mungbean/Trait/GWAS_phenotype_SW_PVE02.txt', sep='\t', stringsAsFactors = F)
pheno$mean=apply(pheno[,-1],1,function(x) mean(x, na.rm=T))
pheno$mean=round(pheno$mean, digit=2)

snp0=read.delim('../8candidates.hmp.txt', sep='\t', stringsAsFactors = F)
snp.n=snp0[,1]
snp0=data.frame(snp1=t(snp0[,-c(1:11)]))
snp0$Genotype=row.names(snp0)
snp0$Genotype=gsub('[.]','-',snp0$Genotype)

i=7
snp=snp0[,c(9,i)]
colnames(snp)[2]='snp1'
snp[which(snp$snp1=='G'),'snp1']='+'
head(snp)
snp=merge(snp,pheno[,c(1,20)], by.x='Genotype', by.y='MMC')

pdf(paste0(snp.n[i],'.pdf'), height = 1.4, width=2)
ggplot(snp, aes(x=snp1, y=mean, fill=snp1, color=snp1)) +
  geom_boxplot() +
  scale_fill_manual('PAV',values=c("orange", "blue")) +
  scale_color_manual('PAV',values=c("orange4", "navyblue")) +
  xlab(paste0(snp.n[i])) +
  ylab('Mean') +
  theme_classic()
dev.off()
