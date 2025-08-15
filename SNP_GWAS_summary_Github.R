library(pheatmap)

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

f2=gsub('MMC_woLDprune_X100SW_','',f1)
f2=gsub('MMC_X100SW_','',f2)
f2=gsub('MMC_X1000SW_','',f2)
f2=gsub('.assoc.txt','',f2)
f2=gsub('NA_','',f2)
f2=gsub('_Minicore.evaluation','',f2)
f2[3]='Bangladesh_Barisal_2017'
f2[4]='Bangladesh_Ishwardi_2017'

#combine GWAS result
d1=read.delim(f1[1], sep='\t', stringsAsFactors = F)
d1$rs=paste(d1$chr, d1$ps, sep='_')
d1[which(d1$p_score< 10^-6),1:3]

library(qqman)

result=data.frame(rs=d1[,2])
result1=data.frame()
for(i in 6:length(f1)){
  d1=read.table(f1[i], sep='\t', header=T, stringsAsFactors = F)
  d1$rs=paste(d1$chr, d1$ps, sep='_')
  
  png(paste0('SW_',f2[i],'.png'), width=480, height= 240)
  manhattan(d1, chr = 'chr', bp = 'ps', p = 'p_score', snp = 'rs', col=c('dodgerblue','gold'), suggestiveline = F, genomewideline = F, main=paste0('SW_',f2[i]))
  dev.off()
  
  colnames(d1)[15]=f2[i]
  result=merge(result,d1[,c(2,15)], by='rs', all.x=T)
  
  d1=d1[which(d1[,15] < 0.000001),]
  colnames(d1)[15]='p_score'
  d1$trait=rep(paste0('SW_',f2[i]))
  result1=rbind(result1,d1) 
}

write.table(result1, 'SW_SNP_sig_summary.txt', quote=F, row.names = F, sep='\t')

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
for(i in 3:ncol(result.g)){
  qtl[[i]]=which(result.g[,i]>1)
}
qtl=unique(unlist(qtl))
qtl=result.g[qtl,'rs']

con=c()
for(i in 1:nrow(result.g)) con[i]=length(which(na.omit(result.g[i,-c(1:2)])>1))
#con.1=result.g[which(con>10),'rs']
con.1=which(con>10)

result.g[con.1,-c(1:2)]=result.g[con.1,-c(1:2)]+3

fig=result.g[which(result.g$rs %in% result.g[con.1,1]),order(colnames(result.g))]
fig=fig[,c(17,13,1:12,14:16,18:20)]
fig=fig[order(fig$rs),]
row.names(fig)=fig$rs
fig$chr=unlist(lapply(strsplit(fig$rs,split='_'),'[[',1))
fig$pos=unlist(lapply(strsplit(fig$rs,split='_'),'[[',2))

library(qqman)
png('SNP_mean_SW.png', width=960, height=480, pointsize = 24)
manhattan(d1, chr = 'chr', bp = 'ps', p = 'p_score', snp = 'rs', col=c('dodgerblue','gold'), suggestiveline = F, genomewideline = F)
dev.off()
