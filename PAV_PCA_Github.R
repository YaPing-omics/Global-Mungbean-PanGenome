
#PAV PCA
d1=read.delim('./5KPAV_dist.txt', sep='\t', stringsAsFactors = F, skip=5, header=F)
d1[,1]=gsub('[.]','-',d1[,1])
d1[1,1]='153_ML1628'

d1=merge(d1, group[,1:2], by='V1')
d1$col=d1$group
d1[which(d1$col=='admix_cultivar'),'col']='grey'
d1[which(d1$col=='admix_wild'),'col']='grey'
d1[which(d1$col=='Rad1'),'col']='skyblue'
d1[which(d1$col=='Rad2'),'col']='blue'
d1[which(d1$col=='SubAF'),'col']='yellow'
d1[which(d1$col=='SubAS_SA'),'col']='peru'
d1[which(d1$col=='SubAS_SEA'),'col']='orange'
d1[which(d1$col=='SubAUe'),'col']='red'
d1[which(d1$col=='SubAUw'),'col']='#900603'
d1[which(d1$col=='SubIDN'),'col']='pink'
d1=d1[order(d1$group),]

pca=prcomp(d1[,-c(1,782,783)],center = T, scale. = T)
x1=pca$sdev[1]^2/sum(pca$sdev^2)
x1=pca$sdev[2]^2/sum(pca$sdev^2)

pdf('PAV_PCA.pdf', width=4, height=4)
par(ps=8)
plot(pca$x[,1], pca$x[,2], col=d1$col, pch=16, xlab='PC1 (96.2%)', ylab='PC2 (2.3%)')
dev.off()
write.table(data.frame(cbind(d1[,1],pca$x)), 'PAV_PCA.txt', sep='\t', quote = F, row.names = F)
