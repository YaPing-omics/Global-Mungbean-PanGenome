pav1=read.delim('./Chr_PAV_sum.txt', sep='\t', stringsAsFactors = F, check.names = F)
pav1=pav1[,c(1:5,order(colnames(pav1)[-c(1:5)])+5)]

core.cry=c()
for(i in 1:nrow(pav1)) core.cry[i]=length(which(pav1[i,-c(1:5)]=='PRESENT'))
core.pan=c()
for(i in 1:nrow(pav2.1)) core.pan[i]=length(which(pav2.1[i,-c(1:5)]=='PRESENT'))

present=c()
for(i in 1:nrow(pav1)){
  present[i]=length(which(pav1[i,-c(1:5)]=='PRESENT'))
}

result=data.frame()
for(i in 6:ncol(pav1)){
  temp=c(colnames(pav1)[i],length(which(pav1[,i]=='LOST')),length(which(pav1[,i]=='PRESENT')))
  result=rbind(result, temp)
}
colnames(result)=c('Accession','LOST','PRESENT')

pav2=read.delim('./pan_PAV_sum.txt', sep='\t', stringsAsFactors = F, check.names = F)
pav2=pav2[,c(1:5,order(colnames(pav2)[-c(1:5)])+5)]
pav2.1[1:5,1:10]
pav2.1=pav2.1[,c(1:5,order(colnames(pav2.1)[-c(1:5)])+5)]
pav2.1=pav2[which(pav2$GeneID %in% clean[,1]),]

tab1=table(group$group)

pav.todo=data.frame(rbind(pav1, pav2.1))
colnames(pav.todo)=gsub('[.]','-',colnames(pav.todo))
colnames(pav.todo)[6]='153_ML1628'
for(i in 6:ncol(pav.todo)) pav.todo[which(pav.todo[,i]=='PRESENT'),i]=1
for(i in 6:ncol(pav.todo)) pav.todo[which(pav.todo[,i]=='LOST'),i]=0
pav.todo[1:5,1:15]
write.table(pav.todo, 'PAV_table.txt', sep='\t', quote=F, row.names = F)

pav.todo2=pav.todo
for(i in 6:ncol(pav.todo2)) pav.todo2[,i]=as.numeric(pav.todo2[,i])
pav.todo2$sum=apply(pav.todo2[,-c(1:5)], 1, sum)
core.gene=data.frame(pav.todo2[which(pav.todo2$sum==780),1])
write.table(core.gene, 'core_gene.txt', sep='\t', quote=F, row.names = F, col.names = F)
dis.gene=data.frame(pav.todo2[-which(pav.todo2$sum==780),1])
write.table(dis.gene, 'dis_gene.txt', sep='\t', quote=F, row.names = F, col.names = F)

temp=data.frame(ID=colnames(pav.todo)[-c(1:5)])
for(i in 6:ncol(pav.todo)) temp$Present[i-5]=length(which(pav.todo[,i]=='PRESENT'))
for(i in 6:ncol(pav.todo)) temp$core_Present[i-5]=length(which(pav.todo[which(pav.todo$Chr!=0),i]=='PRESENT'))
for(i in 6:ncol(pav.todo)) temp$pan_Present[i-5]=length(which(pav.todo[which(pav.todo$Chr==0),i]=='PRESENT'))
temp$ID=gsub('[.]','-',temp$ID)
temp$ID[1]='153_ML1628'

group=read.table('./780_mungbean_k8_Q0.7_group.txt', header=F, stringsAsFactors = F)
group=merge(group, temp, by.x='V1', by.y='ID')

pav.todo1=data.frame(t(pav.todo))
colnames(pav.todo1)=pav.todo1[1,]
pav.todo1=pav.todo1[-c(1:5),]
pav.todo1$ID=row.names(pav.todo1)
group=merge(group, pav.todo1, by.x='V1', by.y='ID')
write.table(group, 'pav_list_group.txt', sep='\t', quote=F, row.names = F)

group1=data.frame(t(group))
group1$Index=row.names(group1)
group1=group1[,c(ncol(group1),1:ncol(group1)-1)]
for(i in 2:ncol(group1)) group1[which(group1[,i]=='PRESENT'),i]=1
for(i in 2:ncol(group1)) group1[which(group1[,i]=='LOST'),i]=0
colnames(group1)[-1]=group1[2,-1]
group1=group1[-c(1:5),]
for(i in 2:ncol(group1)) group1[,i]=as.numeric(group1[,i])

rad.g=group1[,c(1,grep('Rad',colnames(group1)))]
rad.g$sum=apply(rad.g[,-1], 1, function(x) sum(x, na.rm=T))
sub.g=group1[,c(1,grep('Sub',colnames(group1)))]
sub.g$sum=apply(sub.g[,-1], 1, function(x) sum(x, na.rm=T))

rad.gene=rad.g[which(rad.g$sum!=0 & sub.g$sum==0),1]
write.table(data.frame(rad.gene), 'rad_specific_gene.txt', quote=F, row.names = F, col.names = F, sep='\t')

sub.gene=sub.g[which(rad.g$sum==0 & sub.g$sum!=0),1]
write.table(data.frame(sub.gene), 'sub_specific_gene.txt', quote=F, row.names = F, col.names = F, sep='\t')

rad1.g=group1[,c(1,grep('Rad1',colnames(group1)))]
rad1.g$sum=apply(rad1.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(rad1.g[1:33161,'sum']!=0))
length(which(rad1.g[33162:57067,'sum']!=0))
length(which(rad1.g[,'sum']!=0))

rad2.g=group1[,c(1,grep('Rad2',colnames(group1)))]
rad2.g$sum=apply(rad2.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(rad2.g[1:33161,'sum']!=0))
length(which(rad2.g[33162:57067,'sum']!=0))
length(which(rad2.g[,'sum']!=0))

subaf.g=group1[,c(1,grep('SubAF',colnames(group1)))]
subaf.g$sum=apply(subaf.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subaf.g[1:33161,'sum']!=0))
length(which(subaf.g[33162:57067,'sum']!=0))
length(which(subaf.g[,'sum']!=0))

subidn.g=group1[,c(1,grep('SubIDN',colnames(group1)))]
subidn.g$sum=apply(subidn.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subidn.g[1:33161,'sum']!=0))
length(which(subidn.g[33162:57067,'sum']!=0))
length(which(subidn.g[,'sum']!=0))

subaue.g=group1[,c(1,grep('SubAUe',colnames(group1)))]
subaue.g$sum=apply(subaue.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subaue.g[1:33161,'sum']!=0))
length(which(subaue.g[33162:57067,'sum']!=0))
length(which(subaue.g[,'sum']!=0))

subauw.g=group1[,c(1,grep('SubAUw',colnames(group1)))]
subauw.g$sum=apply(subauw.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subauw.g[1:33161,'sum']!=0))
length(which(subauw.g[33162:57067,'sum']!=0))
length(which(subauw.g[,'sum']!=0))

subassa.g=group1[,c(1,grep('SubAS_SA',colnames(group1)))]
subassa.g$sum=apply(subassa.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subassa.g[1:33161,'sum']!=0))
length(which(subassa.g[33162:57067,'sum']!=0))
length(which(subassa.g[,'sum']!=0))

subassea.g=group1[,c(1,grep('SubAS_SEA',colnames(group1)))]
subassea.g$sum=apply(subassea.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subassea.g[1:33161,'sum']!=0))
length(which(subassea.g[33162:57067,'sum']!=0))
length(which(subassea.g[,'sum']!=0))

subadmix_cul.g=group1[,c(1,grep('admix_cultivar',colnames(group1)))]
subadmix_cul.g$sum=apply(subadmix_cul.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subadmix_cul.g[1:33161,'sum']!=0))
length(which(subadmix_cul.g[33162:57067,'sum']!=0))
length(which(subadmix_cul.g[,'sum']!=0))

subadmix_wi.g=group1[,c(1,grep('admix_wild',colnames(group1)))]
subadmix_wi.g$sum=apply(subadmix_wi.g[,-1], 1, function(x) sum(x, na.rm=T))
length(which(subadmix_wi.g[1:33161,'sum']!=0))
length(which(subadmix_wi.g[33162:57067,'sum']!=0))
length(which(subadmix_wi.g[,'sum']!=0))

length(which(subaf.g$sum==0 & subassa.g$sum==0 & subassea.g$sum==0 & subaue.g$sum==0 & subauw.g$sum==0 & subidn.g$sum==0 & rad1.g$sum==0 & rad2.g$sum!=0))

write.table(data.frame(group1[which(rad2.g[,'sum']!=0),1]), 'Rad2_genes.txt', quote=F, row.names = F, col.names = F, sep='\t')

write.table(data.frame(setdiff(group1[which(rad1.g[,'sum']!=0),1],group1[which(rad2.g[,'sum']!=0),1])), 'Rad1_specific_genes.txt', quote=F, row.names = F, col.names = F, sep='\t')

write.table(data.frame(setdiff(group1[which(rad2.g[,'sum']!=0),1],group1[which(rad1.g[,'sum']!=0),1])), 'Rad2_specific_genes.txt', quote=F, row.names = F, col.names = F, sep='\t')

x1=read.table('./Group_gene_list/Rad1.txt', sep='\t', stringsAsFactors = F)
x2=read.table('./Group_gene_list/Rad2.txt', sep='\t', stringsAsFactors = F)
x4=unique(c(x1[,1],x2[,1]))
write.table(x4, 'Rad1_Rad2_genes.txt', sep='\t', quote=F, row.names = F, col.names = F)

x1=read.table('./Group_gene_list/SubAS_SA.txt', sep='\t', stringsAsFactors = F)
x2=read.table('./Group_gene_list/SubAS_SEA.txt', sep='\t', stringsAsFactors = F)
x4=unique(c(x1[,1],x2[,1]))
write.table(x4, 'SubAS_genes.txt', sep='\t', quote=F, row.names = F, col.names = F)

x1=read.table('./Group_gene_list/SubAUe.txt', sep='\t', stringsAsFactors = F)
x2=read.table('./Group_gene_list/SubAUw.txt', sep='\t', stringsAsFactors = F)
x4=unique(c(x1[,1],x2[,1]))
write.table(x4, 'SubAU_genes.txt', sep='\t', quote=F, row.names = F, col.names = F)


private=data.frame(pav1$GeneID)
for(i in 1:10){
  private1=matrix(ncol=2, nrow=nrow(pav1))
  for(j in 1:nrow(pav1)){
    private1[j,1]=length(which(pav1[j,which(colnames(pav1) %in% group[which(group$group==names(tab1)[i]),'ID'])]=='PRESENT'))
    private1[j,2]=length(which(pav1[j,which(colnames(pav1) %in% group[which(group$group!=names(tab1)[i]),'ID'])]=='LOST'))
  }
  colnames(private1)=paste0(names(tab1)[i], c('_PRESENT','_LOST'))
  private=cbind(private, private1)
}

#private_cry=private
#private_pan=private
head(private_cry)
length(which(private_cry$admix_cultivar_PRESENT==0 & private_cry$admix_wild_PRESENT==0 & private_cry$Rad1_PRESENT==0 & private_cry$Rad2_PRESENT==0 & private_cry$SubAF_PRESE==0 &private_cry$SubAS_SA_PRESENT==0 & private_cry$SubAS_SEA_PRESENT==0 & private_cry$SubAUe_PRESENT==0 & private_cry$SubAUw_PRESENT==0 & private_cry$SubIDN_PRESENT==0))

result=data.frame()
for(i in 6:ncol(pav2.1)){
  temp=c(colnames(pav2.1)[i],length(which(pav2.1[,i]=='LOST')),length(which(pav2.1[,i]=='PRESENT')))
  result=rbind(result, temp)
}
colnames(result)=c('Accession','LOST','PRESENT')
hist(as.numeric(result[,3]))

group=merge(group, result, by.x='V1', by.y='Accession')
colnames(group)=c('ID','group','Lost_ref','Present_ref','Lost_pan','Present_pan')
group$Lost_all=as.numeric(group$Lost_pan)+as.numeric(group$Lost_ref)
group$Present_all=as.numeric(group$Present_pan)+as.numeric(group$Present_ref)
for(i in 3:8) group[,i]=as.numeric(group[,i])
aggregate(group[,3:8], list(group[,2]), mean)


temp.group=data.frame(t(group[,-c(1:5)]))
temp=c()
for(i in 1:nrow(group)) temp[i]=length(table(temp.group[,i]))
table(temp)
grep('IC',group[,1])

for(i in 1:ncol(group)) group[which(group[,i]=='N'),i]='LOST'
write.table(group, 'todo_PAV_group_list.txt', sep='\t',quote=F,row.names = F)

library(ggplot2)
library(RColorBrewer)
colnames(group)[2]='group'
pdf('Present_all.pdf', width=6, height=4)
ggplot(group, aes(x=group, y=Present, color=group, fill=group)) +
  geom_boxplot()+
  scale_color_manual(values=c('grey','grey','skyblue','blue','yellow','peru','orange','red','#900603','pink'))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values=c('grey','grey','skyblue','blue','yellow','peru','orange','red','#900603','pink'))+
  theme_minimal()
dev.off()

#clean pan gff3
gff=read.delim('./mungbean_pan_plant_specific.gff3', sep='\t', stringsAsFactors = F, header=F)
gff=read.delim('./final_Pan_plant_specific_genes.gff3', sep='\t', stringsAsFactors = F, header=F)

gff$no=seq(nrow(gff))
gff$gene=unlist(lapply(strsplit(gff$V9, split = ';'),'[[',2))
gff$gene=gsub('Name=','',gff$gene)
gff$gene=gsub('Parent=','',gff$gene)
head(gff)

gff.1=gff[which(gff$gene %in% pav2.1$GeneID),]
gff.1=gff.1[order(gff.1$no),]
gff.1$V9=gsub(';;',';',gff.1$V9)
head(gff.1)
write.table(gff.1[,-c(10:11)], 'final_Pan_plant_specific_genes.gff3', sep='\t', row.names = F, col.names = F, quote=F)

#gene enrichment analysis
cry.gff=read.delim('~/Big_data/Mungbean/Reference/Crystal/Crystal_annotation.gff3', sep='\t', stringsAsFactors = F, header=F)

cry.gff.1=cry.gff[grep('GO:', cry.gff$V9),]
cry.gff.1$ID=unlist(lapply(strsplit(cry.gff.1$V9, split=';'),'[[',1))
cry.gff.1$ID=gsub('ID=','',cry.gff.1$ID)
cry.gff.1$GO=unlist(lapply(strsplit(cry.gff.1$V9, split=';GO:'),'[[',2))
cry.gff.1$GO=paste0('GO:',cry.gff.1$GO)

gff.1.1=gff.1[grep('GO:', gff.1$V9),]
gff.1.1$ID=unlist(lapply(strsplit(gff.1.1$V9, split=';'),'[[',1))
gff.1.1$ID=gsub('ID=','',gff.1.1$ID)
gff.1.1$GO=unlist(lapply(strsplit(gff.1.1$V9, split=';GO:'),'[[',2))
gff.1.1$GO=paste0('GO:',gff.1.1$GO)

go.enrich=data.frame(rbind(cry.gff.1[,10:11], gff.1.1[,12:13]))
write.table(go.enrich, 'GO_enrich_bakcground.txt', sep='\t', quote=F, row.names = F, col.names = F)

rad=c(private_cry[which(private_cry$Rad1_PRESENT>0 | private_cry$Rad2_PRESENT>0),1],private_pan[which(private_pan$Rad1_PRESENT>0 | private_pan$Rad2_PRESENT>0),1])
rad=rad[which(rad %in% go.enrich$ID)]

write.table(data.frame(rad), 'GO_rad.txt', sep='\t', quote=F, row.names = F, col.names = F)

Sub=c(private_cry[which(private_cry$SubAF_PRESENT>0 | private_cry$SubAS_SA_PRESENT>0 | private_cry$SubAS_SEA_PRESENT>0 | private_cry$SubAUe_PRESENT>0 | private_cry$SubAUw_PRESENT>0 | private_cry$SubIDN_PRESENT>0 ),1], private_pan[which(private_pan$SubAF_PRESENT>0 | private_pan$SubAS_SA_PRESENT>0 | private_pan$SubAS_SEA_PRESENT>0 | private_pan$SubAUe_PRESENT>0 | private_pan$SubAUw_PRESENT>0 | private_pan$SubIDN_PRESENT>0 ),1])
Sub=Sub[which(Sub %in% go.enrich$ID)]

write.table(data.frame(Sub), 'GO_sub.txt', sep='\t', quote=F, row.names = F, col.names = F)


clean=read.delim('./cleaned_pangenome_plant_annotation_table.txt', header=T, stringsAsFactors = F)
clean[1:5,]

library(seqinr)
fa=read.fasta('mungbean_pan_plant_specific.fa')
fa=read.fasta('mungbean_pan_plant_specific_inDB.fa')
names(fa)[1:5]
fa1=fa[which(names(fa) %in% clean[,'key'])]
write.fasta(fa1, 'mungbean_pan_plant_specific_inDB.fa', names=names(fa1))

gff=read.delim('./mungbean_pan_plant_specific.gff3', sep='\t', stringsAsFactors = F, header=F)
gff$length=gff$V5-gff$V4
gff$ID=unlist(lapply(strsplit(gff$V9, split=';'),'[[',1))
gff$ID=gsub('ID=','',gff$ID)

gff1=gff[which(gff$ID %in% clean[,1]),]
table(gff1$V3)
hist(gff[which(gff1$V3=='gene'),'length'])
gff1[which(gff1$V3=='gene' & gff1$length > 30000),]
