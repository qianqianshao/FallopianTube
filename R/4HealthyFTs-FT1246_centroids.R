# 7.9.2021 by Qianyi
# compared rank correlation of cluster centroids across 4 healthy FTs: FT1, FT2, FT4 and FT6
# Related to Figure S1C.

FT1: 3 segments of fallopian tube from 1 healthy perimenopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT2: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT4: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla
(3) Isthmus
FT6: 1 whole fallopian tube from Sample3 (2518-AJ)  - included FT6, ovary, 2 Uterus parts from 1 caucasian woman, sequenced in 1 batch
for each fallopian tube, did clustering for each whole organ
merge the 4 healthy individual fallopian tubes together and cross-tabulate the individual cluster centroids
make the clustering order of 4 healthy individual organs to be consistent with FT1


R 

library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

redblue100<-rgb(read.table(paste0('data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette1=gg_color_hue(4)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3","2518-AJ-1","2518-AJ-2","2518-AJ-3","2518-AJ-4","2788-NU-1","2788-NU-2")
dataset1=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3","2518-AJ-1_CACTACGA-ATCAGTCT","2518-AJ-2_CACGGTGA-TGTGACGA","2518-AJ-3_ATGGCTTG-CACAACAT","2518-AJ-4_CCTTCTAG-TCGTTGTA","2788-NU-1_TATCAGCC-AGGACGAA","2788-NU-2_TGGTCCCA-ACGCCAGA")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium1","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3","Fimbria4","Ampulla4","Isthmus4","Uterus2","FT5","Ovary1","FT6","Ovary2","Myometrium2","Endometrium2","Myometrium3","Endometrium3")
run=c("NovaA-223","NovaA-223","NovaA-223","NovaA-237","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-295","NovaA-295","NovaA-295","NovaA-301","NovaA-301","NovaA-301","NovaZ-8","NovaZ-8","NovaZ-8","NovaZ-8","2788-NU","2788-NU")
sample=c("747-YS","747-YS","747-YS","1073-NU","1427-NU","1427-NU","1427-NU","1407-NU","1407-NU","1714-YS","1714-YS","1714-YS","1847-NU","1847-NU","1847-NU","2518-AJ","2518-AJ","2518-AJ","2518-AJ","2788-NU","2788-NU")
part=c("Fimbria","Ampulla","Isthmus","Myometrium","Fimbria","Ampulla","Isthmus","FimAmp","Isthmus","Fimbria","Ampulla","Isthmus","Uterus","FT","Ovary","FT","Ovary","Myometrium","Endometrium","Myometrium","Endometrium")
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus","FallopianTube2","FallopianTube2","FallopianTube2","FallopianTube3","FallopianTube3","FallopianTube4","FallopianTube4","FallopianTube4","Uterus2","FallopianTube5","Ovary1","FT6","Ovary2","Sample3U2","Sample3U2","Uterus3","Uterus3")
menopause=c("peri","peri","peri","post","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","unknown","unknown")
subject=c("Human1","Human1","Human1","Human2","Human3","Human3","Human3","Human4","Human4","Human5","Human5","Human5","Human6","Human6","Human6","Sample3","Sample3","Sample3","Sample3","Uterus3","Uterus3")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)

subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube","Uterus","FallopianTube2","FallopianTube3","FallopianTube4","Uterus2","FallopianTube5","Ovary1","FallopianTube6")

ft=grep("Fallopian|FT",organ)
indivft=grep("Fallopian|FT",organs)


healthyft=indivft[c(1,2,4,6)]
all="FallopianTube1246"



### load object for each of 4 healthy FTs 
dgealllist=list()
for(indiv in healthyft){
  load(file=paste0(subjectorgans[indiv],".Robj")) 
dge=dgeall
dgealllist[[indiv]]=dge
}

i=16
  load(file=paste0(dataset[i],".Robj"))
dgealllist[[healthyft[4]]]=dge
dgelist=dgealllist




# merge all samples for healthy fallopian tubes (all subjects and all parts)
CheckDuplicateCellNames(dgelist)


for(i in healthyft){
print(table(gsub("_.*","",rownames(dgelist[[i]]@meta.data))))
}

dge1=dgelist[[1]]
for(i in healthyft[-c(1,4)]){
dge12=dgelist[[i]]
dge2=merge(dge1,dge12)
dge1=dge2
}
dge12=dgelist[[healthyft[4]]]
dge2=merge(dge1,dge12)
dgeall=dge2

save(dgeall,file=paste0(all,".Robj"))
table(gsub("_.*","",names(Idents(dge))))





## Rank cor

### using either union of HVG or top 50 markers for each cluster
res=c(paste0("RNA_snn_res.0.",c(4,0,2,2,3,0,0,0,6),"ordered") )
reps=healthyft

hvg.union=NULL
for(i in healthyft){
hvg.union=unique(c(hvg.union,VariableFeatures(dgelist[[i]])))
}
length(hvg.union) # 3623

top50markers=NULL
for(i in healthyft){
markers=read.table(paste0("plot/",subjectorgans[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markers=read.table(paste0("plot/",organ[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markers %>% group_by(cluster) %>% top_n(50, avg_logFC)  -> top50
top50markers=unique(c(top50markers,unique(top50$gene)))
}
length(top50markers) #923

genelist=list(hvg.union,top50markers)
genelabels=c("HVG","Top50Markers")

genes=hvg.union
genelabel="HVG"

### order cells by batch first, then by clusters of each batch
dge=dgelist[[1]]
Idents(dge) <- dge$ordered
dgelist[[1]] <- dge
dgealllist=dgelist

blockident=NULL
for(i in healthyft){
  tmp=paste(organs[i],Idents(dgelist[[i]]),sep="_")
  names(tmp)=names(Idents(dgelist[[i]]))
  blockident=c(blockident,tmp)
}
blockident=blockident[names(Idents(dgeall))]

### Clusters ordered first by batches, then by res
batch=organs[healthyft]
nbatch=length(batch)
ncluster=NULL
for(i in healthyft){
  ncluster=c(ncluster,length(unique(Idents(dgelist[[i]]))))
}
ncluster 
clusters=list()
for(i in 1:nbatch){
  clusters[[i]]=rep(1:ncluster[i])
}
levels2=NULL
for(bb in 1:nbatch){
    cluster=clusters[[bb]]
    levels2=c(levels2,paste(batch[bb],cluster,sep="_"))
}
levels2=levels2[which(levels2 %in% unique(blockident))]
levels=levels2

ident=factor(blockident,levels=levels)

### Calculate correlation for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
dge=dgeall
Idents(dge)<-ident
dge$healthyftclusters<-ident

dgeall=dge
table(ident)
save(dgeall,file=paste0(all,".Robj"))




centroid=log(AverageExpression(dge)$RNA+1)
write.table(centroid,paste0("plot/",all,"_healthyft_Centroid.txt"),row.names=T,col.names=T,quote=F,sep="\t")


pdf(file=paste("plot/organs_indiv_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=7.5,width=7)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))

for(g in 1:length(genelist)){
genelabel=genelabels[g]
genes=genelist[[g]]
print(length(genes))

cc=cor(as.matrix(centroid)[genes,],method="spearman")
dim(cc)
min(cc) # 0.036

data.use=cc[levels,levels]
write.table(data.use,paste0("plot/",all,"_healthyft_Centroid_rho_",genelabel,".txt"),row.names=T,col.names=T,quote=F,sep="\t")

### load cluster centroid rank correlation using HVG
data.use=read.table(paste0("plot/",all,"_healthyft_Centroid_rho_",genelabel,".txt"),header=T,row.names=1)
colnames(data.use)=rownames(data.use)
levels=rownames(data.use)
batch=organs[healthyft]

### labeling
colsep.use=cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])-table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))]/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

sidecol=do.call(rbind,strsplit(levels,"_"))
batchtmp=batch[which(batch %in% unique(sidecol[,1]))]
for(rep in 1:length(unique(sidecol[,1]))){
a=batchtmp[rep]
sidecol[which(sidecol[,1]==a),1]<-rep
}

rlab=matrix(0,2,length(levels))
rlab[1,]=rep(c("white","black"),6)[as.numeric(sidecol[,1])]
for(i in 1:nrow(sidecol)){
  rlab[2,i]=myBrewerPalette[as.numeric(sidecol[i,2])]
}
clab=cbind(rlab[2,],rlab[1,])
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

col.use=redblue100

heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.5,cexRow=.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))


}
dev.off()
# saved as Figure S1C.