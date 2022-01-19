# 8.28.2020 by Qianyi
# individual clustering for FT2
# individual clustering for FT3

FT1: 3 segments of fallopian tube from 1 healthy perimenopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT2: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT3: 2 segments of fallopian tube from 1 pre-menopausal women with disease, sequenced in 1 batch
(1) Fimbria and Ampulla
(2) Isthmus
did individual clustering for each segment
for each fallopian tube, directly-merge the segments together and do clustering

R


library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

redblue100<-rgb(read.table(paste0('data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium1","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3")
run=c("NovaA-223","NovaA-223","NovaA-223","NovaA-237","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-277")
sample=c("747-YS","747-YS","747-YS","1073-NU","1427-NU","1427-NU","1427-NU","1407-NU","1407-NU")
part=c("Fimbria","Ampulla","Isthmus","Myometrium","Fimbria","Ampulla","Isthmus","FimAmp","Isthmus")
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus","FallopianTube2","FallopianTube2","FallopianTube2","FallopianTube3","FallopianTube3")
menopause=c("peri","peri","peri","post","pre","pre","pre","pre","pre")
subject=c("Human1","Human1","Human1","Human2","Human3","Human3","Human3","Human4","Human4")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)

subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube","Uterus","FallopianTube2","FallopianTube3")

ft=grep("Fallopian",organ)
indivft=grep("Fallopian",organs)
all="FallopianTube123"


#for each of the fallopian tubes, directly-merge the segments together
# merge all parts for each individual organ
dgealllist=list()

plot2=list()


for(indiv in 3:length(organs)){
reps=grep(organs[indiv],organ)
if(length(reps)==1){
  dgeall=dgelist[[reps]]
} else {

dge1=RenameCells(dgelist[[reps[1]]], add.cell.id = names[reps[1]])
for(i in 2:length(reps)){
dge12=RenameCells(dgelist[[reps[i]]], add.cell.id = names[reps[i]])
dge2=merge(dge1,dge12)
dge1=dge2
}
dge=dge1


### Highly variable genes
dge<-FindVariableFeatures(dge)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge), 10)
plot1 <- VariableFeaturePlot(dge)
plot2[[indiv]] <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000
### Scale data 
dge <- ScaleData(dge,features=rownames(dge))
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)
dgeall=dge

}



dgealllist[[indiv]]=dgeall
save(dgeall,file=paste0(subjectorgans[indiv],".Robj"))
print(table(gsub("_.*","",names(Idents(dge)))))
}

pdf(paste0("plot/organs_HVG_x0.2_y0.2.pdf"),height=6,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(plot2,cols=2)
dev.off()

for(indiv in 3:length(organs)){

dgeall=dgealllist[[indiv]]
### order cells by batch first, then by clusters of each batch
blockident=NULL
for(i in 1:length(grep(organs[indiv],organ))){
  rep=grep(organs[indiv],organ)[i]
  tmp=paste(part[rep],Idents(dgelist[[rep]]),sep="_")
  names(tmp)=paste(names[rep],names(Idents(dgelist[[rep]])),sep="_") # paste(dataset[rep],,sep="_")
  blockident=c(blockident,tmp)
}
blockident=blockident[names(Idents(dgeall))]

### Clusters ordered first by batches, then by res
batch=part[grep(organs[indiv],organ)]
nbatch=length(batch)
ncluster=NULL
for(i in 1:length(grep(organs[indiv],organ))){
  rep=grep(organs[indiv],organ)[i]
  ncluster=c(ncluster,length(unique(Idents(dgelist[[rep]]))))
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

dgeall=dgealllist[[indiv]]
dge=dgeall
Idents(dge)<-ident
dge$indivclusters<-ident

dgeall=dge
dgealllist[[indiv]]=dge


save(dgeall,file=paste0(subjectorgans[indiv],".Robj"))
}



for(i in indivft[-1]){
dge=dgealllist[[i]]
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))
}

###### Determine top PCs
numPCs=c(10,0,12,18)
for(i in indivft[-1]){
dge=dgealllist[[i]]
pdf(paste("plot/",organs[i],"_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(Stdev(dge,reduction="pca")[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
legend("topright",legend=organs[1],cex=1.5)
abline(v=numPCs[i]+0.5,lty=2)
### density plot of Eigenvalue
eigenvalue=Stdev(dge,reduction="pca")[numPCs[i]]
print(eigenvalue)
plot(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(Stdev(dge,reduction="pca")),col="black")
lines(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.55,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=organs[1],cex=1.5)
dev.off()


pdf(paste("plot/",organs[i],"_PCA_Variablel_variation_topPCs_heatmap.pdf",sep=""),height=15,width=10)
DimHeatmap(dge, dims = 1:20,ncol=5,cells=1000,  balanced = TRUE)
dev.off()


###### Using top PCs
### tSNE
	dge <- RunTSNE(dge, dims = 1:numPCs[i])
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,2,by=0.1))
dgeall=dge

dgealllist[[i]]=dgeall

print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))
}


###### order clusters for each individual organ
res=paste0("RNA_snn_res.0.",c(4,0,2,2))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()

for(resi in indivft[-1]){

dge=dgealllist[[resi]]

Idents(dge)<-dge[[res[resi]]]
print(c(organs[resi],length(unique(Idents(dge)))))
levels=levels(Idents(dge))

centroid=log(AverageExpression(dge)$RNA+1)


### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=centroid[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]
if(resi==3){
  levels=levelss[[resi]][c(1,5:9,13:15,10:12,4:2)]
}
if(resi==4){
  levels=levelss[[resi]][c(1,8:9,2:7,16:13,10,11:12)]
}

### Reordered clusters for all cells
cells.use=colnames(dge)
# random shuffling cells within ordered clusters
ident=factor(Idents(dge),levels=levels-1)

cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

### save ordered cluster ID in dge object
which(unique(cells.ident.ordered)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
# save the dge file
dgelist[[resi]]=dge
dgeall=dge
dgealllist[[resi]]=dgeall

}

### save dge with ordered clusters
for(i in indivft[-1]){
  dgeall=dgealllist[[i]]
save(dgeall,file=paste0(subjectorgans[i],".Robj"))
}


### Comparison with clusters of individual segments
for(i in indivft[-1]){
  dge=dgealllist[[i]]
ncellsindivclusters=table(Idents(dge),dge@meta.data$indivclusters)
  write.table(ncellsindivclusters,paste0("plot/",organs[i],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")




### PCA and tSNE for ordered clusters of each individual replicate
plotlist=list()
i=1
dge=dgelist[[i]]
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0.pdf",height=8,width=9)
multiplot(plotlist,cols = 2)
dev.off()


