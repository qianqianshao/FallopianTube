# 5.11.2020 by Qianyi
# individual clustering for FT1


#1 batch: 3 segments of fallopian tube (FT1) from 1 peri-menopausal women
#(1) Fimbria
#(2) Ampulla 
#(3) Isthmus
#directly-merge the 3 segments of fallopian tube (FT1) together and cluster

library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

redblue100<-rgb(read.table(paste0('data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

dataset=c("747-YS-1","747-YS-2","747-YS-3")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1")
run=c("NovaA-223","NovaA-223","NovaA-223")
sample=c("747-YS","747-YS","747-YS")
part=c("Fimbria","Ampulla","Isthmus")
organ=c("FallopianTube","FallopianTube","FallopianTube")
menopause=c("pre","pre","pre")
subject=c("Human1","Human1","Human1")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)
subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube")

### load directly-merged 3 segments of fallopian tube1
indiv=1
  load(file=paste0(subjectorgans[1],".Robj"))


########## Analysis for merged data of individual organ
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))

###### Determine top PCs
dgelist=list()
numPCs=10
i=1

pdf(paste("plot/dge_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
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


pdf(paste("plot/dge_PCA_Variablel_variation_topPCs_heatmap.pdf",sep=""),height=12,width=8)
DimHeatmap(dge, dims = 1:12,ncol=4,cells=1000,  balanced = TRUE)
dev.off()


###### Using top PCs
### tSNE
	dge <- RunTSNE(dge, dims = 1:numPCs[i])
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,2,by=0.1))
dgeall=dge


print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))


###### order clusters for each individual organ
res=paste0("RNA_snn_res.0.",c(6,4))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
for(resi in 1:2){
Idents(dge)<-dge[[res[resi]]]
print(c(dataset[resi],length(unique(Idents(dge)))))
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

### get order of seriation
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]

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
}

### save dge with ordered clusters
dgeall=dge;i=1
save(dgeall,file=paste0(subjectorgans[1],".Robj"))



### Comparison with clusters of individual segments
ncellsindivclusters=table(Idents(dge),dge@meta.data$indivclusters)
  write.table(ncellsindivclusters,paste0("plot/",numPCs[i],"ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")


### flip middle clusters in order to make its order consistent with seriated order of Fimbria1 and Isthmus1
# 10PCs, 15 ordered clusters, cluster 9->5, 8->6, 5->8, 6->9
i=2
dge=dgelist[[i]]
Idents(dge)<-dge[["RNA_snn_res.0.4ordered"]]
ident<-as.numeric(Idents(dge))
names(ident)=names(Idents(dge))
  table(ident)
  ident[which(ident==9)]<-55
  ident[which(ident==8)]<-66
  ident[which(ident==5)]<-8
  ident[which(ident==6)]<-9
  ident[which(ident==55)]<-5
  ident[which(ident==66)]<-6
  ident=factor(ident,levels=1:15)
  table(ident)
  dge[["ordered"]]<-ident
Idents(dge)<-dge[["ordered"]]
  table(Idents(dge))
  dgeall=dge
dgelist[[i]]=dge
save(dgeall,file=paste0(subjectorgans[1],".Robj"))






### PCA and tSNE for ordered clusters
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


### Comparison with clusters of individual segments
ncellsindivclusters=table(Idents(dge),dge@meta.data$indivclusters)
  write.table(ncellsindivclusters,paste0("plot/",numPCs[i],"PCs_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")

### Contribution of each batch to each cluster
for(i in 1:2){
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(Idents(dge),dge@meta.data$Part)
  print(ncellsbatchcluster)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,1)
  print(percentcellsbatch)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,2)
  print(percentcellscluster)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0("plot/",numPCs[i],"PCs_ncellspercluster.txt"),quote=F,row.names=T,col.names=T,sep="\t")
}



### Visualize individual time and treatment
label="Part"
pdf(paste0("plot/label_",label,".pdf"),height=6,width=8.5)
plotlist=list()
for(i in 1:2){
	Idents(dgeall)<-dgeall[[label]]
plotlist[[1]]=PCAPlot(dgeall,pt.size=1,label=TRUE,label.size=6)
plotlist[[2]]=PCAPlot(dgeall,c(1,3),pt.size=1,label=TRUE,label.size=6)
plotlist[[3]]=TSNEPlot(dgeall,pt.size=1,label=TRUE,label.size=6)
plotlist[[4]]=DimPlot(dgeall,pt.size=1,reduction="umap",label=TRUE,label.size=6)
multiplot(plotlist,cols = 2)
}
dev.off()
### Visualize individual batch and subject
for(label in c("Part")){
plotlist=list()
for(indiv in 1:2){
dgeall=dgealllist[[indiv]]
Idents(dgeall)<-dgeall[[label]]
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlim=ylim=list()
pos=c("topleft","topright","bottomright","bottomright")
xlim=list(range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"umap")[,1]),range(Embeddings(dge,"tsne")[,1]))
ylim=list(range(Embeddings(dge,"pca")[,2]),range(Embeddings(dge,"pca")[,3]),range(Embeddings(dge,"umap")[,2]),range(Embeddings(dge,"tsne")[,2]))
### plot PCs and tSNE for each batch using the other batches as background
dge=dgeall
sets=levels(Idents(dgeall))
if(length(sets)>1){
  cols=gg_color_hue(length(sets))
plot2set=plot3set=plot4set=plottset=NULL
dim=list(Embeddings(dge,"pca")[,1:2],Embeddings(dge,"pca")[,c(1,3)],Embeddings(dge,"umap"),Embeddings(dge,"tsne"))
size=length(sets)
if(size>4){
  pdf(paste0("plot/",subjectorgans[1],"_",label,"_indiv.pdf"),height=2.3*round(sqrt(size)),width=2.3*ceiling(sqrt(size)))
  par(mfrow=c(round(sqrt(size)),ceiling(sqrt(size))),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
} else {
  pdf(paste0("plot/",subjectorgans[1],"_",label,"_indiv.pdf"),height=2.5,width=2.5*size)
  par(mfrow=c(1,size),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
}
for(j in 1:4){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(Idents(dgeall))
names(ident)=names(Idents(dgeall))
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(cols[seti],0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
if(size>(round(sqrt(size))*ceiling(sqrt(size)))){
  for(new in (round(sqrt(size))*ceiling(sqrt(size))+1):size){
    plot.new()
  }
}
}
dev.off()
}
}
}



###### Differentially-expressed markers
res="ordered"
markerslist=list()
for(i in 1:2){
	dge=dgelist[[i]]
	Idents(dge)<-dge[["ordered"]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
markerslist[[i]]=markers
write.table(markers,paste0("plot/",subjectorgans[1],"_",numPCs[i],"PCs_",res,"_mindiff0.2_logfc2fold_5.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
}


markerslist=list()
for(i in 1:2){
markers=read.table(paste0("plot/",subjectorgans[1],"_",numPCs[i],"PCs_",res,"_mindiff0.2_logfc2fold_5.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markerslist[[i]]=markers
}

plotlist=list()
for(i in 1:2){
dge=dgelist[[i]]
markers=markerslist[[i]]
print(table(Idents(dge)))
  print(i)
  print(table(markers$cluster))
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
size=sqrt(length(top2$gene))
p <- FeaturePlot(dge, top2$gene, min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p)
}

pdf(paste0("plot/markerstop.pdf"),height=2*round(size),width=1.8*ceiling(size))
plotlist
dev.off()


for(i in 1:2){
pdf(file=paste0("plot/",numPCs[i],"PCs_clusters_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=12)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
}
dev.off()

i=2
plotlist=list()
for(j in 1:3){
dgetmp<-subset(dge,Part==part[j])
print(table(Idents(dgetmp)))
plotlist[[j]]=VlnPlot(dgetmp, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3+j]]=VlnPlot(dgetmp, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[6+j]]=VlnPlot(dgetmp, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
}
pdf(file=paste0("plot/",numPCs[i],"PCs_clusters_Parts_PerCellAttributes_ViolinPlot.pdf"),height=7.5,width=12)
multiplot(plotlist,cols = 3)
dev.off()


# Differential gene expression compared across the fallopian tube segments
  Idents(dge)<-dge[["Part"]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
print(table(Idents(dge)))
  print(table(markers$cluster))
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.5),min.pct=0)
  print(table(markers$cluster))
write.table(markers,paste0("plot/",subjectorgans[1],"_3segments_min0_logfc1.5fold_5.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")

pdf(file=paste0("plot/3Segments_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=7.5)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
dev.off()


# FT1:
#-> 15 re-ordered clusters, flipped the order to be consistent with the order of Fimbria
#subclustering for secretory cluster 2-4, epithelial cluster 1-4, cluster 10, and cluster 13
subsets=list(2:4,1:4,10,13)
subsetsname=c("2-4","1-4",10,13)

for(i in 1:length(subset)){

cc=subsets[[i]]
ccname=subsetsname[i]

dgefile=paste0("plot/C",ccname,"_")

dge=subset(dgeall,ordered %in% cc)
table(Idents(dgeall))
table(Idents(dge))
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))

### Highly variable genes
dge<-FindVariableFeatures(dge)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge), 10)
plot1 <- VariableFeaturePlot(dge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000

pdf(paste0(dgefile,"HVG_x0.2_y0.2.pdf"),height=6,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
print(plot2)
dev.off()

### Scale data 
dge <- ScaleData(dge,features=rownames(dge))
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)

###### Determine top PCs
numPCs=c(5,6,7,7)

pdf(paste(dgefile,"dge_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
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


pdf(paste(dgefile,"dge_PCA_Variablel_variation_topPCs_heatmap.pdf",sep=""),height=12,width=8)
DimHeatmap(dge, dims = 1:12,ncol=4,cells=1000,  balanced = TRUE)
dev.off()


###### Using top PCs
### tSNE
  dge <- RunTSNE(dge, dims = 1:numPCs[i],perplexity=5)
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,2,by=0.1))


print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))



###### order clusters for each individual organ
res=paste0("RNA_snn_res.0.",c(3,1,6,6))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
resi=i
Idents(dge)<-dge[[res[resi]]]
print(c(dataset[resi],length(unique(Idents(dge)))))
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
### get order of seriation
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]

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
dgelist[[resi]]=dge


### save dge with ordered clusters
save(dge,file=paste0(subjectorgans[1],"_C",ccname,".Robj"))
}

pdf(paste0(dgefile,"clusters_new1.pdf"),height=8,width=18)
for(label in c("ordered","new1","RNA_snn_res.0.3ordered")){
Idents(dge)<-dge[[label]]
plotlist=list()
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=FALSE)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=FALSE)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=FALSE)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=FALSE)
for(j in 4:7){
plotlist[[j+1]]=PCAPlot(dge,c(1,j),pt.size=1,label=FALSE)
}
multiplot(plotlist,cols = 4)
}
dev.off()


###### Differentially-expressed markers
res=c("RNA_snn_res.0.3ordered","RNA_snn_res.0.1ordered","RNA_snn_res.0.6ordered","RNA_snn_res.0.6ordered")
markerslist=list()
for(i in 1:2){
  dge=dgelist[[i]]
  Idents(dge)<-dge[[res[i]]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
markerslist[[i]]=markers
write.table(markers,paste0(dgefile,"",subjectorgans[1],"_",numPCs[i],"PCs_",res[i],"_mindiff0.2_logfc2fold_5.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
}

markerslist=list()
for(i in 1:2){
markers=read.table(paste0(dgefile,"",subjectorgans[1],"_",numPCs[i],"PCs_",res[i],"_mindiff0.2_logfc2fold_5.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markerslist[[i]]=markers
}

plotlist=list()
for(i in 1:2){
dge=dgelist[[i]]
markers=markerslist[[i]]
print(table(Idents(dge)))
  print(i)
  print(table(markers$cluster))
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
size=sqrt(length(top2$gene))
p <- FeaturePlot(dge, top2$gene, min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p)
}

pdf(paste0(dgefile,"markerstop.pdf"),height=2*round(size),width=1.8*ceiling(size))
plotlist
dev.off()


pdf(file=paste0(dgefile,"",numPCs[i],"PCs_clusters_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=6)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
dev.off()

# 6.3.2020 notes
### For cluster 13 subclustering, remove the filter of >20% difference in detection rate
i=4
  dge=dgelist[[i]]
  Idents(dge)<-dge[[res[i]]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.pct=0)
  print(table(markers$cluster))
write.table(markers,paste0(dgefile,"_",numPCs[i],"PCs_",res[i],"_logfc2fold_minpct0_5.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.6),min.pct=0)
  print(table(markers$cluster))
write.table(markers,paste0(dgefile,"_",numPCs[i],"PCs_",res[i],"_logfc1.6fold_minpct0_5.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")


### For Clusters1-4 reclustering, check markers 4 Vs 5
## Differential expression between cluster 4 and 5
## Differential expression for cluster 4 and 5 Vs All others
i=2
dge=dgelist[[i]]
cc=subsets[[i]]
ccname=subsetsname[i]
dgefile=paste0("plot/C",ccname,"_")
## Differential expression for cluster 4 and 5 Vs All others
markers=FindMarkers(dge,c(4,5),test.use="bimod",logfc.threshold = log(2),only.pos=TRUE,min.diff.pct=0.2)
dim(markers[markers$avg_logFC>0,]) #0
markers=FindMarkers(dge,c(4,5),test.use="bimod",logfc.threshold = log(2),only.pos=TRUE,min.diff.pct=0)
dim(markers[markers$avg_logFC>0,]) #2
markers=FindMarkers(dge,c(4,5),test.use="bimod",logfc.threshold = log(1.6),only.pos=TRUE,min.diff.pct=0.2)
dim(markers[markers$avg_logFC>0,]) #5
markers=FindMarkers(dge,c(4,5),test.use="bimod",logfc.threshold = log(1.6),only.pos=TRUE,min.diff.pct=0)
dim(markers[markers$avg_logFC>0,]) #14
write.table(markers,paste0(dgefile,"6clusters_Cluster45VsAllOthers_logfc1.6fold_7.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers$gene=rownames(markers)
# visualize markers
markers.use=rownames(markers)
size=sqrt(length(markers.use))
pdf(paste0(dgefile,"markers_Cluster45VsAllOthers_VlnPlot.pdf"),height=2*round(size),width=2.5*ceiling(size))
VlnPlot(dge,markers.use,ncol=ceiling(size),pt.size=-1)
dev.off()
pdf(paste0(dgefile,"markers_Cluster45VsAllOthers_heatmap.pdf"),height=8,width=10)
DoHeatmap(dge,features=rownames(markers))
dev.off()
pdf(paste0(dgefile,"markers_Cluster45VsAllOthers_Feature.pdf"),height=1.5*round(size),width=1.5*ceiling(size))
p <- FeaturePlot(dge,markers.use,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
print(plotlist)
dev.off()

## Differential expression between cluster 4 and 5
markers=FindMarkers(dge,4,5,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2)
dim(markers[markers$avg_logFC>0,]) #0
dim(markers[markers$avg_logFC<0,]) #5
markers=FindMarkers(dge,4,5,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0)
dim(markers[markers$avg_logFC>0,]) #2
dim(markers[markers$avg_logFC<0,]) #5
markers=FindMarkers(dge,4,5,test.use="bimod",logfc.threshold = log(1.6),min.diff.pct=0)
dim(markers[markers$avg_logFC>0,]) #3
dim(markers[markers$avg_logFC<0,]) #39
write.table(markers,paste0(dgefile,"6clusters_Cluster4Vs5_logfc1.6fold_7.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindMarkers(dge,4,5,test.use="bimod",logfc.threshold = log(1.6),min.diff.pct=0,min.pct=0)
dim(markers[markers$avg_logFC>0,]) #8
dim(markers[markers$avg_logFC<0,]) #40
write.table(markers,paste0(dgefile,"6clusters_Cluster4Vs5_logfc1.6fold_min0pct_7.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers$gene=rownames(markers)
### visualize markers
topp=markers[markers$avg_logFC>0,]
topn=markers[markers$avg_logFC<0,]
markers.use=c(topp$gene[1:8],topn$gene)
size=sqrt(length(markers.use))
pdf(paste0(dgefile,"markers_Cluster4Vs5_VlnPlot.pdf"),height=2*round(size),width=2.5*ceiling(size))
VlnPlot(dge,markers.use,ncol=ceiling(size),pt.size=-1)
dev.off()
pdf(paste0(dgefile,"markers_Cluster4Vs5_heatmap.pdf"),height=8,width=10)
DoHeatmap(dge,features=c(topp$gene,topn$gene))
dev.off()
pdf(paste0(dgefile,"markers_Cluster4Vs5_Feature.pdf"),height=1.5*round(size),width=1.5*ceiling(size))
p <- FeaturePlot(dge,markers.use,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
print(plotlist)
dev.off()
markers=FindMarkers(dge,4,5,test.use="bimod",logfc.threshold = -Inf,min.diff.pct=-Inf,min.pct=-Inf,min.cells.feature = -Inf,  min.cells.group = -Inf)
dim(markers) #8
write.table(markers,paste0(dgefile,"6clusters_Cluster4Vs5_allgenes_7.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
write.table(markers[,c(1,2)],paste0(dgefile,"6clusters_Cluster4Vs5_allgenes_LRpath_input.txt"),col.names=F,row.names=T,quote=F,sep="\t")




# 10PCs, 15 ordered clusters, incorporate ReclusteredC1-4 -> 1-6, subclsuteredC10 -> 1-3
dge=dgeall
subsets=list(1:4,10,13)
subsetsname=c("1-4",10,13)
### load subclustering for merged FT1 - 15 re-ordered clusters with 10 PCs
dgelist=list()
for(i in 1:length(subsets)){
cc=subsets[[i]]
ccname=subsetsname[i]
load(file=paste0(subjectorgans[1],"_C",ccname,".Robj"))
dgelist[[i]]=dge
}
ident<-as.numeric(Idents(dgeall))
names(ident)=names(Idents(dgeall))

id1=paste(subsetsname[1],as.numeric(Idents(dgelist[[1]])),sep="_")
id2=paste(subsetsname[2],as.numeric(Idents(dgelist[[2]])),sep="_")
names(id1)=names(Idents(dgelist[[1]]))
names(id2)=names(Idents(dgelist[[2]]))
id=c(id1,id2)
table(ident)
  table(ident[names(id)])
keep=ident[which(!(names(ident) %in% names(id)))]
identnew=c(id,keep)
  identnew=factor(identnew,levels=c(paste("1-4",1:6,sep="_"),5:9,paste("10",1:3,sep="_"),11:15))
  table(identnew)
1-4_1 1-4_2 1-4_3 1-4_4 1-4_5 1-4_6     5     6     7     8     9  10_1  10_2 
  107    31   123  1089   835    44  3512   783   625   752   214    73    27 
 10_3    11    12    13    14    15 
   64  1164   741    78   168    97
dge=dgeall
  dge[["new1"]]<-identnew
Idents(dge)<-dge[["new1"]]
  table(Idents(dge))
  dgeall=dge
  dge2=dge
dgelist[[i]]=dge
save(dgeall,file=paste0(subjectorgans[1],".Robj"))



