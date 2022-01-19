# 11.30.2020 by Qianyi
# Directly-merged 3 large healthy FTs - FT1, FT2 and FT4 and did clustering
# not used in final analysis

FT1: 3 segments of fallopian tube from 1 healthy perimenopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla - where the ovum is fertilized.
(3) Isthmus
FT2: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla - where the ovum is fertilized.
(3) Isthmus
FT3: 2 segments of fallopian tube from 1 pre-menopausal women with disease, sequenced in 1 batch
(1) Fimbria and Ampulla - where the ovum is fertilized.
(2) Isthmus
FT4: 1 batch: 3 segments of fallopian tube from 1 pre-menopausal women, healthy FT, large fibroids
(1) Fimbria
(2) Ampulla - where the ovum is fertilized.
(3) Isthmus
FT5: 1 batch: merged segments of fallopian tube from 1 pre-menopausal women, healthy FT
did individual clustering for each segment
for each of the 5 fallopian tubes, directly-merged the segments together and did clustering
FT3 is with disease, and FT5 has very few cells, so we decided to combine 3 healthy FTs - FT1, FT2, and FT4. 
merge 3 FTs - FT1, FT2, FT4 (9 datasets) together and do clustering
1. directly-merge all 3 healthy fallopian tubes (9 datasets) together and do clustering


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

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium1","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3","Fimbria4","Ampulla4","Isthmus4","Uterus2","FT5","Ovary1")
run=c("NovaA-223","NovaA-223","NovaA-223","NovaA-237","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-295","NovaA-295","NovaA-295","NovaA-301","NovaA-301","NovaA-301")
sample=c("747-YS","747-YS","747-YS","1073-NU","1427-NU","1427-NU","1427-NU","1407-NU","1407-NU","1714-YS","1714-YS","1714-YS","1847-NU","1847-NU","1847-NU")
part=c("Fimbria","Ampulla","Isthmus","Myometrium","Fimbria","Ampulla","Isthmus","FimAmp","Isthmus","Fimbria","Ampulla","Isthmus","Uterus","FT","Ovary")
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus","FallopianTube2","FallopianTube2","FallopianTube2","FallopianTube3","FallopianTube3","FallopianTube4","FallopianTube4","FallopianTube4","Uterus2","FallopianTube5","Ovary1")
menopause=c("peri","peri","peri","post","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre")
subject=c("Human1","Human1","Human1","Human2","Human3","Human3","Human3","Human4","Human4","Human5","Human5","Human5","Human6","Human6","Human6")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)

subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube","Uterus","FallopianTube2","FallopianTube3","FallopianTube4","Uterus2","FallopianTube5","Ovary1")
ft=grep("Fallopian",organ)
indivft=grep("Fallopian",organs)

all="FallopianTube124"
ft=grep("Fallopian",organ)[c(1:6,9:11)]
indivft=grep("Fallopian",organs)[c(1,2,4)]



### load object for all 5 directly-merged FTs
  load(file=paste0("FallopianTube1-5.Robj"))
dge1=dgeall

### Extract FT1, FT2, FT4 - 3 large healthy FTs
dge=subset(dgeall,Organ %in% c("FallopianTube","FallopianTube2","FallopianTube4"))
table(Idents(dgeall))
table(Idents(dge))
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))

dgefile="plot/Directmerge_"

### Highly variable genes
dge<-FindVariableFeatures(dge)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge), 10)
plot1 <- VariableFeaturePlot(dge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000
### Scale data 
dge <- ScaleData(dge,features=rownames(dge))

### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)
dgeall=dge


save(dgeall,file=paste0(all,".Robj"))

pdf(paste0(dgefile,"HVG_x0.2_y0.2.pdf"),height=6,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot2
dev.off()


########## Analysis for merged data of individual organ

###### Determine top PCs
numPCs=10;i=1
pdf(paste("plot/",all,"_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
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
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,1.2,by=0.1))
dgeall=dge

save(dgeall,file=paste0(all,".Robj"))

print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))
}




###### order clusters for each individual organ
res=c(paste0("RNA_snn_res.0.",c(1)),paste0("integrated_snn_res.0.",c(1,2)));resi=1

## order cell by cluster ID and randomly shuffle cells within each batch
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
### get order of seriation
 do=seriate(da,method="OLO")
levelss=get_order(seriate(da,method="OLO"))

levels=levelss

### Reordered clusters for all cells
cells.use=colnames(dge)
ident=factor(Idents(dge),levels=levels-1)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=names(tmp)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

resorig=dge$RNA_snn_res.0.1
indivftcluster=dge@meta.data$indivftcluster
names(resorig)=names(indivftcluster)=rownames(dge@meta.data)
table(indivftcluster,cells.ident.ordered[names(indivftcluster)])

### re-order based on previous ordering
levels=levelss[c(16:12,9:6,11:10,5:1)]

### Reordered clusters for all cells
cells.use=colnames(dge)
Idents(dge)<-dge[[res[resi]]]
ident=factor(Idents(dge),levels=levels-1)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=names(tmp)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

resorig=dge$RNA_snn_res.0.1
indivftcluster=dge@meta.data$indivftcluster
names(resorig)=names(indivftcluster)=rownames(dge@meta.data)
table(indivftcluster,cells.ident.ordered[names(indivftcluster)])


### save ordered cluster ID in dge object
which(unique(cells.ident.ordered)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
dgeall=dge
save(dgeall,file=paste0(all,".Robj"))
write.table(Idents(dge),file=paste0("plot/",all,"_ident.txt"),row.names=T,col.names=F,quote=F,sep="\t")

### merge into 7 relatively-stable major islands
id=as.numeric(Idents(dge))
names(id)=names(Idents(dge))
id[which(id %in% c(2:6))]<-22
id[which(id %in% c(7:10))]<-33
id[which(id == 11)]<-44
id[which(id %in% c(12:13))]<-55
id[which(id == 14)]<-66
id[which(id == 15)]<-77
id[which(id == 22)]<-2
id[which(id == 33)]<-3
id[which(id == 44)]<-4
id[which(id == 55)]<-5
id[which(id == 66)]<-6
id[which(id == 77)]<-7
table(id)
#    1     2     3     4     5     6     7 
# 2629 14729 18514  4996  5231  1616  3779 
table(Idents(dge),id)
dge$major7<-id
dgeall=dge




### Comparison with clusters of individual segments
ncellsindivclusters=table(Idents(dge),dge@meta.data$indivclusters)
  write.table(ncellsindivclusters,paste0("plot/",organs[i],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")

new11=read.table(file=paste0("plot/",all,"_FT1_19clusters_ident.txt"),row.names=1)
new1=new11[,1]
names(new1)=rownames(new11)
table(new1,Idents(dge))

for(label in c("Organ")){
### Contribution of each batch to each cluster
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@meta.data[,label],Idents(dge))
  print(ncellsbatchcluster)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  print(percentcellsbatch)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  print(percentcellscluster)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0("plot/",numPCs[i],"ncellspercluster_",label,".txt"),quote=F,row.names=T,col.names=T,sep="\t")
}

plotlist=list()
i=1
dge=dgelist[[i]]
Idents(dge) <- dge$RNA_snn_res.0.1ordered
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0.pdf",height=8,width=9)
multiplot(plotlist,cols = 2)
dev.off()


###### Differentially-expressed markers
res=c(paste0("RNA_snn_res.0.",c(1),"ordered"));resi=1

markerslist=list()
	dge=dgelist[[i]]
	Idents(dge)<-dge[[res[i]]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
markerslist[[i]]=markers
write.table(markers,paste0("plot/",all,"_directmerge_mindiff0.2_logfc2fold_9.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")



plotlist=list()
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

pdf(paste0("plot/markerstop.pdf"),height=2*round(size),width=1.8*ceiling(size))
plotlist
dev.off()


pdf(file=paste0("plot/",numPCs[i],"PCs_clusters_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=12)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
dev.off()

