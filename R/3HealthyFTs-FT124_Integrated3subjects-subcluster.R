# 1.17.2021 by Qianyi
# did subclustering for ciliated cells and non-ciliated epithelial cells, respectively
# Related to Figure 2-4, Figure S2-5, and Table S3-4

merge 3 healthy FTs - FT1, FT2, FT4 (9 datasets) together and do clustering
CCA-all genes: selected HVG for each subject before integration; integrated all genes
-> subclustering for each stable islands
clusters 1   ciliated cell type
clusters 2-6 secretory (non-ciliated epithelial) cell type
clusters 7-10 myofibroblast cell type
# notes from Jun: keep integrated data for global integration and do not re-do integration for subsets.
cross-tabulate the subcluster centroids across individual fallopian tube



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
library(RColorBrewer)

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


### using original 17 clusters
subsets=list(2:6,7:10,1:6)
subsetsname=c("2-6","7-10","1-6")
# 2:6  secretory subset   -> 6 subclusters
# 7:10 old stromal subset -> 5 subclusters - very similar as global clusters
# 1:6  ciliated cluster 1 + 6 secretory subclusters

### after renaming the 12 clusters to 1:12
subsets=list(1:8,3:8,1,c(1:2,7),c(1:2,6:8))
subsetsname=c("1-8","3-8",1,"127","12678")
# 1:8   epithelial + new stromal, removing immune cells
# 3:8   new stromal including myofibroblast and endothelial
# 1     ciliated only -> 4 ciliated subclusters
# 127   epithelial   + endothelial cluster 7
# 12678 epithelial + pericyte 6 + endothelial 7&8

subsets=list(c(1,"2_2"))
subsetsname=c("1-2_2")



### CCA-3subjects-allgenes to integrate 3 large healthy FTs (FT1, FT2, FT4)
load(file=paste0(all,"_CCA-3subjects-allgenes-NoC13.Robj"))
dgeall

for(i in 1:length(subsets)){

i=1
i=2
i=3
i=4

cc=subsets[[i]]
ccname=subsetsname[i]
dgefile=paste0("plot/C",ccname,"_")
dge=dgeall
# keep cells in the subset
#dge=subset(dgeall,ft124CCA3all %in% cc) # original 17 clusters
dge=subset(dgeall,ft124CCA3all11 %in% cc)  # renamed 12 clusters
#dge=subset(dgeall,ft124CCA3all22 %in% cc)  # renamed 12 clusters + secretory subclusters
table(Idents(dgeall))
table(Idents(dge))
table(dgeall$Organ)
table(dge$Organ)
# removed non-detected genes
nUMIperGene=apply(dge@assays$RNA@data,1,sum)
genes.use=names(nUMIperGene)[which(nUMIperGene>0)]
genes0=names(nUMIperGene)[which(nUMIperGene==0)]
genes0[which(genes0 %in% VariableFeatures(dge))]
dge <- subset(dge, features = genes.use)
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))
genes0[which(genes0 %in% VariableFeatures(dge))] # integer(0)
length(VariableFeatures(dge)) #[1] 1964 1969
dge
dge1=dge

# split and re-do CCA to determine HVG
DefaultAssay(dge1)<-"RNA"
datalists=SplitObject(dge1,split.by="Organ")
names(datalists)=paste0("FT",c(1,2,4))
datalists
datalist <- lapply(X = datalists, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

### Integration 
k.filter <- min(200, min(sapply(datalist, ncol)))
k.filter # revised k.filter based on the number of cells in the smallest dataset
# 112
anchors <- FindIntegrationAnchors(object.list = datalist,k.filter=k.filter, dims = 1:30) # anchor.features =length(genes.use),
integrated <- IntegrateData(anchorset = anchors, features.to.integrate = genes.use, dims = 1:30)

# change HVG to be from integration of the subset cells
# previous HVG were from global CCA integration
aa=VariableFeatures(integrated)
bb=VariableFeatures(dge)
length(aa) # 2000
length(bb) # 1969 1964 1969 1995
length(aa[which(aa %in% bb)]) # 1029 1137 1425 1683
VariableFeatures(dge)=aa
length(VariableFeatures(dge)) # [1] 2000

### Scale data
# Jun asked me to keep using global CCA-integrated data (v4), not subset CCA-integrated data (saved as integrated above)
DefaultAssay(dge)
dge <- ScaleData(dge,features=rownames(dge))
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)

###### Determine top PCs
numPCs=c(6,6,6)
numPCs=c(9,9,8,10,10)
numPCs=5 # C1-1_2
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
text(eigenvalue+0.2,0.25,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=organs[1],cex=1.5)
dev.off()


###### Using top PCs
### tSNE
  dge <- RunTSNE(dge, dims = 1:numPCs[i])
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,1.2,by=0.1))


Idents(dge) <- dge$ft124CCA3all22

save(dge,file=paste0(all,"_C",ccname,".Robj"))




print(c( length(unique(dge$integrated_snn_res.0.1)),length(unique(dge$integrated_snn_res.0.2)),length(unique(dge$integrated_snn_res.0.3)),length(unique(dge$integrated_snn_res.0.4)),length(unique(dge$integrated_snn_res.0.5)),length(unique(dge$integrated_snn_res.0.6)),length(unique(dge$integrated_snn_res.0.7)),length(unique(dge$integrated_snn_res.0.8)),length(unique(dge$integrated_snn_res.0.9)),length(unique(dge$integrated_snn_res.1)),length(unique(dge$integrated_snn_res.1.1)),length(unique(dge$integrated_snn_res.1.2)) ))


###### order clusters for each individual organ
res=c(paste0("integrated_snn_res.0.",c(1,1,1)));resi=i


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
pdf(file=paste0("plot/Centroid_norm_Seriation_C",ccname,"_",res[resi],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))
 hmap(da) # default method="OLO"
dev.off()

### get order of seriation
 do=seriate(da,method="OLO")
levelss=get_order(seriate(da,method="OLO"))

levels=levelss
if(resi==1){
  levels=rev(levelss)[c(1:4,6,5)]
}
if(resi==2){
  levels=levelss[c(2:4,1,5)]
}
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


### save ordered cluster ID in dge object
which(unique(cells.ident.ordered)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
# save the dge file
dgelist[[resi]]=dge

save(dge,file=paste0(all,"_C",ccname,".Robj"))


### Label clustering IDs
dge$ft124CCA3allC2to6=dge$integrated_snn_res.0.1ordered

dge$ft124CCA3allC7to10=dge$integrated_snn_res.0.1ordered



i=resi
### Comparison with clusters of individual segments
ncellsindivclusters=table(Idents(dge),dge@meta.data$indivclusters)
  write.table(ncellsindivclusters,paste0("plot/",organs[indivft[i]],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")

table(dge$new1,Idents(dge))
table(dge$ft124CCA3all,Idents(dge))


### visualize in PCA, tSNE and/or UMAP for ordered clusters 
# C1 ciliated subclusters
myBrewerPalette=gg_color_hue(4)


# C2-6 secretory subclusters
myBrewerPalette=brewer.pal(7,"Set2")[1:6]


# C1-6 epithelial recluster2-6
if(i==3){
  myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2"))
}

dge2=dge
dge1=dgeCiliated
id2=gsub("-6","",as.character(dge2$ft124CCA3all2))
names(id2)=rownames(dge2@meta.data)
id1=as.character(Idents(dge1))
names(id1)=names(Idents(dge1))
id=c(paste0("1_",id1),id2[which(id2!=1)])
names(id)=c(names(id1),names(id2)[which(id2!=1)])
id=id[names(id2)]
id=factor(id,levels=c("1_1","1_2","1_3","1_4","2_1","2_2","2_3","2_4","2_5","2_6"))
dge$ordered <- id
Idents(dge) <- dge$ordered
save(dge,file=paste0(all,"_C",ccname,".Robj"))


# C7-10 using global clusters but zoomed-in projection: re-did UMAP
dge=dgelist[[2]]
Idents(dge) <- dge$ft124CCA3all
myBrewerPalette=brewer.pal(12,"Paired")[-1]


dge=dgelist[[i]]
Idents(dge) <- dge[[ordered]]
plotlist=list()
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0.pdf",height=8,width=9.5)
multiplot(plotlist,cols = 2)
dev.off()
# remove labels
plotlist=list()
png("plot/clusters_ordered1.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()

plotlist=list()
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=PCAPlot(dge,c(1,4),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=PCAPlot(dge,c(1,5),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[5]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[6]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered1.pdf",height=8,width=14)
multiplot(plotlist,cols = 3)
dev.off()

# C1-8: epithelial + new stromal subset
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[c(2:7,9:12)])

# C3-8: new stromal subset
myBrewerPalette=brewer.pal(12,"Paired")[2:7]

# C1,2,7: epithelial subset+ endothelial cluster 7
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[6])

# C1,2,6,7,8: epithelial subset + new stromal clusters 6 - 8
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[c(5:7)])

# C1-2_2: ciliated C1 recluster + secretory C2_2
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[2])

dge2=dge
dge1=dgeCiliated
id2=as.character(dge2$ft124CCA3all22)
names(id2)=rownames(dge2@meta.data)
id1=as.character(Idents(dge1))
names(id1)=names(Idents(dge1))
id=c(paste0("1_",id1),id2[which(id2!=1)])
names(id)=c(names(id1),names(id2)[which(id2!=1)])
id=id[names(id2)]
id=factor(id,levels=c("1_1","1_2","1_3","1_4","2_2"))
dge$ordered <- id
Idents(dge) <- dge$ordered
save(dge,file=paste0(all,"_C",ccname,".Robj"))



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
dgeall=dge
for(label in c("Part","Organ","Name")){
plotlist=list()
Idents(dgeall)<-dgeall[[label]]
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlim=ylim=list()
pos=c("topright","topright","topleft","bottomright")
# dge4 used top 9 PCs for tSNE
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
  pdf(paste0("plot/",all,"_",label,"_indiv.pdf"),height=2.3*round(sqrt(size)),width=2.3*ceiling(sqrt(size)))
  par(mfrow=c(round(sqrt(size)),ceiling(sqrt(size))),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
} else {
  pdf(paste0("plot/",all,"_",label,"_indiv.pdf"),height=2.5,width=2.5*size)
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
if(size>4 & size<(round(sqrt(size))*ceiling(sqrt(size)))){
  for(new in (size+1):(round(sqrt(size))*ceiling(sqrt(size)))){
    plot.new()
  }
}
}
dev.off()
}
}
}




### Per-cell attributes for each cell type
table(dge$integrated_snn_res.0.1ordered,dge$ft124CCA3allC2to6)
levels=levels(Idents(dge))
for(j in levels){
print(c(mean(dge@meta.data$nFeature_RNA[which(dge$ft124CCA3allC2to6==j)]),mean(dge@meta.data$nCount_RNA[which(dge$ft124CCA3allC2to6==j)]),mean(dge@meta.data$percent.mt[which(dge$ft124CCA3allC2to6==j)])))
}

levels=levels(Idents(dge))
for(j in levels){
print(c(mean(dge@meta.data$nFeature_RNA[which(Idents(dge)==j)]),mean(dge@meta.data$nCount_RNA[which(Idents(dge)==j)]),mean(dge@meta.data$percent.mt[which(Idents(dge)==j)])))
}



pdf(file=paste0("plot/clusters_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=12)
for(i in 1:2){
  plotlist=list()
dge=dgelist[[i]]
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
}
dev.off()
# saved as Figure S2B-C.

### Differentially-expressed markers
numPCs=c(6,6)
numPCs=c(9,9,8)
res=rep("integrated_snn_res.0.1ordered",3)
  Idents(dge)<-dge[[res[i]]]

markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
markerslist[[i]]=markers
write.table(markers,paste0(dgefile,numPCs[i],"PCs_",res[i],"_mindiff0.2_logfc2fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
dge1=dge
DefaultAssay(dge1)="RNA"
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
table(Idents(dge))
table(markers$cluster)
markerslist[[i]]=markers
write.table(markers,paste0(dgefile,numPCs[i],"PCs_",res[i],"_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
dge1=dge
DefaultAssay(dge1)="RNA"
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(1.6))
write.table(markers,paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(markers$cluster)
DefaultAssay(dge1)="RNA"
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(1.6),min.pct=0)
write.table(markers,paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
# saved as Table S3 

plotlist=plotlist2=list()
markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top2
size=sqrt(length(top2$gene))
p <- FeaturePlot(dge, top2$gene, min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
p2 <- FeaturePlot(dge, top2$gene, reduction="tsne",min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
  p2[[j]] <- p2[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p,ncol=5)
plotlist2[[i]]=cowplot::plot_grid(plotlist = p2,ncol=5)
}

pdf(paste0("plot/markerstop.pdf"),height=2*4,width=1.8*5)
plotlist
plotlist2
dev.off()

# secretory subclusters
### Find Markers among 3-4-5 for secretory subclusters
dge1=dge
DefaultAssay(dge1)="RNA"
markers3=FindMarkers(dge1,3,c(4,5),logfc.threshold = log(2),min.diff.pct=0.2)
print(c(length(which(markers3$avg_logFC>0)),length(which(markers3$avg_logFC<0))))
write.table(markers3,paste0(dgefile,"3Vs45_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers3=FindMarkers(dge1,3,c(4,5),logfc.threshold = log(1.6))
print(c(length(which(markers3$avg_logFC>0)),length(which(markers3$avg_logFC<0))))
write.table(markers3,paste0(dgefile,"3Vs45_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers4=FindMarkers(dge1,4,c(3,5),logfc.threshold = log(2),min.diff.pct=0.2)
print(c(length(which(markers4$avg_logFC>0)),length(which(markers4$avg_logFC<0))))
write.table(markers4,paste0(dgefile,"4Vs35_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers4=FindMarkers(dge1,4,c(3,5),logfc.threshold = log(1.6))
print(c(length(which(markers4$avg_logFC>0)),length(which(markers4$avg_logFC<0))))
write.table(markers4,paste0(dgefile,"4Vs35_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers5=FindMarkers(dge1,5,c(3,4),logfc.threshold = log(2),min.diff.pct=0.2)
print(c(length(which(markers5$avg_logFC>0)),length(which(markers5$avg_logFC<0))))
write.table(markers5,paste0(dgefile,"5Vs34_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers5=FindMarkers(dge1,5,c(3,4),logfc.threshold = log(1.6))
print(c(length(which(markers5$avg_logFC>0)),length(which(markers5$avg_logFC<0))))
write.table(markers5,paste0(dgefile,"5Vs34_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
# check the number of markers overlapped among each cluster
markers=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),stringsAsFactors=F)
levels=levels(Idents(dge))
n=length(levels)
cc=matrix(0,n,n)
for(i in 1:n){
    for(j in 1:n){
        cc[i,j]=length(intersect(markers$gene[which(markers$cluster==levels[i])],markers$gene[which(markers$cluster==levels[j])]))
    }
}
cc
length(which(rownames(markers3)[which(markers3$avg_logFC>0)] %in% markers$gene[which(markers$cluster==3)]))
length(which(rownames(markers4)[which(markers4$avg_logFC>0)] %in% markers$gene[which(markers$cluster==4)]))
length(which(rownames(markers5)[which(markers5$avg_logFC>0)] %in% markers$gene[which(markers$cluster==5)]))
markerslist=list(markers3,markers4,markers5)
cc=matrix(0,3,3)
for(i in 1:3){
    for(j in 1:3){
        cc[i,j]=length(intersect(rownames(markerslist[[i]])[which(markerslist[[i]]$avg_logFC>0)],markers$gene[which(markers$cluster==j+2)]))
    }
}
cc
markers=read.table(paste0(dgefile,"RNAassay1_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)
markers1=markers
# check the number of unique (non-overlapped) markers for each cluster
markers=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers3=read.table(paste0(dgefile,"3Vs45_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers4=read.table(paste0(dgefile,"4Vs35_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers5=read.table(paste0(dgefile,"5Vs34_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
table(Idents(dge))
table(markers$cluster)
genelist=list(markers$gene[which(markers$cluster==1)],
  markers$gene[which(markers$cluster==2)],
  unique(c(markers$gene[which(markers$cluster==3)],rownames(markers3)[which(markers3$avg_logFC>0)])),
  unique(c(markers$gene[which(markers$cluster==4)],rownames(markers4)[which(markers4$avg_logFC>0)])),
  unique(c(markers$gene[which(markers$cluster==5)],rownames(markers5)[which(markers5$avg_logFC>0)])),
  markers$gene[which(markers$cluster==6)]
)
names(genelist)=1:6
uniquelist=list()
levels=1:6
n=length(levels)
for(i in 1:n){
  list1=genelist[[i]]
  list2=unique(unlist(genelist[-i]))
  unique=list1[which(!(list1 %in% list2))]
  print(length(unique))
  uniquelist[[i]]=unique
}
cc=matrix(0,n,n)
for(i in 1:n){
    for(j in 1:n){
        cc[i,j]=length(intersect(uniquelist[[i]],uniquelist[[j]]))
    }
}
cc



### combine subclusters 3-4-5 together for secretory subclusters
id=as.numeric(Idents(dge))
names(id)=names(Idents(dge))
id[which(id %in% c(3:5))]<-"3-5"
table(id)
id=factor(id,levels=c(1,2,"3-5",6)) 
table(Idents(dge),id)
table(id)
    1     2   3-5     6 
  268   868 12550   316 

dge$ft124CCA3allC2to6merge345<-id
Idents(dge) <- dge$ft124CCA3allC2to6merge345

dge1=dge
DefaultAssay(dge1)="RNA"
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
table(Idents(dge))
table(markers$cluster)
write.table(markers,paste0(dgefile,"merge345_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
DefaultAssay(dge1)="RNA"
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(1.6))
table(Idents(dge))
table(markers$cluster)
write.table(markers,paste0(dgefile,"merge345_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")


FindMarkers(dge1,"3-5",c(1,2,6),only.pos=TRUE,logfc.threshold = log(1.6),min.pct=0)
FindMarkers(dge1,"3-5",c(1,2,6),only.pos=TRUE,logfc.threshold = log(1.4),min.pct=0)
# 1 marker


######## Heatmap for all markers
avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0(dgefile,"Centroid-uncorrected.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0(dgefile,"Centroid-postCCAallgenes.txt"),row.names=T,col.names=T,quote=F,sep="\t")

centroid=read.table(paste0(dgefile,"Centroid-postCCAallgenes.txt"))
colnames(centroid)=gsub("X","",colnames(centroid))
centroid[1:5,1:3]

### Genes Standardized Across Cell Types
# note: used this
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)

### load markers for each cluster
# 2FC, 20%diff, 10%pct
markers=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
table(Idents(dge))
table(markers$cluster)
genes=markers$gene

# 1.6FC, 10%pct
markers=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
table(Idents(dge))
table(markers$cluster)
genes=markers$gene

# 1.6FC only
markers=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
table(Idents(dge))
table(markers$cluster)
genes=markers$gene

# Union of Global 1VsAllOthers: 1.6FC, 10%pct + Local 3,4,5VsNeighbor for secretory subclusters
markers=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers3=read.table(paste0(dgefile,"3Vs45_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers4=read.table(paste0(dgefile,"4Vs35_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers5=read.table(paste0(dgefile,"5Vs34_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
genelist=list(markers$gene[which(markers$cluster==1)],
  markers$gene[which(markers$cluster==2)],
  unique(c(markers$gene[which(markers$cluster==3)],rownames(markers3)[which(markers3$avg_logFC>0)])),
  unique(c(markers$gene[which(markers$cluster==4)],rownames(markers4)[which(markers4$avg_logFC>0)])),
  unique(c(markers$gene[which(markers$cluster==5)],rownames(markers5)[which(markers5$avg_logFC>0)])),
  markers$gene[which(markers$cluster==6)]
)
genes=unlist(genelist)


### Visualize markers in heatmap across all cell types
data.use=centroid.std

levels=colnames(centroid.std)

colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
col.lab=rep("",length(levels))
col.lab=gsub(".*_","",levels)

ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cell Type")
colnames(clab)=c("Cell Type","")

col.use=redblue100

data.use=centroid.std[markers$gene,]
row.lab=rownames(data.use)

write.table(data.use,paste0(dgefile,"1_Centroid-postCCAallgenes_markers_std.txt"),row.names=T,col.names=T,quote=F,sep="\t")


jpeg(file=paste0(dgefile,"centroid_std_markersall.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0(dgefile,"centroid_std_markersall2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
# saved as Figure S2D.

### cell-cell correlation for C7-10 recluster
# check if subcluster 4 cells have higher correlation with two other groups
data=as.matrix(dge@assays$integrated@data[VariableFeatures(dge),])
rho=cor(data,method="sp") # rank correlation using HVG
save(rho,file=paste0(all,"_C",ccname,"_allcells_rho.Robj"))

id=sort(Idents(dge))
cells=names(id)
rho=rho[cells,cells]
summary(c(rho))
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.005621 0.367369 0.412608 0.418969 0.464273 1.000000
png(file=paste0(dgefile,"allcells_rho.png"),res=300,height=5000,width=5000)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
image(rho,col=redblue100,zlim=c(-1,1))
dev.off()
png(file=paste0(dgefile,"allcells_rho0-1.png"),res=300,height=5000,width=5000)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
image(rho,col=redblue100,zlim=c(0,1))
dev.off()
png(file=paste0(dgefile,"allcells_id.png"),res=300,height=800,width=3200)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
image(cbind(id,id))
dev.off()


# visualize known markers:
knownmarkers=c("RUNX3","ACTA2","EPCAM","JCHAIN","KIT","CD3E","CAVIN2","PTPRC")
length(knownmarkers) #8
gene=knownmarkers[which(knownmarkers %in% rownames(dge))]
length(gene) #8
### violin plot
plotlist=list()
plotlist2=list()
for(j in 1:2){
  dge<-dgelist[[j]]
  plotlist[[j]]=VlnPlot(dge,gene,ncol=4,pt.size=-1,cols=myBrewerPalette)
  DefaultAssay(dge) <- "RNA"
  plotlist2[[j]]=VlnPlot(dge,gene,ncol=4,pt.size=-1,cols=myBrewerPalette)
}
pdf(paste0(dgefile,"ordered0_knownmarkers_Violin_integrated.pdf"),height=4,width=10)
print(plotlist)
print(plotlist2)
dev.off()

# 2.2.2021
# 4) visualize known markers for the subclustering of nonciliated epithelial cells
a.  Cluster 1: CAVIN2, VWF, SOX17, DCN, PECAM1, MCAM
b.  Cluster 2: ACTA2, TAGLN, LSGALS1, DES, MYLK, SPARC, S100A4, PDGFRB, PDGFRA, THY1
c.  Cluster 3:
d.  Cluster 4: LAMB3, PLAUR, CXCL8, CCL20, BRIC3, LIF, MSLN
e.  Cluster 5: PGR, CYP1B1, COL1A2, OVGP1, LGR5, POSTN, CD44
f.  Cluster 6: RUNX3, CD3E, TPSAB1, Kit
g.  EMT specific: VIM, TIMP3, SPARC, COL1A1, TWIST1, TWIST2, S100A4, PRRX1, MMP2, ZEB1, ZEB2, FN1, CDH11, CDH2, SNAI1, SNAI2, MUC1, OCLN, LOX, CD44, PROM1 (CD133), POU5F1 (OCT-4), SUSD2
h.  Hormone receptors: ESR, PGR, AR

knownmarkers=c("CAVIN2","VWF","SOX17","DCN","PECAM1","MCAM",
"ACTA2","TAGLN","LGALS1","DES","MYLK","SPARC","S100A4","PDGFRB","PDGFRA","THY1",
"LAMB3","PLAUR","CXCL8","CCL20","BIRC3","LIF","MSLN",
"PGR","CYP1B1","COL1A2","OVGP1","LGR5","POSTN","CD44",
"RUNX3","CD3E","TPSAB1","KIT",
"VIM","TIMP3","SPARC","COL1A1","TWIST1","TWIST2","S100A4","PRRX1","MMP2","ZEB1","ZEB2","FN1","CDH11","CDH2","SNAI1","SNAI2","MUC1","OCLN","LOX","CD44","PROM1","POU5F1","SUSD2",
"ESR1","PGR","AR"
)
length(knownmarkers) # 60
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
# "LSGALS1" "BRIC3"   "ESR"
# LSGALS1 is typo for LGALS1
# BRIC3 is typo for BIRC3
# ESR should be ESR1
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=10,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers1_Violin_RNA.pdf"),height=12,width=25)
print(plotlist)
dev.off()

# 4.13.2021
additional violin plots for the nonciliated epithelial cell subclusters for the FT124 analysis.

knownmarkers=c("ITGA6","MUC16","PTPRC","CD8A","CD4")
length(knownmarkers) # 60
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
# "CD45"
# CD45 should be PTPRC
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=5,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers2_Violin_RNA.pdf"),height=2,width=12.5)
print(plotlist)
dev.off()


4.15.2021
knownmarkers=c("CRISP3","IGHA1","IER2","FOS","ATF3","BHLHE41",
  "KIT","CD34","ACTA2","DES","VIM","CAVIN2","CAVIN1","B3GAT1","NES","CD1A","SELP","CD68"
)
# B3GAT1 (CD57/NK1), NES (nestin), CD1A, SELP (CD62P)

length(knownmarkers) # 18
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
# "CD45"
# CD45 should be PTPRC
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=6,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers3_Violin_RNA.pdf"),height=6,width=15)
print(plotlist)
dev.off()

4.26.2021
# 3)  Non-ciliated secretory population
a.  Violin plots
i.  Figure 2D: 
1.  Comparison of 2-2 and 2-6: Pax 8, KRT7, EPCAM
2.  Cluster 2-5 “stem markers” CD24, TNFRSF19 (TROY)
ii. Figure 2E:
1.  EMT (From Hu Cancer cell EM cluster): A2M, TPM2, MFAP4, CRISPLD2, MYH11, SPARCL1, SLC2A3
2.  WNT: RSPO1- regulator of wnt signalling (Hu, Cancer Cell)

knownmarkers=c("PAX8","KRT7","EPCAM",
"CD24","TNFRSF19",
"RSPO1"
)

length(knownmarkers) # 6
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=6,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers4_Violin_RNA.pdf"),height=2,width=15)
print(plotlist)
dev.off()


# 6)  Ciliated C1 4 subclusters
knownmarkers=c("FOXO4","NR2F2","CAPS","FOXJ1","EPCAM","KRT7","ACTA2","PAX8","CD44","PDGFRB","PDGFRA",
  "KIT","SPARC","TIMP3","TIMP2","PRRX1","MFAP4","CRISPLD2","LMOD1","FBN1","ADAMTS4","MEG3","PMP22","MMP2","SLC2A3","SPARCL1","COL6A3","AEBP1","OGN","CAVIN2","PLPP3","CALD1","C1S","COL6A2","COL3A1","EMP3","FBLN1","COL4A2","LGALS1","SERPINF1","NEXN","LUM","ADAMTS1","DCN","DPT","AKAP12","SFRP4","COL4A1","SFRP1","SELENOM","BST2","FILIP1L","A2M","ACKR3","GEM","TGFB1I1","C11orf96","SYNPO2","TPM2","ADH1B","AXL","CAVIN1","PTGDS","FHL1","TFPI","IGFBP6","PPP1R14A","C1R","MGP","PHLDA1","GPRC5A","FGF7","GYPC","TUBB6","COL1A1","IGF1","CCDC80","IGFBP5","IGFBP4","ATP1B3","VIM","C7","IL6","BAG3","PTMS","TAGLN","MCAM","GNG11","MYC","KLF9","THY1","ACTA2","F3","MYL9","RGS2","NBL1","MT1M","MFGE8","MYH11","CNN1","RAMP1","FLNA","UBB","EMP1","UAP1","MYADM","HSPA6","PPP1R12A","MT1A","TMSB4X","ID3","ADIRF","CAV1","DNAJB1","FTL","ZNF331","DES","HSPA1A","TPM1","ACTG2","MYLK","CREM","FOSB","LDHA","S100A4","RGS16","FSTL1","KLF2","ACTB","KLF4","COL1A2","SRGN","ID2","CSRP1","LPP","APOD","GSN","LMNA","TCEAL4","MT1X","MT2A","SELENOP","ACKR1","TIMP1","ICAM1","MDK","CCL21","CCL14","CCL2","CFD",
  "S100A8","CCNO","NME2","S100A9","EVL","AHI1","LGALS1","CDC20B",
  "TTC29","ANXA13",
  "MSLN","ASS1","TTYH1","FAM107A","LYPD1","S100A16","HMGA1","CXCL3","CRISP3","TACSTD2","CXCL1","CLDN1","TUBB","SERPINA3","CCL20","ID4",
  "ZEB1","ZEB2",
  "ESR1","PGR","AR",
  "CYP1B1","CYP4B1","CYP51A1"
)

length(knownmarkers) # 184
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=8,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers_Violin_RNA.pdf"),height=46,width=20)
print(plotlist)
dev.off()
 
# 5.27.2021 violins for these genes only in the ciliated subtypes
knownmarkers=c("RSPH1","RSPH9","RSPH4A","DNAI1","DNAI2","DNAH11","CCDC103","LRRC6","ZMYND10","FOXJ1"
  ,"MUC1")
length(knownmarkers) # 11
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=4,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers2_Violin_RNA.pdf"),height=6,width=10)
print(plotlist)
dev.off()

# 11.18.2021 violin plots for two markers in ciliated subtypes
knownmarkers=c("CDKN2A","CDKN2B")
length(knownmarkers) # 2
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=2,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkersCDKN2A_Violin_RNA.pdf"),height=2,width=5)
print(plotlist)
dev.off()


# 10.26.2021 violin pots for secretory subtypes
knownmarkers=c("ENG")
length(knownmarkers) #8
gene=knownmarkers[which(knownmarkers %in% rownames(dge))]
length(gene) #8
### violin plot
plotlist=list()
plotlist2=list()
j=1
  plotlist[[j]]=VlnPlot(dge,gene,ncol=1,pt.size=-1,cols=myBrewerPalette)
  DefaultAssay(dge) <- "RNA"
  plotlist2[[j]]=VlnPlot(dge,gene,ncol=1,pt.size=-1,cols=myBrewerPalette)
}
pdf(paste0(dgefile,"ordered0_ENG_Violin_integrated.pdf"),height=2,width=2.5)
print(plotlist)
print(plotlist2)
dev.off()


# 2.14.21 notes
#b.  We will also need a dot plot for the subclustering of epithelial cells. We can start with this order. 
i.  Cluster 1: CAVIN2, VWF, SOX17, PECAM1, MCAM
ii. Cluster 2: DCN, ACTA2, TAGLN, LSGALS1, DES, MYLK, S100A4
iii.  Cluster 3:
iv. Cluster 4: LAMB3, PLAUR, CXCL8, CCL20, BIRC3, LIF, MSLN
v.  Cluster 5: PGR, ESR, AR, CYP1B1, COL1A2, OVGP1, LGR5, POSTN, CD44
vi. Cluster 6: RUNX3, CD3E

knownmarkers=c("CAVIN2","VWF","SOX17","PECAM1","MCAM",
  "DCN","ACTA2","TAGLN","LGALS1","DES","MYLK","S100A4",
  "LAMB3","PLAUR","CXCL8","CCL20","BIRC3","LIF","MSLN",
  "PGR","ESR1","AR","CYP1B1","COL1A2","OVGP1","LGR5","POSTN","CD44",
  "RUNX3","CD3E"
)
length(knownmarkers) # 60
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
# "LSGALS1"   "ESR"
# LSGALS1 is typo for LGALS1
# ESR should be ESR1
gene=knownmarkers

DefaultAssay(dge) <- "integrated"
pdf(paste0(dgefile,"ordered1_knownmarkers_DotPlot_RNA.pdf"),height=4,width=10)
DotPlot(dge, features = gene) + RotatedAxis()
DotPlot(dge, features = gene,cols=c("white","blue")) + RotatedAxis() # use this
DotPlot(dge, features = gene,cols=c("white","red")) + RotatedAxis()
dev.off()

# 4.10.2021
add EPCAM and PAX8 

gene=c("EPCAM","PAX8",knownmarkers)

DefaultAssay(dge) <- "integrated"
pdf(paste0(dgefile,"ordered1_knownmarkers_DotPlot_RNA.pdf"),height=4,width=10)
DotPlot(dge, features = gene) + RotatedAxis()
DotPlot(dge, features = gene,cols=c("white","blue")) + RotatedAxis() # use this
DotPlot(dge, features = gene,cols=c("white","red")) + RotatedAxis()
dev.off()





5.14.2021
# C1-2 epithelial subset: ciliated 4 subclusters + secretory 6 subclusters subset
markers across 3 segments

ccname="1-6"
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6])
DefaultAssay(dge) <- "RNA"      
dge$Part=factor(dge$Part,levels=c("Fimbria","Ampulla","Isthmus"))  

knownmarkers=c("ALDH1A1","ALDH1A2")

### violin plot
DefaultAssay(dge) <- "RNA"
plotlist=list()
for(gene in knownmarkers){
  plotlist[[gene]]=VlnPlot(dge,gene,split.by="Part",pt.size=-1,cols=gg_color_hue(3))
}
pdf(paste0(dgefile,"segment_Violin_RNA.pdf"),height=2.5,width=7)
print(plotlist)
dev.off()

DefaultAssay(dge) <- "RNA"
  for(subject in unique(dge$Organ)){
    dge1=subset(dge,Organ == subject)
plotlist=list()
for(gene in knownmarkers){
  plotlist=VlnPlot(dge1,gene,split.by="Part",pt.size=-1,cols=gg_color_hue(3))
png(paste0(dgefile,subject,gene,"_segment_Violin_RNA.png"),res=300,,height=1000,width=3600)
print(plotlist)
dev.off()
}
}
# saved as Figure S4.

# 3.27.21 notes
#10) Dot plot for new stromal cells (clusters 3-8) – Figure 3
a.  Cluster 3 (former cluster 7): fibroblast: DCN, CD34, LGALS1, COL1A1, thy1, ADH1B, FN1, CYP1B1, POSTN
b.  Cluster 4, myofibroblast (former cluster 8): PRRX1, PGR, PDGFRA, POSTN, DES, MYLK, ACTA2, SLMAP, 
c.  Cluster 5, Smooth muscle (Former cluster 9): DES, ACTG2, TAGLN, MYH11, SLMAP, 
d.  Cluster 6, Pericyte (Former cluster 10): MCAM, PDGFRB, KIT, CSPG$
e.  Cluster7, Blood endothelia (Former cluster 11): SELE, VWF, CD34, VCAM1, CD36, ESAM, NRP1, SOX17, CAVIN2
f.  Cluster8, Lymphatic endothelial (Former cluster 12): PROX1, PDPN, LYVE1, NRP2, MMRN1

knownmarkers=unique(c("DCN","CD34","LGALS1","COL1A1","THY1","ADH1B","FN1","CYP1B1","POSTN",
  "PRRX1","PGR","PDGFRA","POSTN","DES","MYLK","ACTA2","SLMAP",
  "DES","ACTG2","TAGLN","MYH11","SLMAP",
  "MCAM","PDGFRB","KIT","CSPG4",
  "SELE","VWF","CD34","VCAM1","CD36","ESAM","NRP1","SOX17","CAVIN2",
  "PROX1","PDPN","LYVE1","NRP2","MMRN1"
))
length(knownmarkers) # 36
knownmarkers1=knownmarkers
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

DefaultAssay(dge) <- "integrated"
pdf(paste0(dgefile,"ordered1_knownmarkers_DotPlot_RNA.pdf"),height=4,width=10)
DotPlot(dge, features = gene) + RotatedAxis()
DotPlot(dge, features = gene,cols=c("white","blue")) + RotatedAxis() # use this
DotPlot(dge, features = gene,cols=c("white","red")) + RotatedAxis()
dev.off()

#5.2 fix stromal dot plot
a.  Clusters/markers – Move PDGFRA, CD36, CSPG4, ACTG2, CD36; Delete NRP1
b.  New Order: DCN, CD34, LGALS1, COL1A1, THY1, ADH1B, FN1, CYP1B1, PDGFRA, POSTN, 
PRRX1, DES, MYLK, ACTG2, ACTA2, SLMAP, MYH11, TAGLN, 
MCAM, PDGFRB, CSPG4, KIT, 
SELE, VWF, VCAM1, SCAM, SOX17, CD36, CAVIN2, 
PROX1, PDPN, LYVE1, NRP2, MMRN1
# SCAM is typo for ESAM

knownmarkers=unique(c("DCN","CD34","LGALS1","COL1A1","THY1","ADH1B","FN1","CYP1B1","PDGFRA","POSTN","PRRX1","DES","MYLK","ACTG2","ACTA2","SLMAP","MYH11","TAGLN","MCAM","PDGFRB","CSPG4","KIT","SELE","VWF","VCAM1","ESAM","SOX17","CD36","CAVIN2","PROX1","PDPN","LYVE1","NRP2","MMRN1"))
length(knownmarkers) # 34
knownmarkers1[which(!(knownmarkers1 %in% knownmarkers))]
# [1] "PGR"  "NRP1"
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers


DefaultAssay(dge) <- "integrated"
pdf(paste0(dgefile,"ordered1_knownmarkers2_DotPlot_RNA.pdf"),height=4,width=10)
DotPlot(dge, features = gene) + RotatedAxis()
DotPlot(dge, features = gene,cols=c("white","blue")) + RotatedAxis() # use this
DotPlot(dge, features = gene,cols=c("white","red")) + RotatedAxis()
dev.off()



# 8)  Stromal (cluster 3-8) violin plot for just healthy FT124 for genes of interest
knownmarkers=unique(c("TP53","CDKN1A","CDKN2A","COL1A2","ADH1B","COL3A1","THY1","S100A10","LUM","COL6A1","GAS1","PDGFRL","RGCC","CD34","RGS2","CHRDL1","FN1","CITED2","COL15A1","CYP1B1","PLAGL1","COL6A3","MSX1","PIK3R1","HAS1","H19","COL5A2","APOE","TWIST2","IFITM1","SFRP4","ECM1","ALDH1A2","POSTN","DES","ACTG2","MYLK","MYH11","ACTA2","SLMAP","CLDN1","AREG","NOTCH3","CDKN1A","CRIM1","LBH","ALKAL2","PDGFA","ZFHX3","ITGA8","RRAS","RASGRP2","SIVA1","DDIT4"))
#P53: TP53
#P21: CDKN1A
#P16: CDKN2A
length(knownmarkers) # 53
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=10,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers1_Violin_RNA.pdf"),height=12,width=25)
print(plotlist)
dev.off()

4.16.2021
knownmarkers=c(
  "KIT","CD34","ACTA2","DES","VIM","CAVIN2","CAVIN1","B3GAT1","NES","CD1A","SELP","CD68"
)
# B3GAT1 (CD57/NK1), NES (nestin), CD1A, SELP (CD62P)

length(knownmarkers) # 12
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
# "CD1A"
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=6,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers2_Violin_RNA.pdf"),height=4,width=15)
print(plotlist)
dev.off()


# 5.3.2021
knownmarkers=unique(c("KIT","ESR1","PGR","AR","THY1","ADH1B","PRRX1","SOX17","CD44","SUSD2","RUNX3","S100A9",
"KRT7","EPCAM",
"CAPS","FOXJ1",
"PAX8","OVGP1",
"PDGFRA","POSTN","COL1A1",
"DCN","VIM","NR2F2",
"ACTA2","DES","MYH11",
"MCAM","PDGFRB","CSPG4","TAGLN",
"PECAM1","KDR","VWF","PROX1","PDPN",
"KIT","TPSAB1","TPSB2","MS4A2"))

length(knownmarkers) # 39
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=7,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered1_knownmarkers3_Violin_RNA.pdf"),height=12,width=17.5)
print(plotlist)
dev.off()

### adjust Violin plot to be wider for those only showing ticks
cols=gg_color_hue(4) 
data=data.frame(t(dge@assays$RNA@data[knownmarkers,]),CellType=dge$integrated_snn_res.0.1ordered)
for(GENE in knownmarkers){
p1=ggplot(data,aes(x=CellType,y=get(GENE)))+geom_violin(aes(fill=CellType),scale="width",adjust=2)+scale_fill_manual(values = cols)+ theme(legend.position = "none")+
theme_bw()+
theme(axis.text.x = element_text(angle = 45,hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_blank(),
       axis.line    = element_line(color='black'))+
ggtitle(GENE)+
ylab("Expression Level")
jpeg(paste0(dgefile,"ordered1_knownmarkers3_Violin_",GENE,"v2.jpg"),res=300,height=660,width=962.5)
print(p1)
dev.off()
}


5.27.2021
knownmarkers=c(
  "MSLN","MUC1","CYP1B1","OVGP1","S100A4","POU5F1")

length(knownmarkers) # 6
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

### violin plot
DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=3,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers4_Violin_RNA.pdf"),height=4,width=7.5)
print(plotlist)
dev.off()

# saved as Figure 2B-E, 3D-G, S3F

# 11.9.2020
### Visualize known EMT markers - received from Nicole on 10/28/2020
### 18 EMT markers
# SANI2: typo, should be SNAI2
# PPRX1: typo, should be PRRX1
# ER: ESR1
# PR: PGR
EMT=c("ZEB2","TWIST1","SNAI2","SNAI1","S100A4","FN1","GSK3B","MMP2","SUMO2","CDH1","MUC1","PRRX1",
"CD44","ACTA2","RUNX3","ESR1","PGR","EPCAM"  )
EMT[which(!(EMT %in% rownames(dge@assays$RNA@data)))]
length(EMT)  # 18
markerslist=list(c("CD44","ACTA2"),
  c("CD44","RUNX3"),
  c("CD44","ESR1"),
  c("CD44","PGR"),
  c("ACTA2","ESR1"),
  c("ACTA2","PGR"),
  c("RUNX3","ESR1"),
  c("RUNX3","PGR"),
  c("EPCAM","CD44"),
  c("EPCAM","ACTA2"),
  c("EPCAM","RUNX3"),
  c("CD45","RUNX3"),
c("MKI67","CAPS")
 )

knownmarkers=EMT
i=1
plotlist=list()
dge=dgelist[[1]]
gene=knownmarkers[which(knownmarkers %in% rownames(dge@assays$RNA@data))]
length(gene) # 18



### fraction of cells expression the marker in each cluster
library(ineq)
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
tmp=list()
  for(g in gene){
print(g)
tmp[[g]]=list()
for(id in levels(Idents(dge))){
tmp[[g]][[id]]=GetAssayData(dge)[g,which(Idents(dge)==id)]
}
}
    tt=pct=matrix(0,ncol=length(levels(Idents(dge))),nrow=length(gene))
    rownames(tt)=rownames(pct)=gene
    colnames(tt)=colnames(pct)=levels(Idents(dge))
  for(g in gene){
for(id in levels(Idents(dge))){
# Number of Non-0 cells
tt[g,id]=sum(tmp[[g]][[id]]!=0)
# % of Non-0 cells
pct[g,id]=sum(tmp[[g]][[id]]!=0)/length(which(Idents(dge)==id))
}
}
print(tt)
print(pct)
}

### fraction of cells co-expressing 2 markers in each cluster
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
tmp=list()
  for(g in 1:length(markerslist)){
print(g)
genes=markerslist[[g]]
tmp[[g]]=list()
for(id in levels(Idents(dge))){
tmp[[g]][[id]]=GetAssayData(dge)[genes,which(Idents(dge)==id)]
}
}
    tt=pct=matrix(0,ncol=length(levels(Idents(dge))),nrow=length(markerslist))
    colnames(tt)=colnames(pct)=levels(Idents(dge))
    tt0=pct0=tt
  for(g in 1:length(markerslist)){
genes=markerslist[[g]]
for(id in levels(Idents(dge))){
  tmpgenes=tmp[[g]][[id]]
# Number of Non-0 cells
tt[g,id]=length(which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0))
# % of Non-0 cells
pct[g,id]=length(which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0))/length(which(Idents(dge)==id))
# Number of Double-0 cells
tt0[g,id]=length(which(tmpgenes[1,]==0 & tmpgenes[2,]==0))
# % of Double-0 cells
pct0[g,id]=length(which(tmpgenes[1,]==0 & tmpgenes[2,]==0))/length(which(Idents(dge)==id))
}
}
print(tt)
print(pct)
print(tt0)
print(pct0)
}
}

DefaultAssay(dge) <- "RNA"
taball=NULL
  for(g in 1:length(markerslist)){
genes=markerslist[[g]]
tmpgenes=GetAssayData(dge)[genes,]
tab=matrix(c(length(which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0)),
length(which(tmpgenes[1,]!=0 & tmpgenes[2,]==0)),
length(which(tmpgenes[1,]==0 & tmpgenes[2,]!=0)),
length(which(tmpgenes[1,]==0 & tmpgenes[2,]==0)) ), nrow=2)
fisher.test(tab)
od=fisher.test(tab)$estimate
p=fisher.test(tab)$p
gene1=tmpgenes[1,which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0)]
gene2=tmpgenes[2,which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0)]
rho=cor(gene1,gene2,method="sp")
taball=rbind(taball,c(tab,od,p,rho))
}
taball


# average expression
print(c(expMean(tmp1),expMean(tmp2),expMean(tmp3),expMean(tmp4)))
# average expression in Non0 Cells
print(c(expMean(tmp1[tmp1!=0]),expMean(tmp2[tmp2!=0]),expMean(tmp3[tmp3!=0]),expMean(tmp4[tmp4!=0])))
# Gini in all cells
print(c(ineq(tmp1),ineq(tmp2),ineq(tmp3),ineq(tmp4)))
}


### Feature plot 
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
  for(genes in markerslist){
jpeg(paste0(dgefile,paste(genes,collapse="_"),"_",names(datalist)[j],".jpeg"),res=300,height=1000,width=3800)
  p<-FeaturePlot(dge,genes,pt.size=.5, blend=TRUE,blend.threshold=0.25)
  print(p)
dev.off()
}
}

### FeatureScatter
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
  for(genes in markerslist){
jpeg(paste0(dgefile,paste(genes,collapse="_"),"_Scatter_",names(datalist)[j],".jpeg"),res=300,height=1000,width=1200)
 p<-FeatureScatter(dge,genes[1],genes[2],pt.size=.5,cols=myBrewerPalette)
  print(p)
dev.off()
}
}
# saved as Figure S3C