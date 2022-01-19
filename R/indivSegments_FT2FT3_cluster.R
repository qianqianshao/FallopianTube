# 8.24.2020 by Qianyi
# individual clustering for each segment of FT2 and FT3

FT2: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT3: 2 segments of fallopian tube from 1 pre-menopausal women with disease, sequenced in 1 batch
(1) Fimbria and Ampulla 
(2) Isthmus
individual clustering for each part

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




########## Analysis for individual dataset
# Load the raw dataset from CellRanger
runbatch=5:n
rr=length(runbatch)
datalist=list()
for(i in runbatch){
  dir=paste0("/nfs/turbo/umms-hammou/10x/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/filtered_feature_bc_matrix/")
  data <- Read10X(data.dir = dir)
  datalist[[i]]=data
}
setlist=dgelist=list()
countcell=countgene=NULL
for(i in runbatch){
  dgedata=as.matrix(datalist[[i]])
  ### Filter for cells (>500 genes and <10% MT)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^MT-|^mt-", rownames(dgedata.tmp), value = T) 
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1] 
print(c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
countcell=rbind(countcell,c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
### Filter for genes (Non-0 genes)
nCellperGene <- rowSums(dgedata.tmp>0)
dgedata2=dgedata.tmp[which(nCellperGene >= 1),]
print(c(nrow(dgedata),nrow(dgedata2)))
countgene=rbind(countgene,c(nrow(dgedata),nrow(dgedata2)))
print(summary(rowSums(dgedata2)))
setlist[[i]]=dgedata2
  dgedata2=setlist[[i]]
  dge <- CreateSeuratObject(counts = dgedata2, project = dataset[i], min.cells = 0, min.features = 0)
  dge[["percent.mt"]] <- PercentageFeatureSet(dge, pattern = "^MT-")
  dgelist[[i]]=dge
  save(dge,file=paste0(dgefile,dataset[i],".Robj"))
}
for(i in runbatch){
  dgedata2=setlist[[i]]
write.table(dgedata2,file=paste0(dataset[i],"rawexp.txt"),col.names=T,row.names=T,sep="\t",quote=F)
}

### Per-cell attributes: Depth, Average number of genes, UMIs, %MT per cell
for(i in runbatch){
dge=dgelist[[i]]
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))
}


# Visualize QC metrics as a violin plot
plot11=plot22=plot33=plot1=plot2=list()
for(i in runbatch){
dge=dgelist[[i]]
ii=i-min(runbatch)+1
plot11[[ii]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1)+ theme(legend.position = 'none')
plot22[[ii]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+ theme(legend.position = 'none')
plot33[[ii]]=VlnPlot(dge, features = "percent.mt",pt.size=-1)+ theme(legend.position = 'none')
plot1[[ii]] <- FeatureScatter(dge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2[[ii]] <- FeatureScatter(dge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
}
pdf(paste0("plot/FeatureScatter.pdf"),height=4.5,width=rr*5)
multiplot(plot1,cols=rr)
multiplot(plot2,cols=rr)
dev.off()
pdf(paste0("plot/VlnPlot.pdf"),height=4.5,width=rr*3)
multiplot(plot11,cols=rr)
multiplot(plot22,cols=rr)
multiplot(plot33,cols=rr)
dev.off()



runbatch=5:n
rr=length(runbatch)

plot2=list()
for(i in runbatch){
dge=dgelist[[i]]
ii=i-min(runbatch)+1

### Add Batch description
  dge[["Name"]] <- names[i]
  dge[["Part"]] <- part[i]
  dge[["Organ"]] <- organ[i]
  dge[["Run"]] <- run[i]
  dge[["Menopause"]] <- menopause[i]
  dge[["Subject"]] <- subject[i]
### Normalize data
dge<-NormalizeData(dge)
### Highly variable genes
dge<-FindVariableFeatures(dge)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge), 10)
plot1 <- VariableFeaturePlot(dge)
plot2[[ii]] <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000
### Scale data 
dge <- ScaleData(dge,features=rownames(dge))
  print(i)
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)
dgelist[[i]]=dge
save(dge,file=paste0(dataset[i],".Robj"))
}


pdf(paste0("plot/HVG_x0.2_y0.2.pdf"),height=6,width=16)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(plot2,cols=3)
dev.off()


###### Determine top PCs
numPCs= c(6,8,8,8,7,8,8,7,8)

pdf(paste("plot/dge_PCA_Variablel_variation.pdf",sep=""),height=6,width=3*rr)
par(mfrow=c(2,rr),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(i in runbatch){
dge=dgelist[[i]]
plot(Stdev(dge,reduction="pca")[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
legend("topright",legend=dataset[i],cex=1.5)
abline(v=numPCs[i]+0.5,lty=2)
}
### density plot of Eigenvalue
for(i in runbatch){
dge=dgelist[[i]]
eigenvalue=Stdev(dge,reduction="pca")[numPCs[i]]
print(eigenvalue)
plot(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(Stdev(dge,reduction="pca")),col="black")
lines(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.25,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=dataset[i],cex=1.5)
}
dev.off()


###### Using top PCs
for(i in runbatch){
dge=dgelist[[i]]
### tSNE
	dge <- RunTSNE(dge, dims = 1:numPCs[i])
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,3,by=0.1))
dgelist[[i]]=dge
save(dge,file=paste0(dataset[i],".Robj"))
}


for(i in runbatch){
dge=dgelist[[i]]
print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))
}



###### order clusters for each dataset
res=paste0("RNA_snn_res.0.",c(6,4,3,4,1,1,1,1,1))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
for(resi in runbatch){
dge=dgelist[[resi]]
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
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]
if(resi==2){
  levels=c(levelss[[resi]][1:3],rev(levelss[[resi]][4:7]),levelss[[resi]][8:11])
}
if(resi==5 | resi==7){
  levels=rev(levelss[[resi]])
}

reorder=list(
c(1:5,9:6,10:12),
c(1:3,6:4,7:10),
c(1:4,8:5,9:11),
c(1:3,6:4,10:7),
c(1:2,4:6,3,8,7,9,10)
)
if(resi>=5){
  levels=levels[reorder[[resi-4]]]
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
dgelist[[resi]]=dge
}

### save dge with ordered clusters
for(i in runbatch){
dge=dgelist[[i]]
save(dge,file=paste0(dataset[i],".Robj"))
}


### PCA and tSNE for ordered clusters of each individual replicate
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist3[[i]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist4[[i]]=PCAPlot(dge,c(1,4),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist5[[i]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0.pdf",height=4*2,width=4.5*rr)
multiplot(plotlist2,cols = rr)
multiplot(plotlist3,cols = rr)
multiplot(plotlist4,cols = rr)
multiplot(plotlist5,cols = rr)
multiplot(plotlistt,cols = rr)
dev.off()



### Contribution of each batch to each cluster
for(i in runbatch){
dge=dgelist[[i]]
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(Idents(dge))
  print(ncellsbatchcluster)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster)
  print(percentcellscluster)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellscluster)
  write.table(ncellscluster,paste0("ncellspercluster_",dataset[i],".txt"),quote=F,row.names=T,col.names=T,sep="\t")
}

### per-cell attributes for each cluster
pdf(file=paste0("plot/clusters_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=12)
for(indiv in runbatch){
  plotlist=list()
dge=dgelist[[indiv]]
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
}
dev.off()




###### Differentially-expressed markers 
res=paste0("RNA_snn_res.0.",c(6,4,3,4,1,1,1,1,1),"ordered")
markerslist=list()
for(i in runbatch){
dge=dgelist[[i]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
markerslist[[i]]=markers
write.table(markers,paste0("plot/",dataset[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
}

markerslist=list()
for(i in runbatch){
markers=read.table(paste0("plot/",dataset[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markerslist[[i]]=markers
}

plotlist=list()
for(i in runbatch){
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
ii=i-min(runbatch)+1
plotlist[[ii]]=cowplot::plot_grid(plotlist = p)
}

pdf(paste0("plot/markerstop.pdf"),height=2*round(size),width=1.8*ceiling(size))
plotlist
dev.off()

