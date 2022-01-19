# Analyzing each of the 3 individual segment of FT1 by Qianyi in Apr 2020
# 1st batch: 3 segments of fallopian tube from 1 peri-menopausal women
# individual clustering for each segment of FT1

R 

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

redblue100<-rgb(read.table(paste0('data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))

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


# Load the raw dataset from CellRanger
datalist=list()
for(i in 1:n){
	dir=paste0("/nfs/turbo/umms-hammou/10x/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/filtered_feature_bc_matrix/")
	data <- Read10X(data.dir = dir)
	datalist[[i]]=data
}
setlist=dgelist=list()
countcell=countgene=NULL
for(i in 1:n){
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
for(i in 1:length(dataset)){
	dgedata2=setlist[[i]]
write.table(dgedata2,file=paste0(dataset[i],"rawexp.txt"),col.names=T,row.names=T,sep="\t",quote=F)
}

### Per-cell attributes: Depth, Average number of genes, UMIs, %MT per cell
for(i in 1:length(dataset)){
dge=dgelist[[i]]
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))
}


# Visualize QC metrics as a violin plot
plot11=plot22=plot33=plot1=plot2=list()
for(i in 1:n){
dge=dgelist[[i]]
plot11[[i]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1)+ theme(legend.position = 'none')
plot22[[i]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+ theme(legend.position = 'none')
plot33[[i]]=VlnPlot(dge, features = "percent.mt",pt.size=-1)+ theme(legend.position = 'none')
plot1[[i]] <- FeatureScatter(dge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2[[i]] <- FeatureScatter(dge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
}
pdf(paste0("plot/FeatureScatter.pdf"),height=4.5,width=20)
multiplot(plot1,cols=4)
multiplot(plot2,cols=4)
dev.off()
pdf(paste0("plot/VlnPlot.pdf"),height=5,width=12)
multiplot(plot11,cols=4)
multiplot(plot22,cols=4)
multiplot(plot33,cols=4)
dev.off()


plot2=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
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
plot2[[i]] <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000
### Scale data 
dge <- ScaleData(dge,features=rownames(dge))
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)
dgelist[[i]]=dge
save(dge,file=paste0(dataset[i],".Robj"))
}


pdf(paste0("plot/HVG_x0.2_y0.2.pdf"),height=6,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(plot2,cols=2)
dev.off()


###### Determine top PCs
numPCs= c(6,8,8)

pdf(paste("plot/dge_PCA_Variablel_variation.pdf",sep=""),height=6,width=12)
par(mfrow=c(2,4),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plot(Stdev(dge,reduction="pca")[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
legend("topright",legend=dataset[i],cex=1.5)
abline(v=numPCs[i]+0.5,lty=2)
}
### density plot of Eigenvalue
for(i in 1:length(dataset)){
dge=dgelist[[i]]
eigenvalue=Stdev(dge,reduction="pca")[numPCs[i]]
print(eigenvalue)
plot(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(Stdev(dge,reduction="pca")),col="black")
lines(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.55,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=dataset[i],cex=1.5)
}
dev.off()
[1] 4.591819
[1] 3.708508
[1] 4.104249



###### Using top PCs
for(i in 1:length(dataset)){
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


for(i in 1:length(dataset)){
dge=dgelist[[i]]
print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))
}



###### order clusters for each dataset
res=paste0("RNA_snn_res.0.",c(6,4,3))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
for(resi in 1:length(dataset)){
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

### get order of seriation
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]
if(resi==2){
  levels=c(levelss[[resi]][1:3],rev(levelss[[resi]][4:7]),levelss[[resi]][8:11])
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
which(unique(cells.ident)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
# save the dge file
dgelist[[resi]]=dge
}

### save dge with ordered clusters
for(i in 1:length(dataset)){
dge=dgelist[[i]]
save(dge,file=paste0(dataset[i],".Robj"))
}


## double-check if I need to reverse the cluster ID orders
### put SPG (Zbtb16) / Leydig (Star) as beginning, and Myoid (Myh11) as end
# 4.14.2020 known markers from Nicole
knownmarkers=c("KRT7","PAX8",
  "CAPS","CCDC17","CCDC78","FOXJ1",
  "CLU",
  "KRT17","KRT23",
  "HLA-DQA1","HLA-DPA1","HLA-DPB1",
  "ALDH1A1","ALDH3B2","CDKN1A",
  "MKI67","MCM7","TK1","STMN1",
  "FANCD2","FANCL","MSH2",
  "HMGB2","SMC1A",
  "PTBP1","ZPR1","PRPF38A",
  "RGS16",
  "TIMP3","SPARC",
  "NR2F2","PDGFRA","TCF21",
  "ACTA2","MYH11",
  "KIT")
plotlist=list()
for(i in 1:length(dataset)){
gene=knownmarkers[which(knownmarkers %in% rownames(dgelist[[i]]))]
plotlist[[i]]=VlnPlot(dgelist[[i]],gene,ncol=6,pt.size=-1,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0_knownmarkers_Violin.pdf",height=10,width=15)
print(plotlist)
dev.off()
plotlist=list()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
for(i in 1:length(dataset)){
gene=knownmarkers[which(knownmarkers %in% rownames(dgelist[[i]]))]
p <- FeaturePlot(dgelist[[i]],gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p)
}
pdf("plot/clusters_ordered0_knownmarkers_Feature.pdf",height=10,width=10)
print(plotlist)
dev.off()

# 4.29.2020 known markers from Nicole
knownmarkers=c("SLMAP","ACTG2","MYLK","DES",
"CSPG4","PDGFRB","TBX18",
"ACTA2","MYH11",
"GLI1","ITGB1","CD44",
"KDR","PECAM1",
"PLIN2","PPARG","FABP1",
"NFIB","RUNX1","BRAF",
"COL4A1","POSTN",
"COL1A1","COL14A1","DKK3","TBX20","PDGFRA","VIM",
"DCN","SOX9","TCF21",
"CASP3","THY1","CD34",
"BMP4","APOE",
"CD8A","RUNX3","CD2",
"PTPRC","CD68",
"CD163","CD86",
"TINAGL1")
#Gene Symbol
NG2 CSPG4
CD29  ITGB1
VEGFR KDR
CD31  PECAM1
ADRP  PLIN2
HCLS  BRAF
SCA1  CASP3
CD45  PTPRC
GetAssayData(dge)["DES",1:5]
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
# FABP1 TBX20
plotlist=list()
for(i in 1:length(dataset)){
  gene=knownmarkers[which(knownmarkers %in% rownames(dgelist[[i]]))]
  plotlist[[i]]=VlnPlot(dgelist[[i]],gene,ncol=5,pt.size=-1,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0_knownmarkers2_Violin.pdf",height=15,width=12.5)
print(plotlist)
dev.off()
plotlist=list()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
for(i in 1:length(dataset)){
gene=knownmarkers[which(knownmarkers %in% rownames(dgelist[[i]]))]
p <- FeaturePlot(dgelist[[i]],gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p)
}
pdf("plot/clusters_ordered0_knownmarkers2_Feature.pdf",height=10,width=10)
print(plotlist)
dev.off()





## double-check if I need to flip PCs, or flip tSNE
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
pdf("plot/clusters_ordered0.pdf",height=4,width=18)
multiplot(plotlist2,cols = 4)
multiplot(plotlist3,cols = 4)
multiplot(plotlist4,cols = 4)
multiplot(plotlist5,cols = 4)
multiplot(plotlistt,cols = 4)
dev.off()



### Contribution of each batch to each cluster
for(i in 1:length(dataset)){
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
for(indiv in 1:length(dataset)){
  plotlist=list()
dge=dgelist[[indiv]]
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
}
dev.off()



###### Differentially-expressed markers 
res=paste0("RNA_snn_res.0.",c(6,4,3),"ordered")
markerslist=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
markerslist[[i]]=markers
write.table(markers,paste0("plot/",dataset[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
}

markerslist=list()
for(i in 1:length(dataset)){
markers=read.table(paste0("plot/",dataset[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markerslist[[i]]=markers
}

plotlist=list()
for(i in 1:length(dataset)){
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



###### Fraction of cells expression each gene in each cluster
# the counts of Non-0 values for each gene
res=paste0("RNA_snn_res.0.",c(6,4,3),"ordered")
for(resi in 1:length(dataset)){
dge=dgelist[[resi]]
genepct=matrix(,nrow(dge),length(unique(Idents(dge))))
ident=Idents(dge)[colnames(dge)]
    print(which(names(ident) != colnames(dge))) # make sure nothing
rownames(genepct)=rownames(dge)
colnames(genepct)=levels(ident)
for(i in 1:length(unique(ident))){
    id=levels(ident)[i]
    genepct[,i]=apply(GetAssayData(dge)[,which(ident==id)],1,function(x) sum(x!=0))/table(ident)[levels(ident)][i]
    print(i)
}
write.table(genepct,paste0("plot/",dataset[resi],"_",res[resi],"_FracPosCells.txt"),quote=F,row.names=T,col.names=T,sep="\t")
centroid=read.table(paste0("plot/",dataset[resi],"_",res[resi],"_Centroid.txt"),header=T,row.names=1)
colnames(centroid)=gsub("^X","",colnames(centroid))
### find unique genes for each cluster not presenting in other clusters
geneu=matrix(0,0,3)
for(i in 1:length(unique(ident))){
  id=levels(ident)[i]
  tmp1=apply(GetAssayData(dge)[,which(ident==id)],1,function(x) sum(x!=0))/table(ident)[levels(ident)][i]
  tmp2=apply(GetAssayData(dge)[,which(ident!=id)],1,function(x) sum(x!=0))/sum(table(ident)[levels(ident)][-i])
  genetmp=names(tmp1)[which(tmp1>0 & tmp2==0)]
  if(length(genetmp)>1){
  tmppct=genepct[genetmp,id]
  tmppct=sort(tmppct,decreasing=T)
  genetmp=names(tmppct)
  }
  tmpu=cbind(rep(id,length(genetmp)),centroid[genetmp,id],genepct[genetmp,id])
  geneu=rbind(geneu,tmpu)
}
colnames(geneu)=c("Cluster","AvgExp","Frac")
print(table(geneu[,1])[levels(ident)])
write.table(geneu,paste0("plot/",dataset[resi],"_",res[resi],"_uniquegenes_clusters_4.2020.txt"),quote=F,row.names=T,col.names=T,sep="\t")

genel=matrix(0,0,5)
for(i in 1:length(unique(ident))){
  id=levels(ident)[i]
  tmp1=apply(GetAssayData(dge)[,which(ident==id)],1,function(x) sum(x!=0))/table(ident)[levels(ident)][i]
  tmp2=apply(GetAssayData(dge)[,which(ident!=id)],1,function(x) sum(x!=0))/sum(table(ident)[levels(ident)][-i])
  centroid2=apply(GetAssayData(dge)[,which(ident!=id)],1,function(x) ExpMean(x))
  genetmp=names(tmp1)[which(tmp1> tmp2+0.5)]
  if(length(genetmp)>1){
  tmppct=tmp1[genetmp]-tmp2[genetmp]
  tmppct=sort(tmppct,decreasing=T)
  genetmp=names(tmppct)
  }
  tmpu=cbind(rep(id,length(genetmp)),centroid[genetmp,id],centroid2[genetmp],tmp1[genetmp],tmp2[genetmp])
  genel=rbind(genel,tmpu)
}
colnames(genel)=c("Cluster","AvgExp","AvgExpOthers","Frac","FracOthers")
print(table(genel[,1])[levels(ident)])
write.table(genel,paste0("plot/",dataset[resi],"_",res[resi],"_diff50pct_4.2020.txt"),quote=F,row.names=T,col.names=T,sep="\t")

}
 


###### Rank correlation for each normalized centroid using HVG
### Cluster centroids
for(resi in 1:length(dataset)){
dge=dgelist[[resi]]
# for each cluster, calculate average normalized expression of each gene
centroid=log(AverageExpression(dge)$RNA+1)
write.table(centroid,paste0("plot/",dataset[resi],"_",res[resi],"_Centroid.txt"),quote=F,row.names=T,col.names=T,sep="\t")
}


### rank correlation using HVG for each cluster of each batch
pdf(file="plot/indiv_clusters_Centroid_RankedCorrelation_HVG.pdf",height=5.2,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(resi in 1:length(dataset)){
dge=dgelist[[resi]]
centroid=read.table(paste0("plot/",dataset[resi],"_",res[resi],"_Centroid.txt"),header=T,row.names=1)
colnames(centroid)=gsub("^X","",colnames(centroid))

if(length(unique(Idents(dge)))>1){

#dge=SetAllIdent(dge,id=paste0(res[resi],"order"))
print(c(dataset[resi],length(unique(Idents(dge)))))
levels=unique(levels(Idents(dge)))

cc=cor(as.matrix(centroid)[VariableFeatures(dge),],method="spearman")
dim(cc)
min(cc) # 0.036

### labeling
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])
row.lab=col.lab=gsub(".*_","",levels)

#ncluster=9
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep("white",length(sidecol[1,]))
sidecol[2,]=myBrewerPalette[as.numeric(gsub(".*-","",gsub(".*_","",levels)))]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")
col.use=redblue100

heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))

}
}
dev.off()

# merge the 3 segments of fallopian tube together and cross-tabulate the individual cluster centroids
# merge all partss for each individual organ
dgealllist=list()

plot2=list()


for(indiv in c(1:length(organs))){
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

pdf(paste0("plot/organs_HVG_x0.2_y0.2.pdf"),height=6,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(plot2,cols=2)
dev.off()


dgealllist[[indiv]]=dgeall
save(dgeall,file=paste0(subjectorgans[indiv],".Robj"))
table(gsub("_.*","",names(Idents(dge))))
}


## Rank cor
datalist=list()
pdf(file=paste("plot/organ_parts_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))

for(indiv in 1:length(organs)){


### using either union of HVG or top 50 markers for each cluster
res=paste0("RNA_snn_res.0.",c(6,4,3))
reps=grep(organs[indiv],organ)

hvg.union=NULL
for(i in 1:length(reps)){
hvg.union=unique(c(hvg.union,VariableFeatures(dgelist[[i]])))
}
length(hvg.union) # 3039

top50markers=NULL
for(i in 1:length(reps)){
markers=read.table(paste0("plot/",dataset[i],"_",res[i],"_mindiff0.2_logfc2fold_4.2020.txt"),header=T,row.names=1,stringsAsFactors=F)
markers %>% group_by(cluster) %>% top_n(50, avg_logFC)  -> top50
top50markers=unique(c(top50markers,unique(top50$gene)))
}
length(top50markers) # 648

genelist=list(hvg.union,top50markers)
genelabels=c("HVG","Top50Markers")


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

### Calculate correlation for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
dgeall=dgealllist[[indiv]]
dge=dgeall
Idents(dge)<-ident
dge$indivclusters<-ident

dgeall=dge
dgealllist[[indiv]]=dge
save(dgeall,file=paste0(subjectorgans[indiv],".Robj"))


centroid=log(AverageExpression(dge)$RNA+1)


for(g in 1:length(genelist)){
genelabel=genelabels[g]
genes=genelist[[g]]
print(length(genes))

cc=cor(as.matrix(centroid)[genes,],method="spearman")
dim(cc)
min(cc) 

data.use=cc[levels,levels]
datalist[[rep]]=data.use

### load cluster centroid rank correlation using HVG
data.use=read.table(paste0("plot/",subjectorgans[indiv],"_rep_Centroid_rho_",genelabel,".txt"),header=T,row.names=1)
colnames(data.use)=rownames(data.use)
levels=rownames(data.use)
batch=part[grep(organs[indiv],organ)]

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
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=.8,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))


}

dev.off()

