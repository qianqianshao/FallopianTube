# 10.12.2020 by Qianyi
# individual clustering for 3rd large healthy FT: FT4

FT4: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
FT4 subject has large fibroids but fallopian tube is normal
(1) Fimbria
(2) Ampulla 
(3) Isthmus
for each fallopian tube, directly-merge the 3 segments together and do clustering

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

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium2","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3","Fimbria4","Ampulla4","Isthmus4")
run=c("NovaA-223","NovaA-223","NovaA-223","NovaA-237","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-295","NovaA-295","NovaA-295")
sample=c("747-YS","747-YS","747-YS","1073-NU","1427-NU","1427-NU","1427-NU","1407-NU","1407-NU","1714-YS","1714-YS","1714-YS")
part=c("Fimbria","Ampulla","Isthmus","Myometrium","Fimbria","Ampulla","Isthmus","FimAmp","Isthmus","Fimbria","Ampulla","Isthmus")
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus","FallopianTube2","FallopianTube2","FallopianTube2","FallopianTube3","FallopianTube3","FallopianTube4","FallopianTube4","FallopianTube4")
menopause=c("peri","peri","peri","post","pre","pre","pre","pre","pre","pre","pre","pre")
subject=c("Human1","Human1","Human1","Human2","Human3","Human3","Human3","Human4","Human4","Human5","Human5","Human5")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)

subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube","Uterus","FallopianTube2","FallopianTube3","FallopianTube4")

ft=grep("Fallopian",organ)
indivft=grep("Fallopian",organs)


# Load the raw dataset from CellRanger
runbatch=10:n
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
### Add Batch description
  dge[["Name"]] <- names[i]
  dge[["Part"]] <- part[i]
  dge[["Organ"]] <- organ[i]
  dge[["Run"]] <- run[i]
  dge[["Menopause"]] <- menopause[i]
  dge[["Subject"]] <- subject[i]
### Normalize data
dge<-NormalizeData(dge)
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



# merge all parts for each individual organ
dgealllist=list()

plot2=list()


for(indiv in 5:length(organs)){


### simply merge Seurat object instead of merging raw gene expression matrix
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


### Normalize data if I skipped this earlier
dge<-NormalizeData(dge)
GetAssayData(dge)["RPL11",1:5]
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



dgealllist[[indiv]]=dgeall
save(dgeall,file=paste0(subjectorgans[indiv],".Robj"))
print(table(gsub("_.*","",names(Idents(dge)))))

pdf(paste0("plot/organs_HVG_x0.2_y0.2.pdf"),height=6,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(plot2,cols=2)
dev.off()



########## Analysis for merged data of individual organ
for(i in indivft[4]){
dge=dgealllist[[i]]
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))
}

###### Determine top PCs
numPCs=c(10,0,12,18,11,12)
for(i in indivft[4]){
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
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,1.2,by=0.1))
dgeall=dge

dgealllist[[i]]=dgeall

print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))
}

for(i in indivft[4]){
  dgeall=dgealllist[[i]]
save(dgeall,file=paste0(subjectorgans[i],".Robj"))
}

###### order clusters for each individual organ
res=paste0("RNA_snn_res.0.",c(4,0,2,2,3))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()

for(resi in indivft[4]){

dge=dgealllist[[resi]]

Idents(dge)<-dge[[res[resi]]]
print(c(organs[resi],length(unique(Idents(dge)))))
levels=levels(Idents(dge))

aa=AverageExpression(dge)$RNA
centroid=log(AverageExpression(dge)$RNA+1)

dd=GetAssayData(dge)[,which(Idents(dge)==3)]
rr=GetAssayData(dge,slot="counts")[,which(Idents(dge)==3)]

cc=ExpMean(dd)
cc["RPL11"]

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=centroid[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
dim(t(as.matrix(tmp))) #[1]    16 29061
da <- dist(t(as.matrix(tmp))[,1:29061], method = "euclidean")
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
if(resi==5){
  levels=levelss[[resi]][c(16:14,9:13,8,6,7,1:3,5,4)]
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
for(resi in indivft[4]){
  dgeall=dgealllist[[resi]]
save(dgeall,file=paste0(subjectorgans[resi],".Robj"))
}


### Comparison with clusters of individual segments
for(i in indivft[4]){
  dge=dgealllist[[i]]
ncellsindivclusters=table(Idents(dge),dge@meta.data$indivclusters)
  write.table(ncellsindivclusters,paste0("plot/",organs[i],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")



for(label in c("Part")){
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
  write.table(ncellscluster,paste0("plot/",organs[i],"ncellspercluster_",label,".txt"),quote=F,row.names=T,col.names=T,sep="\t")
}


# FT4: 
# -> 16 re-ordered clusters, flipped the order to be consistent with the order of Fimbria
# subclustering for secretory cluster 2-3

subsets=list(2:3)
subsetsname=c("2-3")

i=1

cc=subsets[[i]]
ccname=subsetsname[i]

dgefile=paste0("plot/C",ccname,"_")

dge=subset(dgeall,RNA_snn_res.0.3ordered %in% cc)
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
numPCs=5

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

table(dge@meta.data$RNA_snn_res.0.6,dge@meta.data$RNA_snn_res.0.7)
table(dge@meta.data$RNA_snn_res.0.6,dge@meta.data$RNA_snn_res.0.3)


###### order clusters for each individual organ
res=paste0("RNA_snn_res.0.",1)

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
}

### save dge with ordered clusters

save(dge,file=paste0(subjectorgans[1],"_C",ccname,".Robj"))


# merge the secretory subtype of FT1 and FT4 together and cross-tabulate the individual subcluster centroids
# make the subclustering order of secretory cell type to be consistent with Reclustered C1-4_1-6 of FT1
### load object for each secretory subset of FT1 and FT4
dgelist=list()
load(file=paste0(subjectorgans[1],"_C1-4.Robj")) # included ciliated cluster 1
dgelist[[1]]=dge
dge1=dgelist[[1]]
dge1=subset(dge1,RNA_snn_res.0.1ordered %in% c(2:6))
dge1 # 24508 features across 2122 samples 
dgelist[[1]]=dge1
load(file=paste0(subjectorgans[5],"_C2-3.Robj")) # did not include ciliated cluster 1
dgelist[[5]]=dge

dge1=dgelist[[1]]
dge12=dgelist[[5]]
dge2=merge(dge1,dge12)
dge=dge2
dgeall=dge

### Highly variable genes
dge<-FindVariableFeatures(dge)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge), 10)
plot1 <- VariableFeaturePlot(dge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000

all="FT1FT4_SubsetSecretory"
save(dgeall,file=paste0(all,".Robj"))

res=c(paste0("RNA_snn_res.0.",c(1,0,0,0,1),"ordered") )
reps=c(1,5)

hvg.union=NULL
for(i in indivft){
hvg.union=unique(c(hvg.union,VariableFeatures(dgelist[[i]])))
}
length(hvg.union) # 3149


genes=hvg.union
genelabel="HVG"

### order cells by batch first, then by clusters of each batch
blockident=NULL
for(i in reps){
  tmp=paste(organs[i],Idents(dgelist[[i]]),sep="_")
  names(tmp)=names(Idents(dgelist[[i]]))
  blockident=c(blockident,tmp)
}
blockident=blockident[names(Idents(dgeall))]

### Clusters ordered first by batches, then by res
batch=organs[reps]
nbatch=length(batch)
ncluster=NULL
for(i in reps){
  ncluster=c(ncluster,length(unique(Idents(dgelist[[i]]))))
}
ncluster 
clusters=list()
for(i in reps){
  clusters[[i]]=levels(Idents(dgelist[[i]]))
}
levels2=NULL
for(bb in reps){
    cluster=clusters[[bb]]
    levels2=c(levels2,paste(organs[bb],cluster,sep="_"))
}
levels2=levels2[which(levels2 %in% unique(blockident))]
levels=levels2

ident=factor(blockident,levels=levels)

### Calculate correlation for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
dge=dgeall
Idents(dge)<-ident
dge$indivftsubclusters<-ident

dgeall=dge

save(dgeall,file=paste0(all,".Robj"))

centroid=log(AverageExpression(dge)$RNA+1)

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

### load cluster centroid rank correlation using HVG
data.use=read.table(paste0("plot/",all,"_indivft_Centroid_rho_",genelabel,".txt"),header=T,row.names=1)
colnames(data.use)=rownames(data.use)
levels=rownames(data.use)
batch=organs[indivft]

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
