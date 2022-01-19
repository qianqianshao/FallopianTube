# 7.8.2021 by Qianyi
# individual clustering for 4th large healthy FT: FT6

Sample3 (2518-AJ)  - included FT6, ovary, 2 Uterus parts from 1 caucasian woman, sequenced in 1 batch
  FT6
for each fallopian tube, do clustering


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



### load individual replicate clustering
i=16
  load(file=paste0(dataset[i],".Robj"))
  dge

indiv=9
  load(file=paste0(subjectorgans[indiv],".Robj")) 
  dgeall # same as above given FT6 is a whole FT without separate segments

### Add Batch description
  dge[["Name"]] <- names[i]
  dge[["Part"]] <- part[i]
  dge[["Organ"]] <- organ[i]
  dge[["Run"]] <- run[i]
  dge[["Menopause"]] <- menopause[i]
  dge[["Subject"]] <- subject[i]

  dge=RenameCells(dge, add.cell.id = names[i])


### Normalize data
dge<-NormalizeData(dge)
### Highly variable genes
dge<-FindVariableFeatures(dge)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge), 10)
plot1 <- VariableFeaturePlot(dge)
ii=1;plot2=list()
plot2[[ii]] <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(length(VariableFeatures(dge))) # 2000
### Scale data 
dge <- ScaleData(dge,features=rownames(dge))
  print(i)
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)


pdf(paste0("plot/HVG_x0.2_y0.2.pdf"),height=3,width=8)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(plot2,cols=3)
dev.off()


print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))

###### Determine top PCs
numPCs= c(6,8,8,8,7,8,8,7,8,
  0,0,0,
  11,0,0,
  13,14,  13,14,
  11,9)
pdf(paste("plot/PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
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


pdf(paste("plot/PCA_Variablel_variation_topPCs_heatmap.pdf",sep=""),height=15,width=10)
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


print(c( length(unique(dge$RNA_snn_res.0.1)),length(unique(dge$RNA_snn_res.0.2)),length(unique(dge$RNA_snn_res.0.3)),length(unique(dge$RNA_snn_res.0.4)),length(unique(dge$RNA_snn_res.0.5)),length(unique(dge$RNA_snn_res.0.6)),length(unique(dge$RNA_snn_res.0.7)),length(unique(dge$RNA_snn_res.0.8)),length(unique(dge$RNA_snn_res.0.9)),length(unique(dge$RNA_snn_res.1)),length(unique(dge$RNA_snn_res.1.1)),length(unique(dge$RNA_snn_res.1.2)) ))



###### order clusters for each individual organ
res=paste0("RNA_snn_res.0.",c(6,4,3,4,1,1,1,1,1,
  0,0,0,
  4,0,0,
  6,0,  2,2,
  1,5
))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
resi=i=16
Idents(dge)<-dge[[res[resi]]]
print(c(organ[resi],length(unique(Idents(dge)))))
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
if(resi==16){
  levels=levelss[[resi]][c(1:2,8:5,3:4,16:9)]
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

save(dge, file=paste0(dataset[i],".Robj"))

dgeall=dge
save(dgeall,file=paste0(subjectorgans[indiv],".Robj"))



table(Idents(dge))

   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
 298  281 1184  748   97   50  209   60  265  920  580   50  315  124 1421   57 



### PCA and tSNE for ordered clusters of each individual replicate
plotlist=list()
i=5
dge=dgelist[[i]]
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0.pdf",height=8,width=9)
multiplot(plotlist,cols = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(markers,paste0("plot/",organ[1],"_",numPCs[i],"PCs_",res,"_mindiff0.2_logfc2fold_4.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")


pdf(file=paste0("plot/",numPCs[i],"PCs_clusters_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=12)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
}
dev.off()
