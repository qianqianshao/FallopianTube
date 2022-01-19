# 1.4.2021 by Qianyi
# Corrected batch effect for 3 large healthy FTs - FT1, FT2, and FT4 by CCA integration and did clustering
# Related to Figure 1, Table S3 and Table S6.

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
1. directly-merge 3 healthy fallopian tubes (9 datasets) together and do clustering
2. CCA integrate 3 subjects of FT1, FT2, FT4 and do clustering
CCA-all genes: select HVG for each subject before integration; integrate all genes



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



### 1. directly-merged 3 large healthy FTs (FT1, FT2, FT4)
load(file=paste0(all,".Robj"))
dge1=dgeall

### separate 3 healthy FTs - FT1,2,4
dge=dgeall
# removed non-detected genes
nUMIperGene=apply(dge@assays$RNA@data,1,sum)
genes.use=names(nUMIperGene)[which(nUMIperGene>0)]
length(nUMIperGene) # 31932
length(genes.use)    # 31281

dge@assays$RNA@counts=dge@assays$RNA@counts[genes.use,]
dge@assays$RNA@data=dge@assays$RNA@data[genes.use,]


datalists=SplitObject(dge,split.by="Organ")
names(datalists)=paste0("FT",c(1,2,4))
datalists


datalist <- lapply(X = datalists, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})



dgefile="plot/CCA3allgenes_"


### Integration 
k.filter <- min(200, min(sapply(datalist, ncol)))
k.filter # revised k.filter based on the number of cells in the smallest dataset
anchors <- FindIntegrationAnchors(object.list = datalist,k.filter=k.filter, dims = 1:30) # anchor.features =length(genes.use),
save(anchors,file=paste0(all,"_CCA-3subjects-anchors.Robj"))
integrated <- IntegrateData(anchorset = anchors, features.to.integrate = genes.use, dims = 1:30)


### Run PCA using integrated data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
dge <- integrated
dge <- ProjectDim(dge)

hvg1=VariableFeatures(dge1)
hvg2=VariableFeatures(dge)
length(intersect(hvg1,hvg2)) # [1] 1626
# overlapped 1.6k HVG with directly-merged data

###### Determine top PCs
numPCs=c(10,9,9,12);i=3
pdf(paste("plot/CCA-3subjects-allgenes_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(Stdev(dge,reduction="pca")[1:30],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
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

###### Using top PCs
### tSNE
  dge <- RunTSNE(dge, dims = 1:numPCs[i])
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, reduction.type = "pca", resolution = seq(0.1,1.2,by=0.1))



### Label clustering IDs from previous analysis
id=read.table(paste0("plot/",all,"_ident.txt"),header=T,row.names=1)
id1=id[,1]
names(id1)=rownames(id)
dge$ft124directmerge = id1[names(Idents(dge))]


print(c( length(unique(dge$integrated_snn_res.0.1)),length(unique(dge$integrated_snn_res.0.2)),length(unique(dge$integrated_snn_res.0.3)),length(unique(dge$integrated_snn_res.0.4)),length(unique(dge$integrated_snn_res.0.5)),length(unique(dge$integrated_snn_res.0.6)),length(unique(dge$integrated_snn_res.0.7)),length(unique(dge$integrated_snn_res.0.8)),length(unique(dge$integrated_snn_res.0.9)),length(unique(dge$integrated_snn_res.1)),length(unique(dge$integrated_snn_res.1.1)),length(unique(dge$integrated_snn_res.1.2)) ))



###### order clusters for each individual organ
res=c(paste0("RNA_snn_res.0.",c(1)),paste0("integrated_snn_res.0.",c(3,2)));resi=3


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


table(dge$indivftclusters,cells.ident.ordered[names(dge$indivftclusters)])
table(dge$ft124directmerge,cells.ident.ordered[names(dge$indivftclusters)])


### re-order based on previous ordering
levels=levelss[c(17,16,15,14,12,13,7,6,8,9,10,11,5,1,4,3,2)]

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


table(dge$ft124directmerge,cells.ident.ordered[names(dge$indivftclusters)])


### save ordered cluster ID in dge object
which(unique(cells.ident.ordered)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
# save the dge file
dgelist[[resi]]=dge

dge$ft124CCA3all=dge$integrated_snn_res.0.2ordered

save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes.Robj"))









res=c(paste0("RNA_snn_res.0.",c(1),"ordered"),paste0("integrated_snn_res.0.",c(3,2),"ordered"));resi=1
i=3
Idents(dge)<-dge[[res[i]]]
write.table(Idents(dge),file=paste0("plot/",all,"_CCA3-allgenes_ident.txt"),row.names=T,col.names=F,quote=F,sep="\t")


avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0(dgefile,"Centroid-uncorrected.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0(dgefile,"Centroid-postCCAallgenes.txt"),row.names=T,col.names=T,quote=F,sep="\t")


res=c(paste0("RNA_snn_res.0.",c(1),"ordered"),paste0("integrated_snn_res.0.",c(3,2),"ordered"));resi=1

i=resi
### Comparison with clusters of individual segments
ncellsindivclusters=table(dge@meta.data$indivftclusters,Idents(dge))
  write.table(ncellsindivclusters,paste0("plot/",organs[1],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")

table(dge$indivftclusters,Idents(dge))

dge$Name=factor(dge$Name,levels=c("Fimbria","Ampulla","Isthmus","Fimbria2","Ampulla2","Isthmus2","Fimbria4","Ampulla4","Isthmus4"))
for(label in c("Organ","Name")){
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
  write.table(ncellscluster,paste0(dgefile,"1_ncellspercluster_",label,".txt"),quote=F,row.names=T,col.names=T,sep="\t")
}

### PCA and tSNE for ordered clusters 
plotlist=list()
dge=dgelist[[i]]
Idents(dge) <- dge[[ordered]]
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/clusters_ordered0.pdf",height=8,width=9.5)
multiplot(plotlist,cols = 2)
dev.off()


### merge secretory subset into 1 group colored in grey
dge=dgeall
id=as.numeric(Idents(dge))
names(id)=names(Idents(dge))
id[which(id %in% c(2:6))]<-"2-6"
table(id)
id=factor(id,levels=c(1,"2-6",7:17)) 
table(Idents(dge),id)
table(id)
id
    1   2-6     7     8     9    10    11    12    13    14    15    16    17 
 2510 14002  5327  9250  3642  2043  3785  2352   198   474  5770  1713  2310 
dge$ft124CCA3all1<-id
Idents(dge) <- dge$ft124CCA3all1
dgeall=dge


### label secretory subclusters
load(file=paste0(all,"_C2-6.Robj"))
dge2=dge

dge=dgeall
id=as.character(Idents(dge))
names(id)=names(Idents(dge))
id2=as.numeric(Idents(dge2))
names(id2)=names(Idents(dge2))
id[names(id2)[which(id2 == 1)]]<-"2-6_1"
id[names(id2)[which(id2 == 2)]]<-"2-6_2"
id[names(id2)[which(id2 == 3)]]<-"2-6_3"
id[names(id2)[which(id2 == 4)]]<-"2-6_4"
id[names(id2)[which(id2 == 5)]]<-"2-6_5"
id[names(id2)[which(id2 == 6)]]<-"2-6_6"
id=factor(id,levels=c(1,paste("2-6",1:6,sep="_"),7:17)) 
table(Idents(dge),id)
table(id)
id
    1 2-6_1 2-6_2 2-6_3 2-6_4 2-6_5 2-6_6     7     8     9    10    11    12 
 2510   268   868  6522  5479   549   316  5327  9250  3642  2043  3785  2352 
   13    14    15    16    17 
  198   474  5770  1713  2310 
dge$ft124CCA3all2<-id
Idents(dge) <- dge$ft124CCA3all2
dgeall=dge


# separate cluster 13 into two clusters based on the separation on UMAP
Idents(dge) <- dge$ft124CCA3all1

id=as.character(Idents(dge))
names(id)=names(Idents(dge))
umap=dge@reductions$umap@cell.embeddings
umap13=umap[which(id==13),]
cell1=rownames(umap13)[which(umap13[,2]>=0)]
cell2=rownames(umap13)[which(umap13[,2]<0)]
id[cell1]<-"13_1"
id[cell2]<-"13_2"
id=factor(id,levels=c(1,"2-6",7:8,"13_1",9,10:12,14:15,"13_2",16:17)) 
table(Idents(dge),id)
table(id)
dge$ft124CCA3all13<-id
Idents(dge) <- dge$ft124CCA3all13
dgeall=dge


# remove doublet cluster 13 
dge=subset(dgeall,ft124CCA3all1!=13)
levels(dge$ft124CCA3all1)
dge$ft124CCA3all1=factor(dge$ft124CCA3all1,levels=levels(dge$ft124CCA3all1)[-9])
dge$ft124CCA3all2=factor(dge$ft124CCA3all2,levels=levels(dge$ft124CCA3all2)[-14])
dgeall=dge
save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes-NoC13.Robj"))


### label myofibroblast subclusters
load(file=paste0(all,"_C7-10.Robj"))
dge7=dge

dge=dgeall
Idents(dge) <- dge$ft124CCA3all1
id=as.character(Idents(dge))
names(id)=names(Idents(dge))
id2=as.numeric(Idents(dge7))
names(id2)=names(Idents(dge7))
id[names(id2)[which(id2 == 1)]]<-"7-10_1"
id[names(id2)[which(id2 == 2)]]<-"7-10_2"
id[names(id2)[which(id2 == 3)]]<-"7-10_3"
id[names(id2)[which(id2 == 4)]]<-"7-10_4"
id[names(id2)[which(id2 == 5)]]<-"7-10_5"
id=factor(id,levels=c(1,"2-6",paste("7-10",1:5,sep="_"),11,12,14:17)) 
table(Idents(dge),id)
table(id)
id
     1    2-6 7-10_1 7-10_2 7-10_3 7-10_4 7-10_5     11     12     14     15 
  2510  14002   5235   8898   3764    189   2176   3785   2352    474   5770 
    16     17 
  1713   2310  
dge$ft124CCA3all7<-id
Idents(dge) <- dge$ft124CCA3all7
dgeall=dge

save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes-NoC13.Robj"))



### relabel 12 clusters from the original 17 clusters
dge=dgeall
Idents(dge) <- dge$ft124CCA3all1
id=as.numeric(Idents(dge))
names(id)=names(Idents(dge))
table(Idents(dge))
table(id)
id=factor(id,levels=1:12) 
table(Idents(dge),id)
table(id)
id
    1     2     3     4     5     6     7     8     9    10    11    12 
 2510 14002  5327  9250  3642  2043  3785  2352   474  5770  1713  2310 
dge$ft124CCA3all11<-id
Idents(dge) <- dge$ft124CCA3all11
dgeall=dge
# 1:2 - epithelial (1 - ciliated; 2 - secretory)
# 3:8 - new stromal (3:6 - old stromal, myofibroblast; 7:8 - endothelial)
# 9:12 - immune cells


### relabel 12 clusters with 6 secretory subclusters
dge=dgeall
Idents(dge) <- dge$ft124CCA3all2
id=as.character(Idents(dge))
names(id)=names(Idents(dge))
id[which(id == "2-6_1")]<- "2_1" 
id[which(id == "2-6_2")]<- "2_2" 
id[which(id == "2-6_3")]<- "2_3" 
id[which(id == "2-6_4")]<- "2_4"  
id[which(id == "2-6_5")]<- "2_5" 
id[which(id == "2-6_6")]<- "2_6"
id[which(id == 7)]<- 3
id[which(id == 8)]<- 4
id[which(id == 9)]<- 5 
id[which(id == 10)]<- 6
id[which(id == 11)]<- 7
id[which(id == 12)]<- 8
id[which(id == 14)]<- 9
id[which(id == 15)]<- 10
id[which(id == 16)]<- 11
id[which(id == 17)]<- 12
id=factor(id,levels=c(1,paste("2",1:6,sep="_"),3:12)) 
table(Idents(dge),id)
table(id)
id
   1  2_1  2_2  2_3  2_4  2_5  2_6    3    4    5    6    7    8    9   10   11 
2510  268  868 6522 5479  549  316 5327 9250 3642 2043 3785 2352  474 5770 1713 
  12 
2310 
dge$ft124CCA3all22<-id
Idents(dge) <- dge$ft124CCA3all22
dgeall=dge

save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes-NoC13.Robj"))


### relabel 12 clusters with 4 ciliated subclusters
dge=dgeall
load(file=paste0("FallopianTube124_C1.Robj"))
dge11=dge

dge=dgeall
Idents(dge) <- dge$ft124CCA3all22
id=as.character(Idents(dge))
names(id)=names(Idents(dge))
id2=as.numeric(Idents(dge11))
names(id2)=names(Idents(dge11))
id[names(id2)[which(id2 == 1)]]<-"1_1"
id[names(id2)[which(id2 == 2)]]<-"1_2"
id[names(id2)[which(id2 == 3)]]<-"1_3"
id[names(id2)[which(id2 == 4)]]<-"1_4"
id=factor(id,levels=c(paste("1",1:4,sep="_"),paste("2",1:6,sep="_"),3:12)) 
table(Idents(dge),id)
table(id)
id
 1_1  1_2  1_3  1_4  2_1  2_2  2_3  2_4  2_5  2_6    3    4    5    6    7    8 
1466  387  513  144  268  868 6522 5479  549  316 5327 9250 3642 2043 3785 2352 
   9   10   11   12 
 474 5770 1713 2310  
dge$ft124CCA3all21<-id
Idents(dge) <- dge$ft124CCA3all21
dgeall=dge

save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes-NoC13.Robj"))

write.table(Idents(dgeall),file=paste0("plot/FallopianTube124_CCA3-allgenes_ident21.txt"),row.names=T,col.names=F,quote=F,sep="\t")



### PCA and tSNE for grouped secretory subset 
Idents(dge) <- dge$ft124CCA3all1
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[-1])
plotlist=list()
png("plot/clusters_ordered1.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()
png("plot/clusters_ordered11.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()
# visualize cluster 13 
Idents(dge) <- dge$ft124CCA3all1
png("plot/cluster13.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette[9])
plotlist[[2]]=PCAPlot(dge,c(1,3),cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette[9])
plotlist[[3]]=DimPlot(dge,reduction="umap",cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette[9])
plotlist[[4]]=TSNEPlot(dge,cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette[9])
multiplot(plotlist,cols = 2)
dev.off()
Idents(dge) <- dge$ft124CCA3all13
png("plot/cluster13_2.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6)
plotlist[[2]]=PCAPlot(dge,c(1,3),cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6)
plotlist[[3]]=DimPlot(dge,reduction="umap",cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6)
plotlist[[4]]=TSNEPlot(dge,cells=rownames(dge@meta.data)[which(dge$ft124CCA3all1==13)],pt.size=.8,label=TRUE,label.size=6)
multiplot(plotlist,cols = 2)
dev.off()
# After removing cluster 13
Idents(dge) <- dge$ft124CCA3all1
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])
plotlist=list()
png("plot/clusters_ordered1.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()
png("plot/clusters_ordered11.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()
# saved as Figure 1B.


# for 2 subclusters of cluster 13 subset
dge=subset(dgeall,ft124CCA3all1==13)
Idents(dge) <- dge$ft124CCA3all13
DefaultAssay(dge) <- "RNA"       
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
table(Idents(dge))
table(markers$cluster)
13_1 13_2 
  17    2 
write.table(markers,paste0(dgefile,"RNAassay_cluster13_mindiff0.2_logfc2fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.6))
table(markers$cluster)
write.table(markers,paste0(dgefile,"RNAassay_cluster13_min0.1_logfc1.6fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
gene=top2$gene
pdf("plot/cluster13_markers_Violin.pdf",height=2,width=6)
VlnPlot(dge,gene,ncol=4,pt.size=-1)
dev.off()




Idents(dge) <- dge$ft124CCA3all1
DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
table(markers$cluster)
write.table(markers,paste0(dgefile,"RNAassay1_mindiff0.2_logfc2fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
Sys.time()

table(dge$ft124CCA3all1)
1   2-6     7     8     9    10    11    12    13    14    15    16    17 
2510 14002  5327  9250  3642  2043  3785  2352   198   474  5770  1713  2310 


Idents(dge) <- dge$ft124CCA3all11
DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
table(markers$cluster)
write.table(markers,paste0(dgefile,"RNAassay1_renamedID_mindiff0.2_logfc2fold_6.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
Sys.time()

table(dge$ft124CCA3all11)
    1     2     3     4     5     6     7     8     9    10    11    12 
 2510 14002  5327  9250  3642  2043  3785  2352   474  5770  1713  2310 
table(markers$cluster)
  1   2   3   4   5   6   7   8   9  10  11  12 
318  69 136  85  79 140 144 179  35 166 115 242 


# check the number of markers overlapped among each cluster
levels=levels(dge$ft124CCA3all1)
n=length(levels)
cc=matrix(0,n,n)
for(i in 1:n){
    for(j in 1:n){
        cc[i,j]=length(intersect(markers$gene[which(markers$cluster==levels[i])],markers$gene[which(markers$cluster==levels[j])]))
    }
}
cc

markers=read.table(paste0(dgefile,"RNAassay1_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)
markers1=markers
# check the number of markers unique (non-overlapped) for each cluster
uniquelist=list()
levels=levels(dge$ft124CCA3all1)
n=length(levels)
for(i in 1:n){
  list1=markers$gene[which(markers$cluster==levels[i])]
  list2=markers$gene[which(markers$cluster!=levels[i])]
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
# check the number of unique markers not detected in other clusters
uniquepct=NULL
levels=levels(dge$ft124CCA3all1)
n=length(levels)
for(i in 1:n){
  ii=markers[which(markers$cluster==levels[i]),]
  ii=ii[which(ii$pct.2==0),]
  print(length(which(ii$pct.2==0)))
  uniquepct=rbind(uniquepct,ii)
}
# no >2FC20%diff marker also satisfy this: not detected in other clusters

### find unique genes for each cluster not presenting in other clusters
ident=Idents(dge)[colnames(dge)]
    print(which(names(ident) != colnames(dge))) # make sure nothing
genepct=matrix(,nrow(dge),length(unique(Idents(dge))))
rownames(genepct)=rownames(dge)
colnames(genepct)=levels(ident)
for(i in 1:length(unique(ident))){
    id=levels(ident)[i]
    genepct[,i]=apply(GetAssayData(dge)[,which(ident==id)],1,function(x) sum(x!=0))/table(ident)[levels(ident)][i]
    print(i)
}
write.table(genepct,paste0(dgefile,"1_FracPosCells.txt"),quote=F,row.names=T,col.names=T,sep="\t")
#centroid=read.table(paste0(dgefile,"1_Centroid-uncorrected.txt"))
#colnames(centroid)=gsub("X","",colnames(centroid))
#colnames(centroid)[2]<-"2-6"
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
write.table(geneu,paste0(dgefile,"1_uniquegenes_2.2021.txt"),quote=F,row.names=T,col.names=T,sep="\t")

  1 2-6   7   8   9  10  11  12  14  15  16  17 
384 463 145 105  37  26  69  95  55  59  40 161 

summary(c(as.numeric(geneu[,2])))
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
4.799e-06 6.991e-05 2.091e-04 4.874e-04 5.035e-04 1.066e-02

max(as.numeric(geneu[,2])) # 0.01066262
max(as.numeric(geneu[,3])) # 0.6589641
geneu1=geneu[which(as.numeric(geneu$AvgExp))]

# check overlapped cluster 13 markers
markers=read.table(paste0(dgefile,"RNAassay_cluster13_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)
markers2=markers
length(intersect(markers$gene[which(markers$cluster=="13_1")],markers1$gene[which(markers1$cluster==13)]))
length(intersect(markers$gene[which(markers$cluster=="13_2")],markers1$gene[which(markers1$cluster==13)]))
# 0 0

markers=read.table(paste0(dgefile,"RNAassay13_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)
markers13=markers
length(intersect(markers$gene[which(markers$cluster=="13_1")],markers1$gene[which(markers1$cluster==13)]))
length(intersect(markers$gene[which(markers$cluster=="13_2")],markers1$gene[which(markers1$cluster==13)]))
# 14 11
length(intersect(markers$gene[which(markers$cluster=="13_1")],markers2$gene[which(markers2$cluster=="13_1")]))
length(intersect(markers$gene[which(markers$cluster=="13_1")],markers2$gene[which(markers2$cluster=="13_2")]))
length(intersect(markers$gene[which(markers$cluster=="13_2")],markers2$gene[which(markers2$cluster=="13_1")]))
length(intersect(markers$gene[which(markers$cluster=="13_2")],markers2$gene[which(markers2$cluster=="13_2")]))
# 4 0 0 1

# check cluster 13 subset markers overlapped with unique non-overlapped markers
length(which(markers2$gene %in% unlist(uniquelist)))
geneou2=unlist(uniquelist)[which(unlist(uniquelist) %in% markers2$gene)]
cc=matrix(0,n,2)
for(i in 1:n){
        cc[i,1]=length(intersect(uniquelist[[i]],markers2$gene[which(markers2$cluster=="13_1")]))
        cc[i,2]=length(intersect(uniquelist[[i]],markers2$gene[which(markers2$cluster=="13_2")]))
  )
}
cc

length(which(markers13$gene[grepl(13,markers13$cluster)] %in% unlist(uniquelist)))
geneou13=unlist(uniquelist)[which(unlist(uniquelist) %in% markers13$gene[grepl(13,markers13$cluster)])]
cc=matrix(0,n,2)
for(i in 1:n){
        cc[i,1]=length(intersect(uniquelist[[i]],markers13$gene[which(markers13$cluster=="13_1")]))
        cc[i,2]=length(intersect(uniquelist[[i]],markers13$gene[which(markers13$cluster=="13_2")]))
}
cc
# visualize top2 markers for cluster 13 subset
dge=dgeall
Idents(dge) <- dge$ft124CCA3all13
DefaultAssay(dge) <- "RNA" 
cols=c(myBrewerPalette[1:4],gg_color_hue(2)[1],myBrewerPalette[c(5:8,10,11)],gg_color_hue(2)[2],myBrewerPalette[12:13])
markers2 %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
gene=top2$gene
pdf("plot/cluster13_markers2_Violin.pdf",height=2,width=9)
VlnPlot(dge,gene,ncol=4,pt.size=-1,cols=cols)
dev.off()
markers13[grepl(13,markers13$cluster),] %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
gene=top2$gene
pdf("plot/cluster13_markers13_Violin.pdf",height=2,width=9)
VlnPlot(dge,gene,ncol=4,pt.size=-1,cols=cols)
dev.off()
# visualize top2 markers for cluster 13 subset overlapped with unique genes
pdf("plot/cluster13_markers2_overlappedwithunique_Violin.pdf",height=4,width=14)
VlnPlot(dge,geneou2,ncol=6,pt.size=-1,cols=cols)
dev.off()
pdf("plot/cluster13_markers13_overlappedwithunique_Violin.pdf",height=6,width=14)
VlnPlot(dge,geneou13,ncol=6,pt.size=-1,cols=cols)
dev.off()


# check cluster 7-10_4 subcluster markers overlapped with unique non-overlapped markers
markers=read.table(paste0("plot/C7-10_6PCs_integrated_snn_res.0.1ordered_mindiff0.2_logfc2fold_RNAassay_1.2021.txt"),stringsAsFactors=F)
markers7=markers               
markers2=markers7[which(markers7$cluster==4),]
length(which(markers2$gene %in% unlist(uniquelist))) # 19
geneou2=unlist(uniquelist)[which(unlist(uniquelist) %in% markers2$gene)]
cc=matrix(0,n,5)
for(i in 1:n){
        cc[i,1]=length(intersect(uniquelist[[i]],markers7$gene[which(markers7$cluster==1)]))
        cc[i,2]=length(intersect(uniquelist[[i]],markers7$gene[which(markers7$cluster==2)]))
        cc[i,3]=length(intersect(uniquelist[[i]],markers7$gene[which(markers7$cluster==3)]))
        cc[i,4]=length(intersect(uniquelist[[i]],markers7$gene[which(markers7$cluster==4)]))
        cc[i,5]=length(intersect(uniquelist[[i]],markers7$gene[which(markers7$cluster==5)]))
}
cc
# compare between C7-10_4 and C14-15
dge=dgeall
Idents(dge) <- dge$ft124CCA3all7
DefaultAssay(dge) <- "RNA"        # should use this
markers4=FindMarkers(dge,"7-10_4",c(14,15),logfc.threshold = log(2),min.diff.pct=0.2,min.diff=0.1)
print(c(length(which(markers4$avg_logFC>0)),length(which(markers4$avg_logFC<0))))
geneou42=unlist(uniquelist)[which(unlist(uniquelist) %in% rownames(markers4)[which(markers4$avg_logFC>0)])]
length(geneou42)
cc=matrix(0,n,2)
for(i in 1:n){
        cc[i,1]=length(intersect(uniquelist[[i]],rownames(markers4)[which(markers4$avg_logFC>0)]))
        cc[i,2]=length(intersect(uniquelist[[i]],rownames(markers4)[which(markers4$avg_logFC<0)]))
}
cc
# visualize markers for cluster 7-10_4 subcluster overlapped with unique markers of global clusters
dge=dgeall
Idents(dge) <- dge$ft124CCA3all7
DefaultAssay(dge) <- "RNA" 
cols=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(7,"Set2")[1:5],brewer.pal(12,"Paired")[c(6:7,9:12)])
pdf("plot/C7-10_4_markers2_overlappedwithunique_Violin.pdf",height=8,width=12)
VlnPlot(dge,geneou2,ncol=5,pt.size=-1,cols=cols)
dev.off()
markers2[which(markers2$gene %in% geneou2),] %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
gene41=top5$gene
markers4$gene=rownames(markers4)
markers4[geneou42,] %>% top_n(20, avg_logFC)  -> top20
gene42=top20$gene
pdf("plot/C7-10_4_markers22_overlappedwithunique_Violin.pdf",height=8,width=12)
#VlnPlot(dge,gene41,ncol=5,pt.size=-1,cols=cols)
VlnPlot(dge,gene42,ncol=5,pt.size=-1,cols=cols)
dev.off()
geneou426=uniquelist[[2]][which(uniquelist[[2]] %in% rownames(markers4)[which(markers4$avg_logFC>0)])]
pdf("plot/C7-10_4_markers_overlappedwithunique2-6_Violin.pdf",height=2,width=12)
#VlnPlot(dge,gene41,ncol=5,pt.size=-1,cols=cols)
VlnPlot(dge,geneou426,ncol=5,pt.size=-1,cols=cols)
dev.off()
# Jun asked me to visualize unique markers in C1 and C2-6
dgefile="plot/CCA3allNoC13_"
markers=read.table(paste0(dgefile,"RNAassay1_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)
markers126=markers[which(markers$cluster %in% c(1,"2-6")),]
markers126[which(markers126$gene %in% unlist(uniquelist[1:2])),] %>% group_by(cluster) %>% top_n(20, avg_logFC)  -> top20
geneu126=top20$gene
pdf("plot/global_markers126unique_Violin.pdf",height=16,width=12)
#VlnPlot(dge,gene41,ncol=5,pt.size=-1,cols=cols)
VlnPlot(dge,geneu126,ncol=5,pt.size=-1,cols=cols)
dev.off()


Idents(dge) <- dge$ft124CCA3all1
DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.6))
table(markers$cluster)
write.table(markers,paste0(dgefile,"RNAassay1_min0.1_logfc1.6fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
Sys.time()

# all genes FC Frac
Idents(dge) <- dge$ft124CCA3all1
DefaultAssay(dge) <- "RNA"        # should use this
Sys.time()
markers=FindAllMarkers(dge,only.pos=FALSE,logfc.threshold = -Inf, min.diff.pct=-Inf, min.pct=-Inf,min.cells.feature=-Inf,min.cells.group=-Inf,return.tresh=-Inf)
table(markers$cluster)
write.table(markers,paste0(dgefile,"RNAassay1_all_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
Sys.time()
save(markers,file=paste0(dgefile,"RNAassay1_all_2.2021.Robj"))




######## Heatmap for all markers
dge=dgeall
dgefile="plot/CCA3allgenes_"

Idents(dge) <- dge$ft124CCA3all1
avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0(dgefile,"1_Centroid-uncorrected.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0(dgefile,"1_Centroid-postCCAallgenes.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid[1:5,1:3]

Idents(dge) <- dge$ft124CCA3all2
avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0(dgefile,"2_Centroid-uncorrected.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0(dgefile,"2_Centroid-postCCAallgenes.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid[1:5,1:3]


Idents(dge) <- dge$ft124CCA3all21
avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0(dgefile,"21_Centroid-uncorrected.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0(dgefile,"21_Centroid-postCCAallgenes.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid[1:5,1:3]



centroid=read.table(paste0(dgefile,"1_Centroid-postCCAallgenes.txt"))
colnames(centroid)=gsub("X","",colnames(centroid))

# removed cluster 13
centroid=centroid[,-9]
colnames(centroid)[2]<-"2-6"

### Genes Standardized Across Cell Types
# note: used this
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)


markers=read.table(paste0(dgefile,"RNAassay1_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)
table(Idents(dge))
table(markers$cluster)[levels(Idents(dge))]


}


DefaultAssay(dge) <- "RNA"
data=GetAssayData(dge)
data.use=data[markers$gene,]
save(data.use,file=paste0(dgefile1,"1_normalized_markers.Robj"))
write.table(data.use,paste0(dgefile1,"1_normalized_markers.txt"),row.names=T,col.names=T,quote=F,sep="\t")
DefaultAssay(dge) <- "integrated"
data=GetAssayData(dge)
data.use=data[markers$gene,]
save(data.use,file=paste0(dgefile1,"1_normalized-postCCAallgenes_markers.Robj"))
write.table(data.use,paste0(dgefile1,"1_normalized-postCCAallgenes_markers.txt"),row.names=T,col.names=T,quote=F,sep="\t")


### Visualize markers in heatmap across all cell types
genes=markers$gene
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




jpeg(file=paste0(dgefile,"1_centroid_std_markersall.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0(dgefile,"1_centroid_std_markersall2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
# saved as Figure 1C and Table S3.


# Differential gene expression compared across the fallopian tube segments
  Idents(dge)<-dge[["Part"]]
DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
print(table(Idents(dge)))
  print(table(markers$cluster))
write.table(markers,paste0(dgefile,"3segments_mindiff0.2_logfc2fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.6),min.pct=0)
  print(table(markers$cluster))
write.table(markers,paste0(dgefile,"3segments_min0_logfc1.6fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")

pdf(file=paste0("plot/3Segments_PerCellAttributes_ViolinPlot.pdf"),height=2.5,width=7.5)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
dev.off()

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
gene=top5$gene
plotlist=VlnPlot(dge,gene,ncol=5,pt.size=-1)
}
pdf("plot/segment_markers_Violin.pdf",height=6,width=8)
print(plotlist)
dev.off()

# differential expression across 3 segments for the same cell types
dgeall=dge
  Idents(dgeall)<-dgeall[["Part"]]
DefaultAssay(dgeall) <- "RNA"        # should use this

levels=levels(dgeall$ft124CCA3all1)
nmarker=NULL
plotlist=list()
for(j in levels){
  dge=subset(dgeall, ft124CCA3all1 == j)
  print(dge)
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
print(table(Idents(dge)))
  print(table(markers$cluster))
nmarker=rbind(nmarker,
  table(Idents(dge)),
  table(markers$cluster)
  )
write.table(markers,paste0(dgefile,"3segments_C",j,"_mindiff0.2_logfc2fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.6),min.pct=0)
  print(table(markers$cluster))
nmarker=rbind(nmarker,
  table(markers$cluster)
  )
write.table(markers,paste0(dgefile,"3segments_C",j,"_min0_logfc1.6fold_2.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
gene=top2$gene
plotlist[[j]]=VlnPlot(dge,gene,ncol=6,pt.size=-1)

markers=read.table(paste0(dgefile1,"3segments_C",j,"_min0_logfc1.6fold_2.2021.txt"))
markers %>% group_by(cluster) %>% top_n(15, avg_logFC)  -> top15
gene=top15$gene
plotlist[[j]]=VlnPlot(dge,gene,ncol=7,pt.size=-1)
png(paste0("plot/3segment_C",j,"_top15markers_Violin.png"),res=300,,height=4000,width=4200)
print(plotlist[[j]])
dev.off()

markers=read.table(paste0(dgefile1,"3segments_C",j,"_min0_logfc1.6fold_2.2021.txt"))
markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
gene=top5$gene
plotlist[[j]]=VlnPlot(dge,gene,ncol=5,pt.size=-1)
png(paste0("plot/3segment_C",j,"_top5markers_Violin.png"),res=300,,height=2000,width=3000)
print(plotlist[[j]])
dev.off()

}
}

pdf("plot/3segment_clusters1_markers_Violin.pdf",height=2,width=9)
print(plotlist)
dev.off()

pdf("plot/3segment_clusters1_top15markers_Violin.pdf",height=2,width=9)
print(plotlist)
dev.off()

for(j in levels){
png(paste0("plot/3segment_C",j,"_top15markers_Violin.png"),res=300,,height=4000,width=4200)
print(plotlist[[j]])
dev.off()
}


# differential expression across 3 segments across each of the 3 healthy subjects for the same cell types
dgeall=dge4
  Idents(dgeall)<-dgeall[["Part"]]
DefaultAssay(dgeall) <- "RNA"        # should use this

plotlist=list()
for(j in 1:2){
  dge=subset(dgeall, ft124CCA3all1 == levels[j])
  print(dge)
  gene=markerslist[[j]]
  for(subject in unique(dgeall$Organ)){
    dge1=subset(dge,Organ == subject)
if(j==1){
  plot=VlnPlot(dge1,gene,ncol=5,pt.size=-1)
plot=VlnPlot(dge1,gene,ncol=5,pt.size=-1)
png(paste0("plot/3segment_C",levels[j],"_",subject,"_top5markers_Violin.png"),res=300,,height=1200,width=4200)
print(plot)
dev.off()
} else{
   plot=VlnPlot(dge1,gene,ncol=4,pt.size=-1)
plot=VlnPlot(dge1,gene,ncol=4,pt.size=-1)
png(paste0("plot/3segment_C",levels[j],"_",subject,"_top4markers_Violin.png"),res=300,height=1200,width=3360)
print(plot)
dev.off()
 
}
}
}


levels=levels(dgeall$ft124CCA3all1)
nmarker=NULL
plotlist=list()
for(j in levels){
  dge=subset(dgeall, ft124CCA3all1 == j)
  print(dge)
  for(subject in unique(dgeall$Organ)){
    dge1=subset(dge,Organ == subject)
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
print(table(Idents(dge)))
  print(table(markers$cluster))
nmarker=rbind(nmarker,
  table(Idents(dge)),
  table(markers$cluster)
  )
write.table(markers,paste0(dgefile,"3segments_C",j,"_",subject,"_mindiff0.2_logfc2fold_11.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge1,only.pos=TRUE,logfc.threshold = log(1.6),min.pct=0)
  print(table(markers$cluster))
nmarker=rbind(nmarker,
  table(markers$cluster)
  )
write.table(markers,paste0(dgefile,"3segments_C",j,"_",subject,"_min0_logfc1.6fold_11.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
gene=top5$gene
plotlist[[j]]=VlnPlot(dge,gene,ncol=6,pt.size=-1)
png(paste0("plot/3segment_C",j,"_",subject,"_top5markers_Violin.png"),res=300,,height=2000,width=3000)
print(plotlist[[j]])
dev.off()
}
}
# saved as Table S6.


### Per-cell attributes for each cell type
levels=levels(Idents(dge))
for(j in levels){
print(c(mean(dge@meta.data$nFeature_RNA[which(dge$ft124CCA3all1==j)]),mean(dge@meta.data$nCount_RNA[which(dge$ft124CCA3all1==j)]),mean(dge@meta.data$percent.mt[which(dge$ft124CCA3all1==j)])))
}

Idents(dge) <- dge$ft124CCA3all1
pdf(file=paste0(dgefile,"1_PerCellAttributes_ViolinPlot.pdf"),height=7.5,width=6)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
dev.off()

Idents(dge) <- dge$ft124CCA3all13
pdf(file=paste0(dgefile,"13_PerCellAttributes_ViolinPlot.pdf"),height=7.5,width=6)
  plotlist=list()
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=c(myBrewerPalette[1:4],gg_color_hue(2)[1],myBrewerPalette[c(5:8,10,11)],gg_color_hue(2)[2],myBrewerPalette[12:13]))+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=c(myBrewerPalette[1:4],gg_color_hue(2)[1],myBrewerPalette[c(5:8,10,11)],gg_color_hue(2)[2],myBrewerPalette[12:13]))+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1,cols=c(myBrewerPalette[1:4],gg_color_hue(2)[1],myBrewerPalette[c(5:8,10,11)],gg_color_hue(2)[2],myBrewerPalette[12:13]))+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
dev.off()



### Visualize individual time and treatment
label=c("Organ","Part")
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


  Idents(dge)<-dge[[res[i]]]

DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(markers,paste0(dgefile,numPCs[i],"PCs_",res[i],"_RNAassay_mindiff0.2_logfc2fold_10.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")

dgeall=dge
save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes.Robj"))




# visualize known markers:
knownmarkers=c("CAPS","PAX8")
gene=knownmarkers
### violin plot
plotlist=list()
dge=dgeall
DefaultAssay(dge) <- "RNA"
  plotlist[[j]]=VlnPlot(dge,gene,ncol=2,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered0_knownmarkers_Violin_integrated.pdf"),height=2,width=5)
print(plotlist)
dev.off()
geom_violin
#adjust Violin plot 
knownmarkers="CAPS" 
knownmarkers="PAX8" 
VlnPlot(dge,knownmarkers,pt.size=-1,cols=myBrewerPalette)
data=data.frame(GENE=GetAssayData(dge)[knownmarkers,],ident=dge$ft124CCA3all,Organ=dge$Organ,Segment=dge$Name)
data$Segment=factor(data$Segment,levels=c("Fimbria","Ampulla","Isthmus","Fimbria2","Ampulla2","Isthmus2","Fimbria4","Ampulla4","Isthmus4"))
pdf(paste0(dgefile,knownmarkers,"_Violin_Organ.pdf"),height=6,width=4)
ggplot(data,aes(x=ident,y=GENE))+geom_violin(aes(fill=ident),scale="width")+scale_fill_manual(values = myBrewerPalette)+ theme(legend.position = "none") +
facet_wrap(~Organ,ncol=1)+ggtitle(knownmarkers)+ylab("")+xlab("Identity")
dev.off()
pdf(paste0(dgefile,knownmarkers,"_Violin_Segment.pdf"),height=6,width=10)
ggplot(data,aes(x=ident,y=GENE))+geom_violin(aes(fill=ident),scale="width")+scale_fill_manual(values = myBrewerPalette)+ theme(legend.position = "none") +
facet_wrap(~Segment)+ggtitle(knownmarkers)+ylab("")+xlab("Identity")
dev.off()


# 2.2.2021 visualize known markers
a.  Pan epithelial: KRT7, EPCAM
b.  Ciliated: CAPS, FOXJ1
c.  Non-ciliated: PAX8, OVGP1
d.  Myofibroblast: PDGFRA, POSTN, COL1A1
e.  Msenchymal stoma: DCN, VIM, NR2F2
f.  Smooth muscle: ACTA2, DES, MYH11
g.  Pericyte: MCAM, PDGFRB, CSPG4, TAGLN
h.  Endothelial: PECAM1, KDR, VWF, PROX1, PDPN
i.  Mast: Kit, TPSAB1, TPSB2, MS4A2
j.  Immune: PTPRC, CD3E, RUNX3, CD34, KLRC1, FOLR2, SDC1(CD138), CD68, CD163, NCAM1 (CD56), CD207, CLEC4C (CD303), CEACAM8 (CD66b), ITGAX (CD11c), ELANE
i.  Bcells: MS4A1 (CD20), JCHAIN
k.  Additional: S100A9, THY1, CD44, SUSD2, 

knownmarkers=c("KRT7","EPCAM",
  "CAPS","FOXJ1",
  "PAX8","OVGP1",
  "PDGFRA","POSTN","COL1A1",
  "DCN","VIM","NR2F2",
  "ACTA2","DES","MYH11",
  "MCAM","PDGFRB","CSPG4","TAGLN",
  "PECAM1","KDR","VWF","PROX1","PDPN",
  "KIT","TPSAB1","TPSB2","MS4A2",
  "PTPRC","CD3E","RUNX3","CD34","KLRC1","FOLR2","SDC1","CD68","CD163",
  "NCAM1","CD207","CLEC4C","CEACAM8","ITGAX","ELANE",
  "MS4A1","JCHAIN",
  "S100A9","THY1","CD44","SUSD2"
)
length(knownmarkers) # 49
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers
g1=gene

# re-ordered and changed markers on 2/16/2021
knownmarkers=c("EPCAM","FOXJ1","CAPS","KRT7","PAX8","OVGP1","DCN","COL1A1","CD34","PDGFRA","VIM","POSTN","NR2F2","DES","ACTA2","MYH11","TAGLN","PDGFRB","MCAM","CSPG4","KDR","PECAM1","VWF","PROX1","PDPN","SDC1","JCHAIN","MS4A1","KLRC1","RUNX3","CD3E","PTPRC","KIT","MS4A2","TPSB2","TPSAB1","FOLR2","CD68","CD163","ITGAX","CD207","CLEC4C","CEACAM8","ELANE")
length(knownmarkers) # 44
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers
g2=gene

### violin plot
dge=dgeall
Idents(dge) <- dge$ft124CCA3all1
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[-1])

DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=7,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered1_knownmarkers_Violin_RNA.pdf"),height=14,width=19)
print(plotlist)
dev.off()

DefaultAssay(dge) <- "integrated"
pdf(paste0(dgefile,"ordered1_knownmarkers_DotPlot_RNA.pdf"),height=5,width=15)
DotPlot(dge, features = rev(gene)) + RotatedAxis()
DotPlot(dge, features = rev(gene),cols=c("white","blue")) + RotatedAxis()
DotPlot(dge, features = rev(gene),cols=c("white","red")) + RotatedAxis()
dev.off()
# saved as Figure 1D.

# 2.14.21 notes
#4)  Violin plots for genes of interest:
a.  EMT specific for epithelial clusters: VIM, TIMP3, SPARC, COL1A1, TWIST1, TWIST2, S100A4, PRRX1, MMP2, ZEB1, ZEB2, FN1, CDH11, CDH2, SNAI1, SNAI2, MUC1, OCLN, LOX, CD44, PROM1 (CD133), POU5F1 (OCT-4), SUSD2
b.  Hormone receptors for all clusters, stromal, and epithelial: ESR, PGR, AR
c.  Additional for stromal: S100A9, THY1, CD44, SUSD2,

knownmarkers=c("ESR1","PGR","AR")
length(knownmarkers) # 3
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

dge=dgeall
Idents(dge) <- dge$ft124CCA3all1
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])

DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=3,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered1_hormonemarkers_Violin_RNA.pdf"),height=2,width=8)
print(plotlist)
dev.off()


### Stromal clusters C7-10
dge=subset(dgeall, ft124CCA3all1 %in% c(7:10))
Idents(dge) <- dge$ft124CCA3all1
myBrewerPalette=brewer.pal(12,"Paired")
dge7=dge

dge=dge7
Idents(dge) <- dge$ft124CCA3all1
myBrewerPalette=brewer.pal(12,"Paired")[-1]


plotlist=list()
plotlist[[1]]=PCAPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=1,label=TRUE,label.size=6,cols=myBrewerPalette)
}
pdf("plot/Stromal_clusters_ordered0.pdf",height=8,width=9.5)
multiplot(plotlist,cols = 2)
dev.off()
# remove labels
plotlist=list()
png("plot/Stromal_clusters_ordered1.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()


knownmarkers=c("ESR1","PGR","AR","S100A9", "THY1", "CD44", "SUSD2")
length(knownmarkers) # 7
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers


DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=4,pt.size=-1,cols=myBrewerPalette)
pdf(paste0("plot/Stromal_knownmarkers_Violin_RNA.pdf"),height=4,width=8)
print(plotlist)
dev.off()



# 5.3.2021 
3)  Violin plots for stromal clusters 3-8 (As far as I can tell we have not done violins for just these global clusters)
a.  KIT, ESR1, PGR, AR, THY1, ADH1B, PRRX1, SOX17, CD44, SUSD2, RUNX3, S100A9
b.  Pan epithelial: KRT7, EPCAM
c.  Ciliated: CAPS, FOXJ1
d.  Non-ciliated: PAX8, OVGP1
e.  Myofibroblast: PDGFRA, POSTN, COL1A1
f.  Msenchymal stoma: DCN, VIM, NR2F2
g.  Smooth muscle: ACTA2, DES, MYH11
h.  Pericyte: MCAM, PDGFRB, CSPG4, TAGLN
i.  Endothelial: PECAM1, KDR, VWF, PROX1, PDPN
j.  Mast: Kit, TPSAB1, TPSB2, MS4A2

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

### violin plot
dge=dgeall
Idents(dge) <- dge$ft124CCA3all11
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])

DefaultAssay(dge) <- "RNA"
  plotlist=VlnPlot(dge,gene,ncol=7,pt.size=-1,cols=myBrewerPalette)
pdf(paste0(dgefile,"ordered1_knownmarkers_Violin_RNA.pdf"),height=12,width=19)
print(plotlist)
dev.off()




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
  c("EPCAM","RUNX3")
 )

knownmarkers=EMT
i=1
plotlist=list()
dge=dgelist[[1]]
gene=knownmarkers[which(knownmarkers %in% rownames(dge@assays$RNA@data))]
length(gene) # 18

### separate healthy Vs disease organs
healthy=rownames(dge@meta.data)[which(dge$Organ!="FallopianTube3")]
disease=rownames(dge@meta.data)[which(dge$Organ=="FallopianTube3")]
print(c(length(healthy),length(disease)))
# 20269  5868
status=c(rep("Healthy",length(healthy)),rep("Disease",length(disease)))
names(status)=c(healthy,disease)
status=status[rownames(dge@meta.data)]
status=factor(status,levels=c("Healthy","Disease"))
dge$Status <- status

which(names(Idents(dge)) != names(status))
IdStatus<- paste(Idents(dge),status,sep="_")
levels=names(table(IdStatus))
names(IdStatus) = names(status)
IdStatus=factor(IdStatus,levels=levels)
table(IdStatus)
dge$IdStatus <- IdStatus

dgeall=dge
save(dgeall,file=paste0(all,"_CCA-3subjects-allgenes.Robj"))


dgeall=dge
dge11=dgeall
  DefaultAssay(dge11) <- "RNA"

datalist=SplitObject(dgeall,split.by="Status")

dge=dgeall
levels1=levels(Idents(dge))
Idents(dge) <- dge$IdStatus
levels=c(paste(levels1,"Healthy",sep="_"),paste(levels1,"Disease",sep="_"))
Idents(dge) <- factor(dge$IdStatus,levels=levels)
  DefaultAssay(dge) <- "RNA"
dge22=dge


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
  for(g in 1:length(markerslist)){
genes=markerslist[[g]]
for(id in levels(Idents(dge))){
  tmpgenes=tmp[[g]][[id]]
# Number of Non-0 cells
tt[g,id]=length(which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0))
# % of Non-0 cells
pct[g,id]=length(which(tmpgenes[1,]!=0 & tmpgenes[2,]!=0))/length(which(Idents(dge)==id))
}
}
print(tt)
print(pct)
}
}

# average expression
print(c(expMean(tmp1),expMean(tmp2),expMean(tmp3),expMean(tmp4)))
# average expression in Non0 Cells
print(c(expMean(tmp1[tmp1!=0]),expMean(tmp2[tmp2!=0]),expMean(tmp3[tmp3!=0]),expMean(tmp4[tmp4!=0])))
# Gini in all cells
print(c(ineq(tmp1),ineq(tmp2),ineq(tmp3),ineq(tmp4)))
}


### violin plot
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
  plotlist[[j]]=VlnPlot(dge,gene,ncol=5,pt.size=-1,cols=myBrewerPalette)
}
pdf(paste0(dgefile,"ordered0_EMTmarkers_Violin.pdf"),height=8,width=12.5)
print(plotlist)
dev.off()

jpeg(paste0(dgefile,"ordered0_EMTmarkers_Violin2.jpg"),res=300,height=3500,width=5000)
  VlnPlot(dge22,gene,ncol=5,pt.size=-1,cols=rep(myBrewerPalette[1:length(levels1)],2))
dev.off()

jpeg(paste0(dgefile,"ordered0_EMTmarkers_Violin_split.jpg"),res=300,height=3500,width=5000)
  VlnPlot(dge11,gene,ncol=5,split.by="Status",pt.size=-1)
dev.off()


### Ridge plot
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
  plotlist[[j]]=RidgePlot(dge,gene,ncol=5,cols=myBrewerPalette)
}
pdf(paste0(dgefile,"ordered0_EMTmarkers_Ridge.pdf"),height=8,width=12.5)
print(plotlist)
dev.off()


pdf(paste0(dgefile,"ordered0_EMTmarkers_Ridge2.pdf"),height=16,width=16)
  RidgePlot(dge22,gene,ncol=5,cols=rep(myBrewerPalette[1:length(levels1)],2))
dev.off()


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

### feature plot for both healthy and disease subjects
library(cowplot)
library(patchwork)
  for(genes in markerslist){
jpeg(paste0(dgefile,paste(genes,collapse="_"),".jpeg"),res=300,height=2000,width=3800)
plotlist=p=list()
for(j in 1:2){
  dge<-datalist[[j]]
  DefaultAssay(dge) <- "RNA"
  p[[j]] <-FeaturePlot(dge,genes,pt.size=.5, blend=TRUE,blend.threshold=0.25,combine=FALSE)
}
  plotlist=unlist(p,recursive=FALSE)
  wrap_plots(plotlist,ncol=4)
dev.off()
}

