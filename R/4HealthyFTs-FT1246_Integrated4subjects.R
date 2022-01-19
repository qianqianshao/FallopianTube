# 7.10.2021 by Qianyi
# integrated 4 large healthy FTs: FT1, FT2, FT4 and FT6.

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
merge 4 large healthy FTs - FT1, FT2, FT4, FT6 (10 datasets) together and do clustering
CCA-all genes: select HVG for each subject before integration; integrate all genes


R

library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)
library(RColorBrewer)
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



### 1. directly-merged 4 large healthy FTs (FT1, FT2, FT4, FT6)
load(file=paste0(all,".Robj"))
dge1=dgeall

### separate 4 healthy FTs 
dge=dgeall
genes.use=rownames(dge)
length(genes.use) # [1] 31594

datalists=SplitObject(dge,split.by="Organ")
names(datalists)=paste0("FT",c(1,2,4,6))
datalists


datalist <- lapply(X = datalists, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})




### Integration 
k.filter <- min(200, min(sapply(datalist, ncol)))
k.filter # revised k.filter based on the number of cells in the smallest dataset
anchors <- FindIntegrationAnchors(object.list = datalist,k.filter=k.filter, dims = 1:30) # anchor.features =length(genes.use),
integrated <- IntegrateData(anchorset = anchors, features.to.integrate = genes.use, dims = 1:30)


### Run PCA using integrated data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
dge <- integrated
dge <- ProjectDim(dge)

save(dge,file=paste0(all,"_CCA-4subjects-allgenes.Robj"))

hvg1=VariableFeatures(dge1)
hvg2=VariableFeatures(dge)
length(intersect(hvg1,hvg2)) # [1] 1626
# overlapped 1.6k HVG with directly-merged data

###### Determine top PCs
numPCs=10;i=1
pdf(paste("plot/CCA-4subjects-allgenes_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
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
id=read.table(paste0("plot/FallopianTube124_CCA3-allgenes_ident.txt"),header=F,row.names=1)
id1=id[,1]
names(id1)=rownames(id)
dge$ft124CCA3all = id1[names(Idents(dge))]


id=read.table(paste0("plot/FallopianTube124_CCA3-allgenes_ident21.txt"),header=F,row.names=1)
id1=id[,1]
names(id1)=rownames(id)
dge$ft124CCA3all21 = id1[names(Idents(dge))]



print(c( length(unique(dge$integrated_snn_res.0.1)),length(unique(dge$integrated_snn_res.0.2)),length(unique(dge$integrated_snn_res.0.3)),length(unique(dge$integrated_snn_res.0.4)),length(unique(dge$integrated_snn_res.0.5)),length(unique(dge$integrated_snn_res.0.6)),length(unique(dge$integrated_snn_res.0.7)),length(unique(dge$integrated_snn_res.0.8)),length(unique(dge$integrated_snn_res.0.9)),length(unique(dge$integrated_snn_res.1)),length(unique(dge$integrated_snn_res.1.1)),length(unique(dge$integrated_snn_res.1.2)) ))




###### order clusters for each individual organ
res=c("integrated_snn_res.0.2","integrated_snn_res.0.3");
resi=1
resi=2


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


table(dge$healthyftclusters,cells.ident.ordered[names(dge$healthyftclusters)])
table(dge$ft124CCA3all,cells.ident.ordered[names(dge$healthyftclusters)])


### re-order based on previous ordering
if(resi==1){
levels=levelss[c(16,12:15,11,8:6,10,9,1:3,5,4)]
}
if(resi==2){
levels=levelss[c(19,15:17,14,18,13,8,7,9,10,12,11,6,1,5:2)]
}

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


table(dge$healthyftclusters,cells.ident.ordered[names(dge$healthyftclusters)])
table(dge$ft124CCA3all,cells.ident.ordered[names(dge$healthyftclusters)])


### save ordered cluster ID in dge object
which(unique(cells.ident.ordered)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]

dge$ft1246CCA4all <- Idents(dge)



### Add identify of FT6 supervised assigned to global FT124 clusters
id=read.table(paste0("plot/FallopianTube124_CCA3-allgenes_ident.txt"),header=F,row.names=1)
id1=id[,1]
names(id1)=rownames(id)
dge$ft124CCA3all = id1[names(Idents(dge))]

classify3=read.table("plot/FT6_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
id3=classify3[,19]
names(id3)=rownames(classify3)
id=c(id1,id3)
id=id[names(Idents(dge))]
id=as.factor(id)
id=factor(id,levels=1:17)
dge$ft6assignedtoft124CCA3all = id
dge$ft124CCA3allFT6assign <- dge$ft6assignedtoft124CCA3all
Idents(dge)<-dge$ft124CCA3allFT6assign

dgeall=dge
dge4=dge


dgeall=dge
save(dgeall,file=paste0(all,"_CCA-4subjects-allgenes.Robj"))



res=paste0("integrated_snn_res.0.",c(2,3),"ordered");resi=i=2;

Idents(dge)<-dge[[res[i]]]
write.table(Idents(dge),file=paste0("plot/",all,"_CCA4-allgenes_ident.txt"),row.names=T,col.names=F,quote=F,sep="\t")


avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0(dgefile,"Centroid-uncorrected.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0(dgefile,"Centroid-postCCAallgenes.txt"),row.names=T,col.names=T,quote=F,sep="\t")


res=paste0("integrated_snn_res.0.",c(2,3),"ordered");resi=i=2;

i=resi
### Comparison with clusters of individual segments
ncellsindivclusters=table(dge@meta.data$healthyftclusters,Idents(dge))
  write.table(ncellsindivclusters,paste0("plot/",organs[1],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")

table(dge$healthyftclusters,Idents(dge))

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
dge$ft124CCA3allFT6assign1<-id
Idents(dge) <- dge$ft124CCA3allFT6assign1
dgeall=dge

# remove doublet cluster 13 
dge=subset(dgeall,ft124CCA3allFT6assign !=13)
dgeall=dge
save(dgeall,file=paste0(all,"_CCA-4subjects-allgenes-NoC13.Robj"))


### rename 6 FTs to new names to submit data to public databases
# rename donors
old -> new
FT1,FT2 unchanged
FT4 -> FT3
FT6 -> FT4
FT3 -> FT5
FT7 -> FT6
healthy: old FT1246 -> new FT1-4 
disease: old FT37 -> new FT5-6 

id=as.character(dge$Organ)
names(id)=names(Idents(dge))
table(id)
 FallopianTube FallopianTube2 FallopianTube4 FT6 
         10527          25396          17453           6659 
id[which(id == "FallopianTube")]<-"FT1"
id[which(id == "FallopianTube2")]<-"FT2"
id[which(id == "FallopianTube4")]<-"FT3"
id[which(id == "FT6")]<-"FT4"
id=factor(id,levels=paste0("FT",1:4)) 
table(id)
id
  FT1   FT2   FT3   FT4
10527 25396 17453  6659

dge$Sample<-id

dgeall=dge

# rename dataset names
rename FimAmp3 to FimAmp5, Isthmus3 to Isthmus5
Fimbria4 to Fimbria3, Ampulla4 to Ampulla3, Isthmus4 to Isthmus3
FT6 to FT4
FT7 to FT5

id=as.character(dge$Name)
names(id)=names(Idents(dge))
table(id)
id[which(id == "Fimbria")]<-"Fimbria1"
id[which(id == "Ampulla")]<-"Ampulla1"
id[which(id == "Isthmus")]<-"Isthmus1"
id[which(id == "Fimbria4")]<-"Fimbria3"
id[which(id == "Ampulla4")]<-"Ampulla3"
id[which(id == "Isthmus4")]<-"Isthmus3"
id[which(id == "FT6")]<-"FT4"
id=factor(id,levels=c("Fimbria1","Ampulla1","Isthmus1","Fimbria2","Ampulla2","Isthmus2","Fimbria3","Ampulla3","Isthmus3","FT4")) 
table(id)
 Fimbria1  Ampulla1  Isthmus1 Fimbria2 Ampulla2 Isthmus2 Fimbria3 Ampulla3 
    1861     4944     3722     9332     8021     8043     8200     2181 
Isthmus3      FT4 
    7072     6659 

dge$Dataset<-id

dgeall=dge
save(dgeall,file=paste0(all,"_CCA-4subjects-allgenes.Robj")) 


# rename cell names
rename FimAmp3 to FimAmp5, Isthmus3 to Isthmus5
Fimbria4 to Fimbria3, Ampulla4 to Ampulla3, Isthmus4 to Isthmus3
FT6 to FT4
FT7 to FT5

id <- as.character(dge$Dataset)
cells <- rownames(dge@meta.data)
table(gsub("_.*","",cells))
table(id)
new <- paste0(id,gsub(".*_","_",cells))
table(gsub("_.*","",new))

dge=RenameCells(dge,new.names=new)

dgeall=dge
save(dgeall,file=paste0(all,"_CCA-4subjects-allgenes_Renamed.Robj"))

cbind(1:ncol(dge@meta.data),colnames(dge@meta.data))
tmp=dge@meta.data[,c(1:4,82:83,81)]
colnames(tmp)[ncol(tmp)] <- "CellType"
dge@meta.data <- tmp
dgeall=dge
saveRDS(dgeall,file="FallopianTube1234_global.rds")
# submitted to cellxgene for visualization 




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
#  j=1
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
# saved as Figure S1D.
