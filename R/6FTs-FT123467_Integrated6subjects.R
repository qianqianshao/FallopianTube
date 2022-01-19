# 10.7.2021 by Qianyi
# integrated all 6 FTs
# Related to Figure 6, S1, S6, Table S7-S9. 

FT1: 3 segments of fallopian tube from 1 healthy perimenopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla
(3) Isthmus
FT2: 3 segments of fallopian tube from 1 healthy pre-menopausal women, sequenced in 1 batch
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT3: 2 segments of fallopian tube from 1 pre-menopausal women with disease, sequenced in 1 batch
(1) Fimbria and Ampulla - where the ovum is fertilized.
(2) Isthmus
FT4: 1 batch: 3 segments of fallopian tube from 1 pre-menopausal women, healthy FT, large fibroids
(1) Fimbria
(2) Ampulla 
(3) Isthmus
FT5: 1 batch: merged segments of fallopian tube from 1 pre-menopausal women, healthy FT
FT6: 1 whole fallopian tube from Sample3 (2518-AJ)  - included FT6, ovary, 2 Uterus parts from 1 caucasian woman, sequenced in 1 batch
FT7: 1 whole  fallopian tube from 1 pre-menopausal women with disease, sequenced in 1 batch
did individual clustering for each segment
for each of the 5 fallopian tubes, directly-merge the 3 segments together and do clustering
FT3 and FT7 are with disease, and FT5 has very few cells, so we decided to combine 4 healthy FTs - FT1, FT2, FT4, FT6. 
merge 6 FTs - FT1-7no5 (11 datasets) together and do clustering
CCA-integration: select HVG for each subject before integration; integrate all genes



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

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3","2518-AJ-1","2518-AJ-2","2518-AJ-3","2518-AJ-4","2788-NU-1","2788-NU-2","4360-YS-1")
dataset1=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3","2518-AJ-1_CACTACGA-ATCAGTCT","2518-AJ-2_CACGGTGA-TGTGACGA","2518-AJ-3_ATGGCTTG-CACAACAT","2518-AJ-4_CCTTCTAG-TCGTTGTA","2788-NU-1_TATCAGCC-AGGACGAA","2788-NU-2_TGGTCCCA-ACGCCAGA","Sample_4360-YS-1_GCTACAAA-AGGGCACG")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium1","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3","Fimbria4","Ampulla4","Isthmus4","Uterus2","FT5","Ovary1","FT6","Ovary2","Myometrium2","Endometrium2","Myometrium3","Endometrium3","FT7")
run=c("NovaA-223","NovaA-223","NovaA-223","NovaA-237","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-277","NovaA-295","NovaA-295","NovaA-295","NovaA-301","NovaA-301","NovaA-301","NovaZ-8","NovaZ-8","NovaZ-8","NovaZ-8","2788-NU","2788-NU","4360-YS")
sample=c("747-YS","747-YS","747-YS","1073-NU","1427-NU","1427-NU","1427-NU","1407-NU","1407-NU","1714-YS","1714-YS","1714-YS","1847-NU","1847-NU","1847-NU","2518-AJ","2518-AJ","2518-AJ","2518-AJ","2788-NU","2788-NU","4360-YS")
part=c("Fimbria","Ampulla","Isthmus","Myometrium","Fimbria","Ampulla","Isthmus","FimAmp","Isthmus","Fimbria","Ampulla","Isthmus","Uterus","FT","Ovary","FT","Ovary","Myometrium","Endometrium","Myometrium","Endometrium","FT")
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus","FallopianTube2","FallopianTube2","FallopianTube2","FallopianTube3","FallopianTube3","FallopianTube4","FallopianTube4","FallopianTube4","Uterus2","FallopianTube5","Ovary1","FT6","Ovary2","Sample3U2","Sample3U2","Uterus3","Uterus3","FT7")
menopause=c("peri","peri","peri","post","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","pre","unknown","unknown","pre")
subject=c("Human1","Human1","Human1","Human2","Human3","Human3","Human3","Human4","Human4","Human5","Human5","Human5","Human6","Human6","Human6","Sample3","Sample3","Sample3","Sample3","Uterus3","Uterus3","FT7")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)

subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube","Uterus","FallopianTube2","FallopianTube3","FallopianTube4","Uterus2","FallopianTube5","Ovary1","FallopianTube6","Ovary2","Sample3U2","Uterus3","FallopianTube7")

ft=grep("Fallopian|FT",organ)
indivft=grep("Fallopian|FT",subjectorgans)


ft=indivft[-5]
all="FallopianTube123467"


### 1. directly-merged 6 FTs
  load(file=paste0(all,".Robj"))
genes.use=rownames(dgeall)
table(dgeall$Organ)
FallopianTube FallopianTube2 FallopianTube3 FallopianTube4 FallopianTube6 
         10527          25396          15571          17453           6659 
FallopianTube7 
          2334 

### separate 6 FTs - FT1,2,3,4,6,7
datalists=SplitObject(dgeall,split.by="Organ")
datalists
names(datalists)=paste0("FT",c(1,2,3,4,6,7))
datalists=datalists[c(1,2,4,5,3,6)]


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

save(dge,file=paste0(all,"_CCA-6subjects-allgenes2.Robj"))


###### Determine top PCs
numPCs=11;i=1
pdf(paste("plot/CCA-6subjects-allgenes_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
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

save(dge,file=paste0(all,"_CCA-6subjects-allgenes_",numPCs[i],"PCs.Robj"))


i=1 # decided to use top 11 PCs

print(c( length(unique(dge$integrated_snn_res.0.1)),length(unique(dge$integrated_snn_res.0.2)),length(unique(dge$integrated_snn_res.0.3)),length(unique(dge$integrated_snn_res.0.4)),length(unique(dge$integrated_snn_res.0.5)),length(unique(dge$integrated_snn_res.0.6)),length(unique(dge$integrated_snn_res.0.7)),length(unique(dge$integrated_snn_res.0.8)),length(unique(dge$integrated_snn_res.0.9)),length(unique(dge$integrated_snn_res.1)),length(unique(dge$integrated_snn_res.1.1)),length(unique(dge$integrated_snn_res.1.2)) ))



### Comparison with clusters of individual segments
ncellsindivclusters=table(dge@meta.data$indivftclusters,Idents(dge))
  write.table(ncellsindivclusters,paste0("plot/",organs[1],"_ncellspercluster_indivclusters.txt"),quote=F,row.names=T,col.names=T,sep="\t")

table(dge$indivftclusters,Idents(dge))


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
for(label in c("Part","Organ","Name","Sample")){
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


  Idents(dge)<-dge[[res[i]]]

DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(markers,paste0(dgefile,numPCs[i],"PCs_",res[i],"_RNAassay_mindiff0.2_logfc2fold_10.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")

dgeall=dge
save(dgeall,file=paste0(all,"_CCA-6subjects-allgenes.Robj"))




############ Add FT3,FT6,FT7 identity supervised assigned to FT124 CCA3all to the FT1234 CCA4all object
dge4=dge

### Add identity of FT124 to FT1234 object
dge=dge4

id=read.table(paste0("plot/FallopianTube124_CCA3-allgenes_ident.txt"),header=F,row.names=1)
id1=id[,1]
names(id1)=rownames(id)
id1=id1[names(Idents(dge))]
dge$ft124CCA3all = id1
dge4=dge

### Add identify of FT3,FT6,FT7 supervised assigned to global FT124 clusters
id=as.numeric(dge4$ft124CCA3all)
names(id)=rownames(dge4@meta.data)
table(gsub("_.*","",names(which(is.na(id)))))
# FimAmp3      FT6      FT7 Isthmus3 
#    7634     6659     2334     7937 
grep("Fimbria1",names(which(is.na(id))),value=T)

classify0=read.table("plot/FT3_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
classify0[1:2,]
id0=classify0[,19]
names(id0)=rownames(classify0)
table(gsub("_.*","",names(id0))) # make sure the same as which(is.na(id))
id[names(id0)] <- id0
table(gsub("_.*","",names(which(is.na(id)))))

classify0=read.table("plot/FT6_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
classify0[1:2,]
id0=classify0[,19]
names(id0)=rownames(classify0)
table(gsub("_.*","",names(id0))) # make sure the same as which(is.na(id))
id[names(id0)] <- id0
table(gsub("_.*","",names(which(is.na(id)))))

classify0=read.table("plot/FT7_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
classify0[1:2,]
id0=classify0[,19]
names(id0)=rownames(classify0)
table(gsub("_.*","",names(id0))) # make sure the same as which(is.na(id))
id[names(id0)] <- id0
table(gsub("_.*","",names(which(is.na(id)))))

table(id)
    1     2     3     4     5     6     7     8     9    10    11    12    13 
 3500   172 14844  1651    78   254  8391 13903  4125  2337  4539  2888   404 
   14    15    16    17 
  846 10251  2686  7071 

dge$ft124CCA3allFT367assign <- id

dgeall=dge
dge4=dge


### Add identity of 6 secretory subclusters for C2-6
# use subclusters of FT124, supervised assigned FT3 to FT124 subclusters, and the FT1234 object 
### subclusters of FT124
i=1
ccname="2-6"

load(file=paste0("FallopianTube124_C",ccname,".Robj"))
dge
dge2=dge

dge=dge2
cells124=colnames(dge)
id124=as.numeric(Idents(dge))
names(id124)=cells124
table(id124)

### supervised assigned FT3 to FT124 subclusters
classify3=read.table(paste0("plot/FT3_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id3=classify3[,ncol(classify3)]
names(id3)=rownames(classify3)
cells3=names(id3)
table(id3)
### supervised assigned FT6 to FT124 subclusters
classify3=read.table(paste0("plot/FT6_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id6=classify3[,ncol(classify3)]
names(id6)=rownames(classify3)
cells6=names(id6)
table(id6)
### supervised assigned FT7 to FT124 subclusters
classify3=read.table(paste0("plot/FT7_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id7=classify3[,ncol(classify3)]
names(id7)=rownames(classify3)
cells7=names(id7)
table(id7)

id2=c(id124,id3,id6,id7)
cells=names(id2)
table(id2)
   1    2    3    4    5    6 
 341 1182 8207 5707 1022  540

### Add identity of 6 secretory subclusters for C2-6 to FT1234 object
dge=dge4
Idents(dge) <- dgeall$ft124CCA3allFT367assign

id=as.character(Idents(dge))
names(id)=names(Idents(dge))
id[names(id2)[which(id2 == 1)]]<-"2-6_1"
id[names(id2)[which(id2 == 2)]]<-"2-6_2"
id[names(id2)[which(id2 == 3)]]<-"2-6_3"
id[names(id2)[which(id2 == 4)]]<-"2-6_4"
id[names(id2)[which(id2 == 5)]]<-"2-6_5"
id[names(id2)[which(id2 == 6)]]<-"2-6_6"
id=factor(id,levels=c(1,paste("2-6",1:6,sep="_"),7:17)) 
table(Idents(dge),id)
table(id)
    1 2-6_1 2-6_2 2-6_3 2-6_4 2-6_5 2-6_6     7     8     9    10    11    12 
 3500   341  1182  8207  5707  1022   540  8391 13903  4125  2337  4539  2888 
   13    14    15    16    17 
  404   846 10251  2686  7071 
dge$ft124CCA3allFT367assign2<-id
Idents(dge) <- dge$ft124CCA3allFT367assign2
dgeall=dge
dge4=dge

### rename 6 FTs to new names
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
 FallopianTube FallopianTube2 FallopianTube3 FallopianTube4 FallopianTube6 
         10527          25396          15571          17453           6659 
FallopianTube7 
          2334 
id[which(id == "FallopianTube")]<-"FT1"
id[which(id == "FallopianTube2")]<-"FT2"
id[which(id == "FallopianTube3")]<-"FT5"
id[which(id == "FallopianTube4")]<-"FT3"
id[which(id == "FallopianTube6")]<-"FT4"
id[which(id == "FallopianTube7")]<-"FT6"
id=factor(id,levels=paste0("FT",1:6)) 
table(id)
id
  FT1   FT2   FT3   FT4   FT5   FT6 
10527 25396 17453  6659 15571  2334

dge$Sample<-id

dgeall=dge
dge4=dge


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
id[which(id == "FimAmp3")]<-"FimAmp5"
id[which(id == "Isthmus3")]<-"Isthmus5"
id[which(id == "Fimbria4")]<-"Fimbria3"
id[which(id == "Ampulla4")]<-"Ampulla3"
id[which(id == "Isthmus4")]<-"Isthmus3"
id[which(id == "FT6")]<-"FT4"
id[which(id == "FT7")]<-"FT6"
id=factor(id,levels=c("Fimbria1","Ampulla1","Isthmus1","Fimbria2","Ampulla2","Isthmus2","Fimbria3","Ampulla3","Isthmus3","FT4","FimAmp5","Isthmus5","FT6")) 
table(id)
Fimbria1 Ampulla1 Isthmus1 Fimbria2 Ampulla2 Isthmus2 Fimbria3 Ampulla3 
    1856     4901     3700     9314     8020     8033     8152     2180 
Isthmus3      FT4  FimAmp5 Isthmus5      FT6 
    7022     6560     7573     7899     2326 
 

dge$Dataset<-id

dgeall=dge


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
save(dgeall,file="FallopianTube123467_CCA-6subjects-allgenes_No13_Renamed.Robj")

cbind(1:ncol(dge@meta.data),colnames(dge@meta.data))
tmp=dge@meta.data[,c(1:4,82:83,87,85)]
colnames(tmp)[ncol(tmp)] <- "CellType"
dge@meta.data <- tmp
dgeall=dge
saveRDS(dgeall,file="FallopianTube123456_HealthyVsDiseased.rds")




### separate healthy Vs disease organs
healthy=rownames(dge@meta.data)[which(dge$Sample!="FT5" & dge$Sample!="FT6")]
disease=rownames(dge@meta.data)[which(dge$Sample=="FT5" | dge$Sample=="FT6")]
print(c(length(healthy),length(disease)))
#[1] 60035 17905
status=c(rep("Healthy",length(healthy)),rep("Disease",length(disease)))
names(status)=c(healthy,disease)
status=status[rownames(dge@meta.data)]
status=factor(status,levels=c("Healthy","Disease"))
dge$Status <- status

healthy=rownames(dge@meta.data)[which(dge$Sample!="FT5" & dge$Sample!="FT6")]
disease1=rownames(dge@meta.data)[which(dge$Sample=="FT5")]
disease2=rownames(dge@meta.data)[which(dge$Sample=="FT6")]
print(c(length(healthy),length(disease1),length(disease2)))
#[1] 60035 17905
status=c(rep("Healthy",length(healthy)),rep("Disease1",length(disease1)),rep("Disease2",length(disease2)))
names(status)=c(healthy,disease1,disease2)
status=status[rownames(dge@meta.data)]
status=factor(status,levels=c("Healthy","Disease1","Disease2"))
dge$Status2 <- status


dgeall=dge
dge4=dge


### relabel 12 clusters with 6 secretory subclusters
dge=dgeall
Idents(dge) <- dge$ft124CCA3allFT367assign2
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
    1   2_1   2_2   2_3   2_4   2_5   2_6     3     4     5     6     7     8 
 3500   341  1182  8207  5707  1022   540  8391 13903  4125  2337  4539  2888 
    9    10    11    12 
  846 10251  2686  7071 
dge$ft124CCA3allFT367assign22<-id
Idents(dge) <- dge$ft124CCA3allFT367assign22

dgeall=dge
dge4=dge


### Add identity of 4 ciliated subclusters for C1
# use subclusters of FT124, supervised assigned FT3 to FT124 subclusters, and the FT1234 object 
### subclusters of FT124
i=1
ccname="1"

load(file=paste0("FallopianTube124_C",ccname,".Robj"))
dge
dge2=dge

cells124=colnames(dge)
id124=as.numeric(Idents(dge))
names(id124)=cells124
table(id124)

### supervised assigned FT3 to FT124 subclusters
classify3=read.table(paste0("plot/FT3_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id3=classify3[,ncol(classify3)]
names(id3)=rownames(classify3)
cells3=names(id3)
table(id3)
### supervised assigned FT6 to FT124 subclusters
classify3=read.table(paste0("plot/FT6_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id6=classify3[,ncol(classify3)]
names(id6)=rownames(classify3)
cells6=names(id6)
table(id6)
### supervised assigned FT7 to FT124 subclusters
classify3=read.table(paste0("plot/FT7_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id7=classify3[,ncol(classify3)]
names(id7)=rownames(classify3)
cells7=names(id7)
table(id7)

id2=c(id124,id3,id6,id7)
cells=names(id2)
table(id2)

   1    2    3    4 
2297  400  580  223


### Add identity of 4 ciliated subclusters for C1 to FT1234 object
dge=dge4
Idents(dge) <- dgeall$ft124CCA3allFT367assign22

id=as.character(Idents(dge))
names(id)=names(Idents(dge))
id[names(id2)[which(id2 == 1)]]<-"1_1"
id[names(id2)[which(id2 == 2)]]<-"1_2"
id[names(id2)[which(id2 == 3)]]<-"1_3"
id[names(id2)[which(id2 == 4)]]<-"1_4"
id=factor(id,levels=c(paste("1",1:4,sep="_"),paste("2",1:6,sep="_"),3:12)) 
table(Idents(dge),id)
table(id)
  1_1   1_2   1_3   1_4   2_1   2_2   2_3   2_4   2_5   2_6     3     4     5 
 2297   400   580   223   341  1182  8207  5707  1022   540  8391 13903  4125 
    6     7     8     9    10    11    12 
 2337  4539  2888   846 10251  2686  7071 
dge$ft124CCA3allFT367assign21<-id
Idents(dge) <- dge$ft124CCA3allFT367assign21
dgeall=dge
dge4=dge

save(dgeall,file="FallopianTube123467_CCA-6subjects-allgenes_11PCs.Robj")


# remove doublet cluster 13
dge=dgeall
table(dgeall$ft124CCA3allFT367assign)
dge=subset(dgeall,ft124CCA3allFT367assign != 13)
dgeNo13=dge

save(dgeNo13,file="FallopianTube123467_CCA-6subjects-allgenes_No13.Robj")


Idents(dge) <- dge$ft124CCA3allFT367assign21
IdStatus<- paste(Idents(dge),dge$Status,sep="_")
names(IdStatus) = names(dge$Status)
levels=names(table(IdStatus))
levels=paste(rep(levels(Idents(dge)),each=2),c("Healthy","Disease"),sep="_")
IdStatus=factor(IdStatus,levels=levels)
table(IdStatus)
dge$IdStatus <- IdStatus
dgeNo13=dge

save(dgeNo13,file="FallopianTube123467_CCA-6subjects-allgenes_No13.Robj")




# C1-2 epithelial subset for Healthy Vs Disease: ciliated 4 subclusters + secretory 6 subclusters subset
ccname="1-2"
dge=subset(dgeall,ft124CCA3allFT367assign %in% c(1:6))
Idents(dge) <- dge$ft124CCA3allFT367assign21
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6])
DefaultAssay(dge) <- "RNA"        

dge12=dge

Idents(dge) <- dge$ft124CCA3allFT367assign21
IdStatus<- paste(Idents(dge),dge$Status,sep="_")
names(IdStatus) = names(dge$Status)
levels=names(table(IdStatus))
levels=paste(rep(levels(Idents(dge)),each=2),c("Healthy","Disease"),sep="_")
IdStatus=factor(IdStatus,levels=levels)
table(IdStatus)
dge$IdStatus <- IdStatus

dge12=dge
save(dge,file=paste0(all,"_C",ccname,".Robj"))

Idents(dge)=dge$IdStatus
centroid12=centroid=log(AverageExpression(dge)$RNA+1)
dge=dge12


ccname="1-2"
load(file=paste0(all,"_C",ccname,".Robj"))
dge12=dge
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6])

dge=dge12
centroid=centroid12


# C1-8 epithelial subclusters + new stromal 
ccname="1-8"
dge=subset(dgeall,ft124CCA3allFT367assign %in% c(1:12))
Idents(dge) <- dge$ft124CCA3allFT367assign21
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[2:7])
DefaultAssay(dge) <- "RNA"        

dge8=dge

Idents(dge) <- dge$ft124CCA3allFT367assign21
IdStatus<- paste(Idents(dge),dge$Status,sep="_")
names(IdStatus) = names(dge$Status)
levels=names(table(IdStatus))
levels=paste(rep(levels(Idents(dge)),each=2),c("Healthy","Disease"),sep="_")
IdStatus=factor(IdStatus,levels=levels)
table(IdStatus)
dge$IdStatus <- IdStatus

Idents(dge) <- dge$ft124CCA3allFT367assign21
IdStatus<- paste(Idents(dge),dge$Status2,sep="_")
names(IdStatus) = names(dge$Status2)
levels=names(table(IdStatus))
levels=paste(rep(levels(Idents(dge)),each=3),c("Healthy","Disease1","Disease2"),sep="_")
IdStatus=factor(IdStatus,levels=levels)
table(IdStatus)
dge$IdStatus2 <- IdStatus


dge8=dge
save(dge,file=paste0(all,"_C",ccname,".Robj"))

Idents(dge)=dge$IdStatus
centroid8=centroid=log(AverageExpression(dge)$RNA+1)
write.table(centroid8,file=paste0("plot/",all,"_C",ccname,"_IdStatus_centroid.txt"),quote=F,row.names=T,col.names=T)
Idents(dge)=dge$IdStatus2
centroid8=centroid=log(AverageExpression(dge)$RNA+1)
write.table(centroid8,file=paste0("plot/",all,"_C",ccname,"_IdStatus2_centroid.txt"),quote=F,row.names=T,col.names=T)
dge=dge8

ccname="1-8"
load(file=paste0(all,"_C",ccname,".Robj"))
dge8=dge
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[2:7])

dge=dge8
centroid8=read.table(file=paste0("plot/",all,"_C",ccname,"_IdStatus_centroid.txt"))
colnames(centroid8)=gsub("X","",colnames(centroid8))
centroid=centroid8
centroid8=read.table(file=paste0("plot/",all,"_C",ccname,"_IdStatus2_centroid.txt"))
colnames(centroid8)=gsub("X","",colnames(centroid8))
centroid=centroid8


# C1 ciliated 4 subclusters for Healthy Vs Disease:
ccname=1
dge=subset(dgeall,ft124CCA3allFT367assign %in% 1)
Idents(dge) <- dge$ft124CCA3allFT367assign21
myBrewerPalette=gg_color_hue(4)
DefaultAssay(dge) <- "RNA"        

dge1=dge

Idents(dge) <- dge$ft124CCA3allFT367assign21
IdStatus<- paste(Idents(dge),dge$Status,sep="_")
names(IdStatus) = names(dge$Status)
levels=names(table(IdStatus))
levels=paste(rep(levels(Idents(dge)),each=2),c("Healthy","Disease"),sep="_")
IdStatus=factor(IdStatus,levels=levels)
table(IdStatus)
dge$IdStatus <- IdStatus

dge1=dge
save(dge,file=paste0(all,"_C",ccname,".Robj"))

ccname=1
load(file=paste0(all,"_C",ccname,".Robj"))
dge1=dge
myBrewerPalette=gg_color_hue(4)


### calculate ECM scores as fraction of ECM genes expressed
# 3 lists of ECM genes
ECM=c("EMILIN1","FN1","LAMA4","LAMA5","LAMB2","LAMC1","NID1","NID2","TGFBI","TNC","VTN","VWA1","LAMA2","FBN1","FBN2","FGA","FGG","LAMB1","TINAGL1","TNXB","AGRN","LAMA3","FGB","IGFBP7","DPT","POSTN","PXDN","FBLN5","EFEMP1","MFAP2","ECM1","LTBP2","ELN","EMILIN2","FBLN2","LTBP1","LTBP4","LAMA1","PAPLN","VWF","AEBP1","HMCN1","THSD4","SRPX","MFGE8","LAMC2","THBS1","MFAP5","VWA5A","MATN2","HMCN2","MFAP4","EMILIN3","LAMB3","SBSPON","LAMB4","EMID1","CILP","MMRN2","FBLN1","NPNT","COLQ","ABI3BP","EFEMP2","PCOLCE","LTBP3","IGFBP3","MXRA5","TINAG","LAMC3","SPON1","SRPX2","MGP","GLDN","FNDC1","SNED1","LGI4","MMRN1","FGL2","THBS2","FRAS1","NTN4","IGFALS","DMBT1","FBN3","OTOG","THBS3","SPARC","WISP2","SVEP1","INTS14","NTN3","FBLN7","MATN4","CILP2","BMPER","IGFBP5","IGFBP4","TNN","CTGF","CYR61","FNDC8","MFAP3","RELN","VWCE","LGI1","MEPE","MFAP1","PXDNL","SPARCL1","IGFBP6","KCP","MATN1","NDNF","NTN1","PCOLCE2","THBS4","CRELD1","DSPP","ECM2","VWA5B1","COMP","SPON2","SPP1","TNFAIP6","VWDE","VWA5A")
# VWA9 should be INTS14
# AW551984 should be VWA5A
ECM[which(!(ECM %in% rownames(dge@assays$RNA@data)))]
#[1] "TINAG"  "OTOG"   "FNDC8"  "VWA5B1"

length(ECM)  # 127
ECM=ECM[which(ECM %in% rownames(dge@assays$RNA@data))]
gene=ECM
length(gene) # 127 


cc1=centroid[gene,]

exp=colSums(as.matrix(expm1(GetAssayData(dge)[gene, ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
#colSums(as.matrix(expm1(GetAssayData(dge)))) # 10000 repeat for each cell
dge$ECM1=exp

expeach=NULL
for(id in levels(Idents(dge))){
tmpcells=names(Idents(dge))[which(Idents(dge)==id)]
print(mean(exp[tmpcells]))
}


COL=c("COL11A2","COL14A1","COL18A1","COL1A1","COL1A2","COL22A1","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3","COL4A5","COL5A1","COL5A2","COL6A2","COL6A3","COL4A4","COL4A6","COL9A3","COL12A1","COL5A3","COL6A1","COL7A1","COL8A1","COL10A1","COL11A1","COL19A1","COL16A1","COL15A1","COL9A1","COL27A1","COL24A1","COL9A2","COL23A1","COL13A1","COL8A2","COL6A5","COL6A6","COL25A1","COL28A1","COL21A1","COL17A1","COL26A1","COL20A1","COL6A4P1")
# COL6A4 should be COL6A4P1
COL[which(!(COL %in% rownames(dge@assays$RNA@data)))]
#[1] "COL6A5" 

length(COL)  # 45
COL=COL[which(COL %in% rownames(dge@assays$RNA@data))]
gene=COL
length(gene) # 44

cc2=centroid[gene,]

exp=colSums(as.matrix(expm1(GetAssayData(dge)[gene, ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
dge$COL=exp

for(id in levels(Idents(dge))){
tmpcells=names(Idents(dge))[which(Idents(dge)==id)]
print(mean(exp[tmpcells]))
}


PRO=c("BGN","HSPG2","LUM","PRELP","ASPN","DCN","VCAN","OGN","PRG4","PRG2","PRG3","FMOD","CHADL","OPTC","IMPG1","PODN","HAPLN1","ACAN","CHAD","SPOCK3")
# VWA9 should be INTS14
# AW551984 should be VWA5A
PRO[which(!(PRO %in% rownames(dge@assays$RNA@data)))]
#[1] "PRG3"

length(PRO)  # 20
PRO=PRO[which(PRO %in% rownames(dge@assays$RNA@data))]
gene=PRO
length(gene) # 19 

cc3=centroid[gene,]

exp=colSums(as.matrix(expm1(GetAssayData(dge)[gene, ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
dge$PRO=exp

for(id in levels(Idents(dge))){
tmpcells=names(Idents(dge))[which(Idents(dge)==id)]
print(mean(exp[tmpcells]))
}


pdf(paste0(dgefile,"ECMscore_Violin.pdf"),height=2,width=5)
VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,pt.size=-1,cols=myBrewerPalette)
dev.off()

png(paste0(dgefile,"ECMscore_Violin_split.png"),res=300,height=1000,width=3000)
  VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,split.by="Status",pt.size=-1)
dev.off()

Idents(dge) <- dge$ft124CCA3allFT367assign21
png(paste0(dgefile,"C1-8_ECMscore_Violin_split.png"),res=300,height=1000,width=4000)
  VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
png(paste0(dgefile,"C1-8_ECMscore_Violin2_split.png"),res=300,height=1000,width=4500)
  VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,split.by="Status2",cols=cols,pt.size=-1)
dev.off()

Idents(dge) <- dge$ft124CCA3allFT367assign21
png(paste0(dgefile,"C1-2_ECMscore_Violin_split.png"),res=300,height=1000,width=3600)
  VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

# C1-2
png(paste0(dgefile,"C",ccname,"_ECMscore_Violin_split.png"),res=300,height=1000,width=3600)
  VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
# C1-8
png(paste0(dgefile,"C",ccname,"_ECMscore_Violin_split.png"),res=300,height=1000,width=4000)
  VlnPlot(dge,c("ECM1","COL","PRO"),ncol=3,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
# saved as Figure 6C. 


# heatmap for each individual ECM gene standardized across healthy and diease
labels=c("ECM","COL","PRO")
geneslist=list(ECM,COL,PRO)

### Genes Standardized Across Cell Types
# note: used this
centroid=centroid[which(apply(centroid,1,sum)!=0),]
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)

### Visualize markers in heatmap across all cell types
data.use=centroid.std

levels=colnames(centroid.std)

colsep.use=cumsum(table(gsub("_Disease.*","",gsub("_Healthy","",levels))))
col.lab=rep("",length(levels))
col.lab=levels

ncluster=length(levels)
sidecol=matrix(0,2,length(levels))

# not missing any column
sidecol[1,]=rep(myBrewerPalette,each=2)
sidecol[2,]=rep(cols,length(levels)/2)

# missed 1_2_Disease column
sidecol[1,]=c(rep(myBrewerPalette[1],2),myBrewerPalette[2],rep(myBrewerPalette[-(1:2)],each=2))
sidecol[2,]=c(cols,cols[1],rep(cols,length(myBrewerPalette)-2))

# missed columns 1_2_Disease1, 1_4_Disease2, 2_4_Disease2
sidecol[1,]=c(rep(myBrewerPalette[1],3),rep(myBrewerPalette[2],2),rep(myBrewerPalette[3],3),rep(myBrewerPalette[4],2),rep(myBrewerPalette[5:7],each=3),rep(myBrewerPalette[8],2),rep(myBrewerPalette[-c(1:8)],each=3))
sidecol[2,]=c(cols,cols[c(1,3)],cols,cols[c(1,2)],rep(cols,3),cols[c(1,2)],rep(cols,length(myBrewerPalette)-8))

clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("Cluster","Status")
colnames(clab)=c("Status","Cluster")

library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
col.use=redblue100

for(g in 1:3){
  genes=geneslist[[g]]
  label=labels[g]
data.use=centroid.std[which(rownames(centroid.std) %in% genes ),]
write.table(centroid[which(rownames(centroid) %in% genes),],paste0(dgefile,"C",ccname,"_centroid_",label,".txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(data.use,paste0(dgefile,"C",ccname,"_centroid_std_",label,".txt"),row.names=T,col.names=T,quote=F,sep="\t")

row.lab=rownames(data.use)
if(g==1){
jpeg(file=paste0(dgefile,"C",ccname,"_centroid_std_",label,".jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="row",Colv=FALSE,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
}
if(g==2){
  jpeg(file=paste0(dgefile,"C",ccname,"_centroid_std_",label,"2.jpeg"),res=300,height=2000,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="row",Colv=FALSE,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.8,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5))
dev.off()
}
if(g==3){
jpeg(file=paste0(dgefile,"C",ccname,"_centroid_std_",label,"3.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="row",Colv=FALSE,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.8,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5))
dev.off()
}

}
# saved as Figure 6D and Table S9.



# visualize FT3,6,7 clusters supervised assigned to FT124 clusters
# use FT124,6 healthy as grey background, color FT3,7 disease clusters in UMAP 
dge=dgeNo13
DefaultAssay(dge) <- "RNA"
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[c(2:7,9:12)])

library(scales)
xlim=ylim=list()
pos=c("topleft","topright","bottomright","bottomright")
xlim=list(range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"umap")[,1]),range(Embeddings(dge,"tsne")[,1]))
ylim=list(range(Embeddings(dge,"pca")[,2]),range(Embeddings(dge,"pca")[,3]),range(Embeddings(dge,"umap")[,2]),range(Embeddings(dge,"tsne")[,2]))
dim=list(Embeddings(dge,"pca")[,1:2],Embeddings(dge,"pca")[,c(1,3)],Embeddings(dge,"umap"),Embeddings(dge,"tsne"))

id1=dge$Status
names(id1)=rownames(dge@meta.data)
id1=sort(id1)
names1=names(id1[which(id1 == "Healthy")])

ident=as.character(dge$ft124CCA3allFT367assign22)
names(ident)=rownames(dge@meta.data)
ident[names1] <- "Healthy"
ident=factor(ident,levels=c("Healthy",levels((dge$ft124CCA3allFT367assign22))))
ident=sort(ident)

j=3
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(myBrewerPalette,0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])

png("plot/FT37assignedclusters_v8UMAP.png",res=300,height=2000,width=2000)
par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(myBrewerPalette,0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
dev.off()
# saved as Figure 6A.

### visualize 4 healthy subjects in UMAP 
dge1=subset(dge,Status=="Healthy")

Idents(dge1) <- dge1$ft124CCA3allFT367assign22
id=as.character(Idents(dge1))
id[grepl("2_",id)] <- 2
id=as.numeric(id)
names(id)=names(Idents(dge1))
id=factor(id,levels=1:12)
Idents(dge1) <- id

myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])
plotlist=list()
png("plot/FT1246_v8UMAP_2.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge1,pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge1,c(1,3),pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge1,reduction="umap",pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge1,pt.size=.8,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()
png("plot/FT1246_v8UMAP.png",res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge1,pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge1,c(1,3),pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge1,reduction="umap",pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge1,pt.size=.8,label=TRUE,label.size=6,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()




### differentially-expressed markers between disease and healthy for each cell type
Idents(dgeall) <- dgeall$ft124CCA3allFT367assign22
dgefile="plot/"


markerslist=plotlist=list()
for(label in levels(dgeall$ft124CCA3allFT367assign22)){
  dge=subset(dgeall,ft124CCA3allFT367assign22 == label)
  Idents(dge)=dge$Status
### visualize cell size factors for each cell type
pdf(file=paste0("plot/C",label,"_PerCellAttributes_ViolinPlot3.pdf"),height=2.5,width=7.5)
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[3]]=VlnPlot(dge, features = "percent.mt",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 3)
dev.off()
### p-value and fold change for markers
DefaultAssay(dge) <- "RNA"        # should use this
gene=FindMarkers(dge,"Disease",only.pos=FALSE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(gene,paste0(dgefile,"C",label,"_RNAassay_Disease_mindiff0.2_logfc2fold_3.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
print(table(Idents(dge)))
print(c(length(which(gene$avg_logFC<0)),length(which(gene$avg_logFC>0))))
gene1=rownames(gene[which(gene$avg_logFC>0),])
gene2=rownames(gene[which(gene$avg_logFC<0),])
genes1=c(gene1,gene2)
gene=FindMarkers(dge,"Disease",only.pos=FALSE,logfc.threshold = log(1.5))
write.table(gene,paste0(dgefile,"C",label,"_RNAassay_Disease_min0.1_logfc1.5fold_3.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
print(c(length(which(gene$avg_logFC<0)),length(which(gene$avg_logFC>0))))
gene1=rownames(gene[which(gene$avg_logFC>0),])
gene2=rownames(gene[which(gene$avg_logFC<0),])
genes2=c(gene1,gene2)
genes2=genes2[which(!(genes2 %in% genes1))]
### p-value and fold change for all genes and volcano plot
DefaultAssay(dge) <- "RNA"        # should use this
markers1=FindMarkers(dge,"Disease",only.pos=FALSE,logfc.threshold = -Inf, min.diff.pct=-Inf, min.pct=-Inf,min.cells.feature=-Inf,min.cells.group=-Inf)
a=markers1
jpeg(file=paste0("plot/HealthyVsDisease_C",label,"_Volcano.jpeg"),res=300,height=1600,width=1600)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(a$avg_logFC,-log10(a$p_val),pch=16,cex=0.8,col=rgb(0,0,0,0.3),xlab="log Fold Change",ylab="-log10(p-value)")
if(length(genes1)!=0){
points(a[genes1,]$avg_logFC,-log10(a[genes1,]$p_val),cex=0.8,pch=16,col=rgb(1,0,0,0.8))
}
if(length(genes2)!=0){
points(a[genes2,]$avg_logFC,-log10(a[genes2,]$p_val),cex=0.8,pch=16,col=rgb(1,0,1,0.8))
}
plot((a$pct.1-a$pct.2)*100,-log10(a$p_val),pch=16,cex=0.8,col=rgb(0,0,0,0.3),xlab="%Cells Diff",ylab="-log10(p-value)")
if(length(genes1)!=0){
points((a[genes1,]$pct.1-a[genes1,]$pct.2)*100,-log10(a[genes1,]$p_val),pch=16,cex=0.8,col=rgb(1,0,0,0.8))
}
if(length(genes2)!=0){
points((a[genes2,]$pct.1-a[genes2,]$pct.2)*100,-log10(a[genes2,]$p_val),pch=16,cex=0.8,col=rgb(1,0,1,0.8))
}
plot(a$avg_logFC,(a$pct.1-a$pct.2)*100,pch=16,cex=0.8,col=rgb(0,0,0,0.3),xlab="log Fold Change",ylab="%Cells Diff")
if(length(genes1)!=0){
points(a[genes1,]$avg_logFC,(a[genes1,]$pct.1-a[genes1,]$pct.2)*100,pch=16,cex=0.8,col=rgb(1,0,0,0.8))
}
if(length(genes2)!=0){
points(a[genes2,]$avg_logFC,(a[genes2,]$pct.1-a[genes2,]$pct.2)*100,pch=16,cex=0.8,col=rgb(1,0,1,0.8))
}
dev.off()
### calculate average expression for each gene
centroid=log(AverageExpression(dge)$RNA+1)
print(which(!(rownames(centroid) %in% rownames(markers1)))) # integer(0)
centroid=data.frame(gene=rownames(centroid),centroid)
markers2=data.frame(gene=rownames(markers1),markers1)
markers=merge(markers2,centroid,by="gene",all=TRUE)
rownames(markers)=markers$gene
markers$diff=markers$pct.1-markers$pct.2
markers1=markers
write.table(markers1,paste0("plot/HealthyVsDisease_C",label,"_allgenes_de_centroid_03.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markerslist[[label]]=markers1
}

### organize markers per cluster in one file
markerslist=list()
for(label in levels(dgeall$ft124CCA3allFT367assign22)){
markers1=read.table(paste0("plot/HealthyVsDisease_C",label,"_allgenes_de_centroid_03.2021.txt"),header=T,stringsAsFactors=F)
markerslist[[label]]=markers1
}

for(label in levels(dgeall$ft124CCA3allFT367assign22)){
print(which(markerslist[[1]][,1] != markerslist[[label]][,1]))
}

for(label in levels(dgeall$ft124CCA3allFT367assign22)){
markerslist[[label]]=markerslist[[label]][,-1]
}
markersall=do.call(cbind,markerslist)
write.table(markersall,paste0("plot/HealthyVsDisease_allclusters_allgenes_de_centroid_03.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")


for(label in levels(dgeall$ft124CCA3allFT367assign22)){
markers1=read.table(paste0("plot/HealthyVsDisease_C",label,"_allgenes_de_centroid_03.2021.txt"),header=T,stringsAsFactors=F)
write.table(markers1[,1:3],paste0("plot/HealthyVsDisease_C",label,"_LRpath_input_03.2021.txt"),col.names=F,row.names=F,quote=F,sep="\t")
}

### organize all LRpath results in one file
levels=c(1,paste0("2-6_",1:6),7:12,14:17)
markerslist=list()
for(label in levels){
markers1=read.table(paste0("plot/C",label,".txt"),quote="\n",fill=TRUE,sep="\t",stringsAsFactors=F)
markers=markers1[,c(1:8)]
markers=markers[which(markers[,1]!="Id"),]
markers=unique(markers)
colnames(markers)=c("Id","Name","Type","nGene",paste0("C",label,c("Coeff","OR","P","FDR")))
markers[,1]=gsub("^O","GO",markers[,1])
#markers$GO=apply(markers[,1:4],1,function(x)paste(x,collapse="-"))
#markers=markers[,c(9,5:8)]
markerslist[[label]]=markers
}

markersall=Reduce(function(x,y) merge(x,y,all=TRUE,by=c("Id","Name","Type","nGene")), markerslist)
unique(markersall[,1])
dim(markersall)
length(unique(markersall[,1])) # 7984

write.table(markersall,paste0("plot/HealthyVsDisease_allclusters_LRpath_03.2021.txt"),col.names=T,row.names=F,quote=F,sep="\t")



# 3.27.21 notes
# 8)  New Stromal (cluster 3-8) violin plot for comparison between disease (FT3) and healthy (FT124) 
knownmarkers=unique(c("TP53","CDKN1A","CDKN2A","COL1A2","ADH1B","COL3A1","THY1","S100A10","LUM","COL6A1","GAS1","PDGFRL","RGCC","CD34","RGS2","CHRDL1","FN1","CITED2","COL15A1","CYP1B1","PLAGL1","COL6A3","MSX1","PIK3R1","HAS1","H19","COL5A2","APOE","TWIST2","IFITM1","SFRP4","ECM1","ALDH1A2","POSTN","DES","ACTG2","MYLK","MYH11","ACTA2","SLMAP","CLDN1","AREG","NOTCH3","CDKN1A","CRIM1","LBH","ALKAL2","PDGFA","ZFHX3","ITGA8","RRAS","RASGRP2","SIVA1","DDIT4"))
#P53: TP53
#P21: CDKN1A
#P16: CDKN2A
length(knownmarkers) # 53
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

dge=subset(dgeall,ft124CCA3allFT367assign22 %in% c(3:8))
Idents(dge) <- dge$ft124CCA3allFT367assign22
myBrewerPalette=brewer.pal(12,"Paired")[2:7]
DefaultAssay(dge) <- "RNA"

dge33=dge

jpeg(paste0(dgefile,"DiseaseVsHealthy_C3-8_Violin_split.jpg"),res=300,height=5000,width=10000)
  VlnPlot(dge,gene,ncol=10,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

jpeg(paste0(dgefile,"DiseaseVsHealthy_C3-8_Violin_split2.jpg"),res=300,height=3800,width=7900)
  VlnPlot(dge,gene,ncol=10,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()


human=read.table("/home/qzm/SCENIC/FT124C3-8/output/regulonAUC.txt",row.names=1,header=T,sep="\t")
tf=rownames(human)
tf=unique(gsub(" .*", "", gsub("_.*","",tf)) )
tf[which(!(tf %in% rownames(dge)))]
length(tf) # 68
gene=tf

jpeg(paste0(dgefile,"DiseaseVsHealthy_C3-8_tf_Violin_split2.jpg"),res=300,height=4600,width=7900)
  VlnPlot(dge,gene,ncol=10,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()


gene="PRRX1"
jpeg(paste0(dgefile,"DiseaseVsHealthy_C3-8_Violin_split_PRRX1.jpg"),res=300,height=660,width=1100)
  VlnPlot(dge,gene,ncol=1,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

gene="MKI67"
jpeg(paste0(dgefile,"DiseaseVsHealthy_C3-8_Violin_split_MKI67.jpg"),res=300,height=660,width=1100)
  VlnPlot(dge,gene,ncol=1,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

jpeg(paste0(dgefile,"DiseaseVsHealthy_C3-8_Violin_split1_MKI67.jpg"),res=300,height=660,width=1100)
  VlnPlot(dge,gene,ncol=1,split.by="Status",pt.size=-1)
dev.off()



# 9)  We will also need violin plots comparing disease and normal for the subclustering of epithelial cells for the markers below
a.  Cluster 1: CAVIN2, VWF, SOX17, DCN, PECAM1, MCAM
b.  Cluster 2: ACTA2, TAGLN, LSGALS1, DES, MYLK, SPARC, S100A4, PDGFRB, PDGFRA, THY1
c.  Cluster 3: JUN, FOS
d.  Cluster 4: LAMB3, PLAUR, CXCL8, CCL20, BRIC3, LIF, MSLN
e.  Cluster 5: CYP1B1, COL1A2, OVGP1, LGR5, POSTN, CD44
f.  Cluster 6: RUNX3, CD3E, TPSAB1, Kit
g.  EMT specific: VIM, TIMP3, SPARC, COL1A1, TWIST1, TWIST2, S100A4, PRRX1, MMP2, ZEB1, ZEB2, FN1, CDH11, CDH2, SNAI1, SNAI2, MUC1, OCLN, LOX, CD44, PROM1 (CD133), POU5F1 (OCT-4),
h.  Hormone receptors: ESR, PGR, AR
i.  Additional Genes of interest for hydrosalpinx disease process: LIF, SUSD2, TPSAB1, MSLN

knownmarkers=c("CAVIN2","VWF","SOX17","DCN","PECAM1","MCAM",
"ACTA2","TAGLN","LGALS1","DES","MYLK","SPARC","S100A4","PDGFRB","PDGFRA","THY1",
"JUN","FOS",
"LAMB3","PLAUR","CXCL8","CCL20","BIRC3","LIF","MSLN",
"CYP1B1","COL1A2","OVGP1","LGR5","POSTN","CD44",
"RUNX3","CD3E","TPSAB1","KIT",
"VIM","TIMP3","SPARC","COL1A1","TWIST1","TWIST2","S100A4","PRRX1","MMP2","ZEB1","ZEB2","FN1","CDH11","CDH2","SNAI1","SNAI2","MUC1","OCLN","LOX",
"CD44","PROM1","POU5F1",
"ESR1","PGR","AR",
"LIF","SUSD2","TPSAB1","MSLN")
length(knownmarkers) # 64
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

dge=subset(dgeall,ft124CCA3allFT367assign %in% c(1:6))
Idents(dge) <- dge$ft124CCA3allFT367assign22
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6])
DefaultAssay(dge) <- "RNA"

dge11=dge

jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_split.jpg"),res=300,height=5800,width=10000)
  VlnPlot(dge,gene,ncol=10,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_split2.jpg"),res=300,height=4600,width=8000)
  VlnPlot(dge,gene,ncol=10,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

# RSPO1 – regulator of wnt signaling (Hu, Cancer Cell)
gene="RSPO1"
jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_splitRSPO1.jpg"),res=300,height=660,width=1200)
  VlnPlot(dge,gene,ncol=1,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

gene="MKI67"
jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_split_MKI67.jpg"),res=300,height=660,width=1200)
  VlnPlot(dge,gene,ncol=1,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_split1_MKI67.jpg"),res=300,height=660,width=1200)
  VlnPlot(dge,gene,ncol=1,split.by="Status",pt.size=-1)
dev.off()

human=read.table("/home/qzm/SCENIC/FT124C2-6/output/regulonAUC.txt",row.names=1,header=T,sep="\t")
tf=rownames(human)
tf=unique(gsub(" .*", "", gsub("_.*","",tf)) )
tf[which(!(tf %in% rownames(dge)))]
length(tf) # 63
gene=tf

jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_tf_Violin_split2.jpg"),res=300,height=4600,width=8000)
  VlnPlot(dge,gene,ncol=10,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()




dge=subset(dgeall,ft124CCA3allFT367assign %in% c(1:6))
Idents(dge) <- dge$ft124CCA3allFT367assign22
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6])
DefaultAssay(dge) <- "RNA"

dge11=dge

knownmarkers=c("DCN","MMP2","PRRX1","SPARC","LGALS1","POSTN","ZEB1","LIF")
length(knownmarkers) # 8
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers
jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_splitFigS6B1.jpg"),res=300,height=3280,width=2000)
  VlnPlot(dge,gene,ncol=2,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
knownmarkers=c("THY1","LGR5","PDGFRA","PDGFRB")
length(knownmarkers) # 4
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers
jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_splitFigS6B2.jpg"),res=300,height=1640,width=2000)
  VlnPlot(dge,gene,ncol=2,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
knownmarkers=c("MUC1","OVGP1","CYP1B1","POU5F1","S100A4")
length(knownmarkers) # 6
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers
jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_splitFigS6C.jpg"),res=300,height=2460,width=2000)
  VlnPlot(dge,gene,ncol=2,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
# saved as FIgure S6B-C.


5.14.2021
# Could you compare (as violin plots) these gene between disease and healthy in the epithelial subclusters 1-1 to 1-4 and 2-1 to 2-6, as well as the global clusters 3-12?
a.    TNFSF13B (BAFF), TNFRSF13B (TACI), TNFRSF17 (BCMA), CXCL8, CCL20, CXCL1
b.       ALDH1A1 and ALDH1A2

knownmarkers=c("TNFSF13B","TNFRSF13B","TNFRSF17","CXCL8","CCL20","CXCL1","ALDH1A1","ALDH1A2")
length(knownmarkers) # 8
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

# Fig6EF
knownmarkers=c("TNF","TNFSF13B","TGFB1","YAP1","CXCL8","TNFRSF17","PDGFB","CCN2")
length(knownmarkers) # 8
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers

5.26.2021
#  TGFB1, PDGFA, PDGFB, IL4, IL13, CCN2 (CTGF), YAP1, TAFAZZIN (TAZ), TNF, IL1B
knownmarkers=c("TGFB1","PDGFA","PDGFB","IL4","IL13","CCN2","YAP1","TAZ","TNF","IL1B")

length(knownmarkers) # 10
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers



# all epithelial subclusters and global clusters
dge=dgeall
dge=subset(dgeall,ft124CCA3allFT367assign != 13)
dgeNo13=dge

dge=dgeNo13
DefaultAssay(dge) <- "RNA"
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[c(2:7,9:12)])

# added an additional cell for the missing category before running below
jpeg(paste0(dgefile,"DiseaseVsHealthy_All_Violin_split514.jpg"),res=300,height=1400,width=6250)
  VlnPlot(dge,gene,ncol=4,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()
jpeg(paste0(dgefile,"DiseaseVsHealthy_All_Violin_split526.jpg"),res=300,height=2100,width=6250)
  VlnPlot(dge,gene,ncol=4,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

### adjust Violin plot to be wider for those only showing ticks
cols=c("#F8766D", "#078992") 
data=data.frame(t(dge@assays$RNA@data[knownmarkers,]),CellType=dge$ft124CCA3allFT367assign21,Status=dge@meta.data$Status)
for(GENE in knownmarkers){
p1=ggplot(data,aes(x=CellType,y=get(GENE)))+geom_violin(aes(fill=Status),scale="width",adjust=2)+scale_fill_manual(values = c(cols,cols[1],rep(cols,12)))+ theme(legend.position = "none")+
theme_bw()+
theme(axis.text.x = element_text(angle = 45,hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_blank(),
       axis.line    = element_line(color='black'))+
ggtitle(GENE)+
ylab("Expression Level")
jpeg(paste0(dgefile,paste0("DiseaseVsHealthy_All_Violin_split_",GENE,"_v2.jpg")),res=300,height=660,width=2500)
print(p1)
dev.off()
}




ccname="3-8"
load(file=paste0(all,"_C",ccname,".Robj"))
dge33=dge
myBrewerPalette=brewer.pal(12,"Paired")[2:7]

jpeg(paste0(dgefile,"DiseaseVsHealthy_Stromal_Violin_split514.jpg"),res=300,height=1200,width=2400)
  VlnPlot(dge,gene,ncol=4,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()


ccname="1-2"
load(file=paste0(all,"_C",ccname,".Robj"))
dge12=dge
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6])

jpeg(paste0(dgefile,"DiseaseVsHealthy_Epithelial_Violin_split514.jpg"),res=300,height=1400,width=3200)
  VlnPlot(dge,gene,ncol=4,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

ccname="1-8"
load(file=paste0(all,"_C",ccname,".Robj"))
dge8=dge
myBrewerPalette=c(gg_color_hue(4),brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[2:7])

jpeg(paste0(dgefile,"DiseaseVsHealthy_C",ccname,"_Violin_split514.jpg"),res=300,height=1400,width=5000)
  VlnPlot(dge,gene,ncol=4,split.by="Status",split.plot=TRUE,pt.size=-1)
dev.off()

# 5.27.2021
ciliated markers in ciliated cells in Healthy Vs Disease


knownmarkers=c("RSPH1","RSPH9","RSPH4A","DNAI1","DNAI2","DNAH11","CCDC103","LRRC6","ZMYND10","FOXJ1"
  ,"MUC1")
length(knownmarkers) # 11
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
gene=knownmarkers


ccname=1
load(file=paste0(all,"_C",ccname,".Robj"))
dge1=dge
myBrewerPalette=gg_color_hue(4)

### violin plot
DefaultAssay(dge) <- "RNA"

jpeg(paste0(dgefile,"DiseaseVsHealthy_Ciliated_Violin_split.jpg"),res=300,height=2100,width=2100)
  VlnPlot(dge,gene,ncol=4,split.by="Status",split.plot=TRUE,pt.size=-1)
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
save(dgeall,file=paste0(all,"_CCA-6subjects-allgenes.Robj"))


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
