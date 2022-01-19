# 11.24.2021 by Qianyi
# visualized subclusters for 4 healthy FTs
# Related to Figure 2-4, S2, S3 and S5. 

R


library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
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


subsets=list(1,2:6,3:8)
subsetsname=c("1","2-6","3-8")
subsetslabel=c("ciliated","NCE","stromal")
### using FT124 original 17 clusters
# 1    ciliated cluster   -> 4 subclusters
# 2:6  secretory subset   -> 6 subclusters
### after renaming the 12 clusters to 1:12
# 3:8   new stromal including myofibroblast and endothelial




# color each subject in UMAP
cols=list(
  gg_color_hue(4),
brewer.pal(7,"Set2")[1:6],
  brewer.pal(12,"Paired")[-1]
)

# C3-8: new stromal subset
myBrewerPalette=brewer.pal(12,"Paired")[2:7]
cols=list( brewer.pal(12,"Paired")[2:7] )


### CCA-4subjects-allgenes to integrate 4 large healthy FTs (FT1, FT2, FT4, FT6)
load(file=paste0(all,"_CCA-4subjects-allgenes.Robj")) # used this
dgeall


### add identity of renamed 12 clusters + epithelial subclusters
# load All 6 subjects - No Cluster 13
load(file="FallopianTube123467_CCA-6subjects-allgenes_No13.Robj")
dgeNo13
id=dgeNo13$ft124CCA3allFT367assign21
cells=rownames(dgeall@meta.data)
dgeall$ft124CCA3allFT6assign21=id[cells]
Idents(dgeall) <- dgeall$ft124CCA3allFT6assign21
save(dgeall,file=paste0(all,"_CCA-4subjects-allgenes.Robj")) # used this




for(i in 1:length(subsets)){

i=1
i=2
i=3


cc=subsets[[i]]
ccname=subsetsname[i]
dgefile=paste0("plot/C",ccname,"_")

myBrewerPalette=cols[[i]]

# Use HVG of FT124 subset
load(file=paste0("FallopianTube124_C",ccname,".Robj"))
dge1=dge
aa=VariableFeatures(dge1)

# keep FT1246 cells in the subset
dge=subset(dgeall,ft124CCA3allFT6assign %in% cc) # original 17 clusters
dge=subset(dgeall,ft124CCA3allFT6assign21 %in% cc) # renamed 12 clusters + epithelial subclusters
table(Idents(dgeall))
table(Idents(dge))
table(dgeall$Organ)
table(dge$Organ)

VariableFeatures(dge)=aa

### Scale data
DefaultAssay(dge)
dge <- ScaleData(dge,features=rownames(dge))
### PCA
Sys.time()  
dge <- RunPCA(dge, features = VariableFeatures(dge),npcs = 50,  ndims.print = 5, nfeatures.print = 5)
dge <- ProjectDim(dge)

###### Determine top PCs
numPCs=c(10,9,8) # original 17 clusters
numPCs=9 # C3-8 new stromal
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
UMAPPlot(dge,group="ft124CCA3allFT6assign21",cols=myBrewerPalette)
UMAPPlot(dge,group="Organ")

save(dge,file=paste0(all,"_FT124C",ccname,".Robj"))


for(i in 1:3){
ccname=subsetsname[i]
dge=dgeall=dgelist[[i]]
for(label in c("Organ")){
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
  pdf(paste0("plot/",all,"_FT124C",ccname,"_",label,"_indiv.pdf"),height=2.3*round(sqrt(size)),width=2.3*ceiling(sqrt(size)))
  par(mfrow=c(round(sqrt(size)),ceiling(sqrt(size))),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
} else {
  pdf(paste0("plot/",all,"_FT124C",ccname,"_",label,"_indiv.pdf"),height=2.5,width=2.5*size)
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
# saved as Figure S2A, S3A, and S5A.  


### add identity of renamed 12 clusters
# All 6 subjects - No Cluster 13
load(file="FallopianTube123467_CCA-6subjects-allgenes_No13.Robj")
dgeNo13

id=dgeNo13$ft124CCA3allFT367assign21
for(i in 1:3){
ccname=subsetsname[i]
dge=dgelist[[i]]
cells=rownames(dge@meta.data)
dge$ft124CCA3allFT6assign21=id[cells]
Idents(dge) <- dge$ft124CCA3allFT6assign21
dgelist[[i]]=dge
}

for(i in 1:3){
ccname=subsetsname[i]
dge=dgelist[[i]]

print(dge)
print(table(Idents(dge)))
save(dge,file=paste0(all,"_FT124C",ccname,".Robj"))
}


### rename 6 FTs to new names
for(i in 1:3){
ccname=subsetsname[i]
dge=dgeall=dgelist[[i]]
# rename donors
id=as.character(dge$Organ)
names(id)=names(Idents(dge))
table(id)
id[which(id == "FallopianTube")]<-"FT1"
id[which(id == "FallopianTube2")]<-"FT2"
id[which(id == "FallopianTube4")]<-"FT3"
id[which(id == "FT6")]<-"FT4"
id=factor(id,levels=paste0("FT",1:4)) 
dge$Sample<-id
# rename dataset names
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
dge$Dataset<-id
# rename cell names
id <- as.character(dge$Dataset)
cells <- rownames(dge@meta.data)
table(gsub("_.*","",cells))
table(id)
new <- paste0(id,gsub(".*_","_",cells))
table(gsub("_.*","",new))
dge=RenameCells(dge,new.names=new)
dgelist[[i]]=dge
save(dge,file=paste0(all,"_FT124C",ccname,"_Renamed.Robj"))
}


dgelist=list()
for(i in 1:length(subsets)){
cc=subsets[[i]]
ccname=subsetsname[i]
load(file=paste0(all,"_FT124C",ccname,"_Renamed.Robj"))
dgelist[[i]]=dge
}


for(i in 1:length(subsets)){
cc=subsets[[i]]
ccname=subsetsname[i]
cclabel=subsetslabel[i]
dge=dgelist[[i]]
cbind(1:ncol(dge@meta.data),colnames(dge@meta.data))
tmp=dge@meta.data[,c(1:4,82:83,81)]
colnames(tmp)[ncol(tmp)] <- "CellType"
dge@meta.data <- tmp
dgelist[[i]]=dge
saveRDS(dge,file=paste0("FallopianTube1234_",cclabel,".rds"))
}



# color each subject in UMAP
cols=list(
  gg_color_hue(4),
brewer.pal(7,"Set2")[1:6],
  brewer.pal(12,"Paired")[-1]
)

for(i in 1:3){
ccname=subsetsname[i]
dge=dgelist[[i]]
Idents(dge) <- dge$ft124CCA3allFT6assign21
plotlist=list()
myBrewerPalette=cols[[i]]
png(paste0("plot/",all,"_FT124C",ccname,"_cluster1.png"),res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.8,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.8,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()
}

# C3-8: new stromal subset
myBrewerPalette=brewer.pal(12,"Paired")[2:7]
cols=list( brewer.pal(12,"Paired")[2:7] )
myBrewerPalette=cols[[i]]
Idents(dge) <- dge$ft124CCA3allFT6assign21
png(paste0("plot/",all,"_FT124C",ccname,"_cluster1.png"),res=300,height=2000,width=2330)
plotlist[[1]]=PCAPlot(dge,pt.size=.4,label=FALSE,cols=myBrewerPalette)
plotlist[[2]]=PCAPlot(dge,c(1,3),pt.size=.4,label=FALSE,cols=myBrewerPalette)
plotlist[[3]]=DimPlot(dge,reduction="umap",pt.size=.4,label=FALSE,cols=myBrewerPalette)
plotlist[[4]]=TSNEPlot(dge,pt.size=.4,label=FALSE,cols=myBrewerPalette)
multiplot(plotlist,cols = 2)
dev.off()

# saved as Figure 2A, 3A, and 4A.


