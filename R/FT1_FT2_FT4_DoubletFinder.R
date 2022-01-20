# Apply DoubletFinder to each individual FT by Qianyi on 10.4.2021

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
myBrewerPalette1=gg_color_hue(4)
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

dgefile="plot/"

# load individual DoubletFinder Seurat object
seulist=list()


# global clusters and subclusters for healthy FTs
id=read.table(file=paste0("plot/",all,"_CCA3-allgenes_ident.txt"),row.names=1)
table(id$V2)

id21=read.table(file=paste0("plot/",all,"_CCA3-allgenes_ident21.txt"),row.names=1)
table(id21$V2)

# load individual sample before aggregation
dgealllist=list()

indiv=1
indiv=3
indiv=5

for(indiv in indivft){

# merged segments, individual subject
  load(file=paste0(subjectorgans[indiv],".Robj")) 


dgealllist[[indiv]]=dgeall
}




library(DoubletFinder)

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu=dgeall
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:10)


save(seu,file=paste0(organs[indiv],"_seu.Robj"))

load(file=paste0(organs[indiv],"_seu.Robj"))

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE) # SCT-transform in Seurat
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

png(file=paste0(dgefile,organs[indiv],"_pK.png"),res=300,height=900,width=1400)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
bcmvn <- find.pK(sweep.stats)
dev.off()
pK=as.numeric(as.character(bcmvn[which(bcmvn$BCmetric==max(bcmvn$BCmetric)),2]))
print(pK)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
seu$id <- id[rownames(seu@meta.data),]
Idents(seu) <- seu$id

seu$id21 <- id21[rownames(seu@meta.data),]
Idents(seu) <- seu$id21

annotations=as.character(seu$id)  # used this

### Assuming 7.5% doublet rate
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pANN=grep("DF",colnames(seu@meta.data),value=T)[1]
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)
pANN2=grep("DF",colnames(seu@meta.data),value=T)[2]

save(seu,file=paste0(organs[indiv],"_seu.Robj"))


[1] "Creating 3509 artificial doublets..."
[1] "Creating 8465 artificial doublets..."
[1] "Creating 5818 artificial doublets..."


table(seu[[pANN]])
table(seu[[pANN2]])
table(seu$id,seu[[pANN]][,1])



# visualize Doublet Vs Singlet in Cluster 1
which(colnames(dgeall) != colnames(seu))
#integer(0)
Idents(seu) <- seu$DF.classifications_0.25_0.005_790
PCAPlot(seu)
### Visualize individual batch and subject
label="DF.classifications_0.25_0.005_790"
dgeall$DF.classifications_0.25_0.005_790 <- seu$DF.classifications_0.25_0.005_790
dgeall$id <- id[rownames(dgeall@meta.data),]
Idents(dgeall) <- seu$DF.classifications_0.25_0.005_790
dge=subset(dgeall,id == 1)
pdf(paste0("plot/",organs[indiv],"_PC1-2.pdf"),height=4,width=4)
par(mfrow=c(1,2),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
PCAPlot(dgeall,group="id")
dev.off()
pdf(paste0("plot/",organs[indiv],"_doublet790_PC1-2.pdf"),height=3,width=4)
par(mfrow=c(1,2),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
PCAPlot(dgeall)
PCAPlot(dge)
dev.off()

### Assuming 0% doublet rate
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu@meta.data$ClusteringResults
nExp_poi <- round(0*nrow(seu@meta.data))  ## Assuming 0% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi
# [1] 0
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pANN="DF.classifications_0.25_0.005_0"
table(seu$id21,seu[[pANN]][,1])
table(seu$id,seu[[pANN]][,1])
         
             1    2    3    4    5    6    7    8    9   10   11   12   13   14
  Doublet    0    0    0    1    0    0    0    0    0    0    0    0    0    0
  Singlet  112    6 1907  166    4   29  836 2820  742  522  750  216   70   73
         
            15   16   17
  Doublet    0    0    0
  Singlet 1910   81  282


### Assuming 1% doublet rate
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu@meta.data$ClusteringResults
nExp_poi <- round(0.01*nrow(seu@meta.data))  ## Assuming 1% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi
# [1] 105  254  175
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pANN="DF.classifications_0.25_0.005_105" # FT1
pANN="DF.classifications_0.25_0.005_254" # FT2
pANN="DF.classifications_0.25_0.005_175" # FT4

table(seu$id,seu[[pANN]][,1])
table(seu$id21,seu[[pANN]][,1])

             1    2    3    4    5    6    7    8    9   10   11   12   13   14
  Doublet    0    0    0   53    0    0    4    0    5    6    2    0   34    0
  Singlet  112    6 1907  114    4   29  832 2820  737  516  748  216   36   73
         
            15   16   17
  Doublet    0    0    1
  Singlet 1910   81  281


### Assuming 2% doublet rate
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu@meta.data$ClusteringResults
nExp_poi <- round(0.02*nrow(seu@meta.data))  ## Assuming 2% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi
# [1] 211 508 349
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pANN="DF.classifications_0.25_0.005_211"
pANN="DF.classifications_0.25_0.005_508"
pANN="DF.classifications_0.25_0.005_349"
table(seu$id,seu[[pANN]][,1])
table(seu$id21,seu[[pANN]][,1])


table(seu[[pANN]][,1],seu$id)
             1    2    3    4    5    6    7    8    9   10   11   12   13   14
  Doublet    0    0    0   80    1   17    7    9   11   11    6    1   66    0
  Singlet  112    6 1907   87    3   12  829 2811  731  511  744  215    4   73
         
            15   16   17
  Doublet    1    0    1
  Singlet 1909   81  281
save(seu,file=paste0(organs[indiv],"_seu.Robj"))
seulist[[indiv]]=seu
}

# save the singlets vs doublets identified by DoubletFinder
pANNs=list(c("DF.classifications_0.25_0.005_105","DF.classifications_0.25_0.005_254","DF.classifications_0.25_0.005_175"),
c("DF.classifications_0.25_0.005_211","DF.classifications_0.25_0.005_508","DF.classifications_0.25_0.005_349") )
for(j in 1:2){
celld=NULL;i=1
for(indiv in indivft){
    seu=seulist[[indiv]]
    pANN=pANNs[[j]][i]
    i=i+1
    tmp=seu[[pANN]]
names(tmp)=paste0(j,"pct")
celld=rbind(celld,tmp)
}
table(gsub("_.*","",celld))

write.table(celld,file=paste0("plot/FT124_DoubletFinder_",j,"pct.txt"),row.names=T,col.names=T,quote=F,sep="\t")
}


# save the FT124 doublets identified by DoubletFinder assuming 1% doublet rate
celld=NULL;i=1
pANNs=c("DF.classifications_0.25_0.005_105","DF.classifications_0.25_0.005_254","DF.classifications_0.25_0.005_175")
for(indiv in indivft){
    seu=seulist[[indiv]]
    pANN=pANNs[i]
    i=i+1
celld=c(celld,rownames(seu[[pANN]])[which(seu[[pANN]][,1]=="Doublet")])
}
table(gsub("_.*","",celld))

Ampulla1 Ampulla2 Ampulla4 Fimbria1 Fimbria2 Fimbria4 Isthmus1 Isthmus2 
      62       96        8       10       51      106       33      107 
Isthmus4 
      61 

write.table(celld,file=paste0("plot/FT124_DoubletFinder_1pct_doublets.txt"),row.names=T,col.names=F,quote=F,sep="\t")

# save the FT124 doublets identified by DoubletFinder assuming 2% doublet rate
celld=NULL;i=1
pANNs=c("DF.classifications_0.25_0.005_211","DF.classifications_0.25_0.005_508","DF.classifications_0.25_0.005_349")
for(indiv in indivft){
    seu=seulist[[indiv]]
    pANN=pANNs[i]
    i=i+1
celld=c(celld,rownames(seu[[pANN]])[which(seu[[pANN]][,1]=="Doublet")])
}
table(gsub("_.*","",celld))

write.table(celld,file=paste0("plot/FT124_DoubletFinder_2pct_doublets.txt"),row.names=T,col.names=F,quote=F,sep="\t")

pdf(paste0("plot/",organs[indiv],"_PC1-2.pdf"),height=4,width=4)
par(mfrow=c(1,2),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
PCAPlot(seu,group="id")
dev.off()
pdf(paste0("plot/",organs[indiv],"_doublet_PC1-2.pdf"),height=3,width=4)
par(mfrow=c(1,2),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
for(nd in c(0,105,211,790)){
plot<-PCAPlot(seu,group=paste0("DF.classifications_0.25_0.005_",nd))
print(plot)
}
dev.off()

dge=dgeall=seu
sets=c("Doublet","Singlet")
plotlist=list()
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlim=ylim=list()
pos=c("bottomright","bottomright","topleft","bottomright")
# dge4 used top 9 PCs for tSNE
xlim=list(range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"umap")[,1]))
ylim=list(range(Embeddings(dge,"pca")[,2]),range(Embeddings(dge,"pca")[,3]),range(Embeddings(dge,"umap")[,2]))
### plot PCs and tSNE for each batch using the other batches as background
dge=dgeall
if(length(sets)>1){
  cols=gg_color_hue(length(sets))
plot2set=plot3set=plot4set=plottset=NULL
dim=list(Embeddings(dge,"pca")[,1:2],Embeddings(dge,"pca")[,c(1,3)],Embeddings(dge,"umap"))
size=length(sets)
if (size>4) {
  pdf(paste0("plot/",organs[indiv],"_doublet_indiv_pN",pN,".pdf"),height=2.3*round(sqrt(size)),width=2.3*ceiling(sqrt(size)))
  par(mfrow=c(round(sqrt(size)),ceiling(sqrt(size))),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
} else {
  pdf(paste0("plot/",organs[indiv],"_doublet_indiv_pN",pN,".pdf"),height=2.5,width=2.5*size)
  par(mfrow=c(1,size),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
}
for(nd in c(0,105,211,790)){
    label=paste0("DF.classifications_0.25_0.005_",nd)
Idents(dgeall)<-dgeall[[label]]
for(j in 1){
#  j=1
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(Idents(dgeall))
names(ident)=names(Idents(dgeall))
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.6,col=c("grey70",alpha(cols[seti],0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
if(label != "label"){
if(size<(round(sqrt(size))*ceiling(sqrt(size)))){
  for(new in (size+1):(round(sqrt(size))*ceiling(sqrt(size)))){
    plot.new()
  }
}
}
}
}
dev.off()
}
}
}
