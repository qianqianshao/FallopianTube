# 3.12.2021 convert to scVelo input files for 3 healthy FTs by Qianyi


R
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(RColorBrewer)

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3")
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium1","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3","Fimbria4","Ampulla4","Isthmus4","Uterus2","FT5","Ovary1")
loomnames=paste0("Sample_",dataset)
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus","FallopianTube2","FallopianTube2","FallopianTube2","FallopianTube3","FallopianTube3","FallopianTube4","FallopianTube4","FallopianTube4","Uterus2","FallopianTube5","Ovary1")
subject=c("Human1","Human1","Human1","Human2","Human3","Human3","Human3","Human4","Human4","Human5","Human5","Human5","Human6","Human6","Human6")
subjectorgans=unique(paste(subject,organ,sep="-"))
subjectorgans=c("FT1","FT2","FT4")

indir="/nfs/turbo/umms-hammou/10x/"
file="FT124"
sets=c("1","_C2-6","_C7-10","_C1-6","_C1-8","_C3-8","_C1-2_2","_C127","_C12678")
indiv=1
all="FallopianTube124"
ss=all[indiv]
filename<-paste0(indir,file[indiv],".loom")


ldat <- ReadVelocity(file = filename)
bm <- as.Seurat(x = ldat)
bm[["RNA"]] <- bm[["spliced"]]

oldnames <- Cells(bm)
newnames=gsub(":","_",oldnames)
newnames=gsub("x","",newnames)

  for(ff in 1:length(names)){
   newnames=gsub(loomnames[ff],names[ff],newnames)
  }

bm=RenameCells(bm, new.names=newnames)
bmall = bm


### v4. CCA-3subjects-allgenes to integrate 3 large healthy FTs (FT1, FT2, FT4)
# remove doublet cluster 13
load(file=paste0(ss,"_CCA-3subjects-allgenes-NoC13.Robj")) # used this
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])
myBrewerPalette2=c(brewer.pal(12,"Paired")[1],brewer.pal(8,"Set2")[1:6],brewer.pal(12,"Paired")[c(2:7,9:12)])
dge=dgeall
write.csv(dgeall@reductions$pca@cell.embeddings[,1:2],"plot/pcaFT124.csv")

cells=colnames(dgeall)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dgeall@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dgeall@reductions$umap@cell.embeddings
bm$ft124CCA3all1 <- dgeall$ft124CCA3all1
bm$ft124CCA3all2 <- dgeall$ft124CCA3all2
bm$seurat_clusters <- dgeall$ft124CCA3all1
Idents(bm) <- dgeall$ft124CCA3all1
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124.h5Seurat")
Convert("FT124.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clusters.csv")
table(Idents(bm))
    1   2-6     7     8     9    10    11    12    14    15    16    17 
 2510 14002  5327  9250  3642  2043  3785  2352   474  5770  1713  2310 
length(table(Idents(bm)))
#[1] 12
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])
myBrewerPalette[1:12]
 [1] "#A6CEE3" "grey"    "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C"
 [8] "#FDBF6F" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

set=sets[2]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC2-6.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunTSNE(bm, dims = 1:20)
bm <- RunUMAP(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3allC2to6 <- dge$ft124CCA3allC2to6
bm$seurat_clusters <- dge$ft124CCA3allC2to6
Idents(bm) <- dge$ft124CCA3allC2to6
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C2-6.h5Seurat")
Convert("FT124_C2-6.h5Seurat", dest = "h5ad")

set=sets[3]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC7-10.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3all <- dge$ft124CCA3all # keep global clusters
bm$seurat_clusters <- dge$ft124CCA3all
Idents(bm) <- dge$ft124CCA3all
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C7-10.h5Seurat")
Convert("FT124_C7-10.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC7-10.csv")
table(Idents(bm))
   7     8     9    10 
5327  9250  3642  2043  
length(table(Idents(bm)))
#[1] 4
myBrewerPalette=brewer.pal(12,"Paired")[2:5]
myBrewerPalette
 [1] "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" 

set=sets[4]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC1-6.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3all2 <- dge$ft124CCA3all2
bm$seurat_clusters <- dge$ft124CCA3all2
Idents(bm) <- dge$ft124CCA3all2
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C1-6.h5Seurat")
Convert("FT124_C1-6.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC1-6.csv")
table(Idents(bm))
    1 2-6_1 2-6_2 2-6_3 2-6_4 2-6_5 2-6_6 
 2510   268   868  6522  5479   549   316  
length(table(Idents(bm)))
#[1] 4
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2"))
myBrewerPalette
1] "#A6CEE3" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"



set=sets[5]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC1-8.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3all22 <- dge$ft124CCA3all22
bm$ft124CCA3all11 <- dge$ft124CCA3all11
bm$seurat_clusters <- dge$ft124CCA3all22
Idents(bm) <- dge$ft124CCA3all22
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C1-8.h5Seurat")
Convert("FT124_C1-8.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC1-8.csv")
table(Idents(bm))
   1  2_1  2_2  2_3  2_4  2_5  2_6    3    4    5    6    7    8 
2510  268  868 6522 5479  549  316 5327 9250 3642 2043 3785 2352 
length(table(Idents(bm)))
#[1] 13
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[c(2:7,9:12)])
myBrewerPalette
 [1] "#A6CEE3" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"
 [8] "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#CAB2D6"
[15] "#6A3D9A" "#FFFF99" "#B15928"

set=sets[6]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC3-8.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3all22 <- dge$ft124CCA3all22
bm$ft124CCA3all11 <- dge$ft124CCA3all11
bm$seurat_clusters <- dge$ft124CCA3all22
Idents(bm) <- dge$ft124CCA3all22
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C3-8.h5Seurat")
Convert("FT124_C3-8.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC3-8.csv")
table(Idents(bm))
length(table(Idents(bm)))
#[1] 4
myBrewerPalette=brewer.pal(12,"Paired")[2:7]
myBrewerPalette
1] "#A6CEE3" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"


set=sets[7]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC1-2_2.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ordered <- dge$ordered
bm$seurat_clusters <- dge$ordered
Idents(bm) <- dge$ordered
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C1-2_2.h5Seurat")
Convert("FT124_C1-2_2.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC1-2_2.csv")
table(Idents(bm))
 1_1  1_2  1_3  1_4  2_2 
1466  387  513  144  868 
length(table(Idents(bm)))
myBrewerPalette1=c(gg_color_hue(4),brewer.pal(7,"Set2")[2])
myBrewerPalette1
[1]  "#F8766D" "#7CAE00" "#00BFC4" "#C77CFF" "#FC8D62"


set=sets[8]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC127.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3all22 <- dge$ft124CCA3all22
bm$seurat_clusters <- dge$ft124CCA3all22
Idents(bm) <- dge$ft124CCA3all22
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C127.h5Seurat")
Convert("FT124_C127.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC127.csv")
table(Idents(bm))
    1 2-6_1 2-6_2 2-6_3 2-6_4 2-6_5 2-6_6 
 2510   268   868  6522  5479   549   316  
length(table(Idents(bm)))
#[1] 8
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[6])
myBrewerPalette
[1] "#A6CEE3" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"
[8] "#E31A1C"


set=sets[9]
load(file=paste0(ss,set,".Robj"))
write.csv(dge@reductions$pca@cell.embeddings[,1:2],"plot/pcaC12678.csv")
cells=colnames(dge)
bm = subset(bmall,cells=cells)
bm<-FindVariableFeatures(bm)
bm <- ScaleData(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
bm@reductions$pca@cell.embeddings <- dge@reductions$pca@cell.embeddings
bm@reductions$umap@cell.embeddings <- dge@reductions$umap@cell.embeddings
bm@reductions$tsne@cell.embeddings <- dge@reductions$tsne@cell.embeddings
bm$ft124CCA3all22 <- dge$ft124CCA3all22
bm$seurat_clusters <- dge$ft124CCA3all22
Idents(bm) <- dge$ft124CCA3all22
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "FT124_C12678.h5Seurat")
Convert("FT124_C12678.h5Seurat", dest = "h5ad")

tmp=data.frame(CellID=rownames(bm@meta.data),Seurat_clusters=bm$seurat_clusters)
write.csv(tmp,"plot/clustersC12678.csv")
table(Idents(bm))
   1  2_1  2_2  2_3  2_4  2_5  2_6    6    7    8 
2510  268  868 6522 5479  549  316 2043 3785 2352 

length(table(Idents(bm)))
#[1] 10
myBrewerPalette=c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[c(5:7)])
myBrewerPalette
[1] "#A6CEE3" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"
[8] "#FB9A99" "#E31A1C" "#FDBF6F"


