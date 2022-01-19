###### Correlation between our fallopian tubes (healthy FT1,2,4) and cell type centroids for FT from Huy Dinh 2021 paper
# by Qianyi on 7/20/2021
# Related to Figure S7B.

### Integrated  3 large healthy FTs by Qianyi in Dec 2020
# 10x
12 global cell types, 10 epithelial subtypes (4 ciliated subtypes, 6 secretory subtypes)


### 10 clusters for FT from Huy Dinh 2021 paper
# centroids of 10 clusters:
secretory 1-3, ciliated 1-4, and transitional 1-3 (1 called early secretory)
avg_expr_epithelials_July-2021.csv
data format: log(mean(counts-per-10k)+1), 52 markers, 5 cluster centroids
# markers for 10 clusters:
secretory 1-3, ciliated 1-4, and transitional 1-3 (called unclassified)
mmc4.xlsx
Table S3



###### load our epithelial subset of healthy FT data analyzed by Qianyi in 2021
library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

dgefile="plot/CCA3allNoC13_"
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])


load(file=paste0("FallopianTube124_CCA-3subjects-allgenes-NoC13.Robj")) # used this
dge4=dgeall
avg=AverageExpression(dge4)
ft=log(avg$RNA+1)

# load centroid of our FT
ft=read.table(paste0(dgefile,"21_Centroid-uncorrected.txt"))

### load markers for 12 global cell types of our FT
markers1=read.table(paste0(dgefile,"RNAassay1_mindiff0.2_logfc2fold_2.2021.txt"),stringsAsFactors=F)


### load 10 epithelial cluster centroids for published FT from Dinh
# data format: log(mean(counts-per-10k)+1), 52 genes, 5 cluster centroids (columns)
dinh=read.csv("avg_expr_epithelials_July-2021.csv",row.names=1)
#ft=read.table("ftGonad4ClusterCentroids_RPKM.txt",header=T)
dim(dinh) #  45037    10

### load markers for 10 epithelial clusters of published FT from Dinh
dinhmarkers=read.table("TableS3.txt",header=T,stringsAsFactors=F,sep="\t")
#ft=read.table("ftGonad4ClusterCentroids_RPKM.txt",header=T)
dim(dinhmarkers) # 1467     7
### number of markers for each cluster
table(dinhmarkers$Cluster)
Ciliated_1      Ciliated_2      Ciliated_3      Ciliated_4     Secretory_1 
             88             108             139             100             129 
    Secretory_2     Secretory_3 Unclassified _1  Unclassified_2  Unclassified_3 
            107             137             186             193             280 



### Cross-tabulate using rank correlation 
genesall=rownames(dinh)[which(rownames(dinh) %in% rownames(ft))]
length(genesall) # 27816

  all=cbind(ft[genesall,],dinh[genesall,])

# using HVG from our epithelial subset
hvg=VariableFeatures(dge)
hvg=hvg[which(hvg %in% genesall)]
  rho=cor(all[hvg,],method="spearman")
rho[1:20,21:30]



# using markers from our paper
markers=markers1$gene
markers=markers[which(markers %in% genesall)]
  rho=cor(all[markers,],method="spearman")
rho[1:20,21:30]


######## Heatmap for all markers
dge=dgeall
dgefile="plot/CCA3allNoC13_"
redblue100<-rgb(read.table(paste0(home,'data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}


centroid=dinh[genesall,]
colnames(centroid)=gsub("X","",colnames(centroid))
colnames(centroid)=gsub("_Epithelials","",colnames(centroid))

### Genes Standardized Across Cell Types
# note: used this
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
centroid.std=centroid.std[,c(4,3,10,8,2,1,5,7,6,9)]


### Visualize markers in heatmap across all cell types
# global markers
dgefile="plot/OurGlobal_"
markers=markers1

# C2-6 markers # Union of Global 1VsAllOthers: 1.6FC, 10%pct + Local 3,4,5VsNeighbor
dgefile=paste0("plot/C2-6_")
numPCs=c(6,6,6);res=c(paste0("integrated_snn_res.0.",c(1,1,1),"ordered"));resi=i=1
markers2=read.table(paste0(dgefile,numPCs[i],"PCs_",res[i],"_logfc1.6fold_min0.1_RNAassay_1.2021.txt"),header=T,row.names=1,stringsAsFactors=F)
markers=markers2
markers$cluster=paste0("2_",markers$cluster)

### markers for 10 epithelial clusters of published FT from Dinh
dinhmarkers=read.table("TableS3.txt",header=T,stringsAsFactors=F,sep="\t")
dgefile=paste0("plot/Dinh_")
markers=dinhmarkers
markers$gene=markers$Gene
markers$cluster=markers$Cluster


# 
markers=markers$gene[which(markers$gene %in% genesall),]
genes=markers$gene
data.use=centroid.std

levels=colnames(centroid.std)

colsep.use=cumsum(table(levels))
col.lab=rep("",length(levels))
col.lab=levels

clusters=unique(markers$cluster)
rowsep.use=cumsum(table(markers$cluster)[clusters])
row.lab=1:length(clusters)
row.lab=rep("",nrow(markers))
row.lab[round(cumsum(table(markers$cluster)[clusters])-table(markers$cluster)[clusters]/2)+15]=clusters


ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=gg_color_hue(sum(ncluster))
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cell Type")
colnames(clab)=c("Cell Type","")

col.use=redblue100

data.use=centroid.std[genes,]
row.lab=rownames(data.use)

write.table(data.use,paste0(dgefile,"1_markers_DinhCentroidstd.txt"),row.names=T,col.names=T,quote=F,sep="\t")


jpeg(file=paste0(dgefile,"1_Dinh_centroid_std_markersall.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep = rowsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0(dgefile,"1_Dinh_centroid_std_markersall2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep = rowsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()

# saved as Figure S7B.