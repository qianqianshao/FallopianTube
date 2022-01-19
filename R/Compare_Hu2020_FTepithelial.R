###### Correlation between our fallopian tubes (healthy FT1,2,4) and 5 cell type centroids for FT in ovarian cancer
# by Qianyi on 11/27/2021
# Related to Figure S7C.

### Integrated 3 large healthy FTs by Qianyi in Dec 2020
# 10x 


### 5 clusters for FT in Ovarian cancer
# centroid of 5 cluster centroids
centroid-mmc7.txt
data format: log(mean(counts-per-10k)+1), 52 markers, 5 cluster centroids
# markers of secretory vs ciliate cell types
markers-mmc5.txt
# markers of secretory cell subtypes
markers-mmc6.txt



###### load our epithelial subset of healthy FT data analyzed by Qianyi in 2021
library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)
all="FallopianTube1234"

load(file=paste0(all,"_CCA-4subjects-allgenes.Robj")) # used this
dgeall

dge1=subset(dgeall,ft124CCA3allFT3assign %in% c(1:6))
dge=subset(dge1,Status == "Healthy")
Idents(dge) <- dge$ft124CCA3allFT3assign21
dge12=dge

ft=log(AverageExpression(dge12)$RNA+1)
write.table(ft,file="FallopianTube124_C1-6_centroid.txt",quote=F,row.names=T,col.names=T)


load(file=paste0("FallopianTube124_C1-6.Robj"))

ft=read.table("FallopianTube124_C1-6_centroid.txt")
colnames(ft)=gsub("X","",colnames(ft))

### load markers for epithelial subclusters
ftmarker1=read.table("plot/C1_8PCs_integrated_snn_res.0.1ordered_logfc1.6fold_min0.1_RNAassay_1.2021.txt",header=T,stringsAsFactors=F,sep="\t")
ftmarker2=read.table("plot/C2-6_6PCs_integrated_snn_res.0.1ordered_logfc1.6fold_min0.1_RNAassay_1.2021.txt",header=T,stringsAsFactors=F,sep="\t")
dim(ftmarker1) # 465     7
dim(ftmarker2) # 574     7
ftmarker1$cluster = paste0("1_",ftmarker1$cluster)
ftmarker2$cluster = paste0("2_",ftmarker2$cluster)
ftmarkers=rbind(ftmarker1,ftmarker2)
### number of markers for each cluster
table(ftmarkers$cluster)
1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 2_5 2_6 
137  30 156 142 115  90  13  27 287  42

### load 5 cluster centroids for published ovarian cancer FT
# data format: log(mean(counts-per-10k)+1), 52 genes, 5 cluster centroids (columns)
oc=read.table("centroid-mmc7.txt",header=T)
dim(oc)   # [1] 52 5
ocmarker1=read.table("markers-mmc5.txt",header=T)
ocmarker2=read.table("markers-mmc6.txt",header=T)




### Cross-tabulate using rank correlation 
markers=rownames(oc)[which(rownames(oc) %in% rownames(ft))]
length(markers) # 51
  all=cbind(ft[markers,],oc[markers,])
  rho=cor(all,method="spearman")
rho[1:10,11:15]
             C3           C4         EMT          C9     Ciliated
1_1  0.06873459 -0.004977488 -0.06408835 -0.19344034  0.834464661
1_2 -0.13498043 -0.121857961 -0.06259477 -0.26202497  0.586930946
1_3  0.12683545  0.024842191 -0.06870488 -0.15722404  0.705131689
1_4 -0.19353379 -0.175071834  0.27590578 -0.57982306  0.346185177
2_1  0.42100500  0.443041698  0.75186135  0.08248264 -0.125803241
2_2  0.25263015  0.335075454  0.81871056 -0.07229680 -0.274866505
2_3  0.53756872  0.545668454  0.43504041  0.35845089 -0.004796814
2_4  0.57232064  0.670965407  0.26717057  0.27850339  0.190605485
2_5  0.54503496  0.682504129  0.30826677  0.39584422  0.148610735
2_6  0.52720650  0.596981832  0.41114306  0.33830557  0.011132229
# saved as Figure S7C.

### check if markers overlap
for(i in 1:10){
  genes=ftmarkers[which(ftmarkers$cluster==unique(ftmarkers$cluster)[i]),]$gene
  print(length(which(genes %in% rownames(oc)))
    )
}
