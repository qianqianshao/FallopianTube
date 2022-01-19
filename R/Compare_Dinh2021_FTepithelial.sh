###### Correlation between our fallopian tubes (healthy FT1,2,4) and cell type centroids for FT from Huy Dinh 2021 paper
# by Qianyi on 7/19/2021
# Related to Figure S7A.

### Integrated 3 large healthy FTs by Qianyi in Dec 2020
# 10x 


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
all="FallopianTube1234"

load(file=paste0(all,"_CCA-4subjects-allgenes.Robj")) # used this
dgeall

dge1=subset(dgeall,ft124CCA3allFT3assign %in% c(1:6))
dge=subset(dge1,Status == "Healthy")
Idents(dge) <- dge$ft124CCA3allFT3assign21
dge12=dge

ft=log(AverageExpression(dge12)$RNA+1)

load(file=paste0("FallopianTube124_C1-6.Robj"))


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

### load 10 cluster centroids for published Dinh FT
# data format: log(mean(counts-per-10k)+1), 10 cluster centroids (columns)
dinh=read.csv("avg_expr_epithelials_July-2021.csv",row.names=1)
#ft=read.table("ftGonad4ClusterCentroids_RPKM.txt",header=T)
dim(dinh) #  45037    10

### load markers for 10 clusters of Dinh FT data
dinhmarkers=read.table("TableS3.txt",header=T,stringsAsFactors=F,sep="\t")
#ft=read.table("ftGonad4ClusterCentroids_RPKM.txt",header=T)
dim(dinhmarkers) # 1467     7
### number of markers for each cluster
table(dinhmarkers$Cluster)
Ciliated_1      Ciliated_2      Ciliated_3      Ciliated_4     Secretory_1 
             88             108             139             100             129 
    Secretory_2     Secretory_3 Unclassified _1  Unclassified_2  Unclassified_3 
            107             137             186             193             280 

# check the number of markers overlapped for each cluster
markerlist1=markerlist2=list()
n=10
for(i in 1:n){
  markerlist1[[i]]=ftmarkers$gene[which(ftmarkers$cluster==names(table(ftmarkers$cluster))[i])]
  markerlist2[[i]]=dinhmarkers$Gene[which(dinhmarkers$Cluster==names(table(dinhmarkers$Cluster))[i])]
}
cc=matrix(0,n,n)
for(i in 1:n){
    for(j in 1:n){
        cc[i,j]=length(intersect(markerlist1[[i]],markerlist2[[j]]))
    }
}

### Cross-tabulate using rank correlation 
genes=rownames(dinh)[which(rownames(dinh) %in% rownames(ft))]
length(genes) # 28241

  all=cbind(ft[genes,],dinh[genes,])

# using HVG from our epithelial subset
hvg=VariableFeatures(dge)
hvg=hvg[which(hvg %in% genes)]
  rho=cor(all[hvg,],method="spearman")
rho[11:20,1:10]



# using markers from Dinh paper
markers=dinhmarkers$Gene
markers=markers[which(markers %in% genes)]
length(markers) # 1,442
  rho=cor(all[markers,],method="spearman")
rho[11:20,1:10]
                                  1_1       1_2       1_3       1_4       2_1
X2_Secretory_Epithelials    0.6942236 0.7376651 0.8416488 0.7377539 0.8645224
X1_Secretory_Epithelials    0.6815364 0.7278704 0.8141873 0.7189718 0.7957961
X2_Ciliated_Epithelials     0.9507778 0.8797735 0.8317918 0.8481561 0.5162866
X1_Ciliated_Epithelials     0.9081996 0.8509415 0.7760591 0.8085931 0.4359168
X3_Secretory_Epithelials    0.6482817 0.7074291 0.8257318 0.7106532 0.8814043
X1_Transitional_Epithelials 0.6358768 0.6248453 0.5950078 0.6033609 0.4206668
Early_Secretory             0.6977588 0.7343850 0.7329674 0.8723864 0.6882158
X4_Ciliated_Epithelials     0.9222530 0.8500560 0.7779337 0.8087853 0.4277237
X2_Transitional_Epithelials 0.2784917 0.3629164 0.4184073 0.3308890 0.5746281
X3_Ciliated_Epithelials     0.8339618 0.7938415 0.7315550 0.7801790 0.4348259
                                  2_2       2_3       2_4       2_5       2_6
X2_Secretory_Epithelials    0.8459767 0.8907274 0.8942928 0.8987023 0.8603897
X1_Secretory_Epithelials    0.7731956 0.8232447 0.8238091 0.8093182 0.8417958
X2_Ciliated_Epithelials     0.4889690 0.5646498 0.5780590 0.5537098 0.4794281
X1_Ciliated_Epithelials     0.4088137 0.4808657 0.4914257 0.4686169 0.4267109
X3_Secretory_Epithelials    0.8651689 0.8971049 0.9151770 0.8912281 0.8868689
X1_Transitional_Epithelials 0.3791467 0.4255729 0.4338247 0.4127900 0.5759487
Early_Secretory             0.7468043 0.6583295 0.6577118 0.6241493 0.5925128
X4_Ciliated_Epithelials     0.4038949 0.4807033 0.4908731 0.4758300 0.3939784
X2_Transitional_Epithelials 0.5317800 0.5720477 0.5508882 0.5431160 0.7275858
X3_Ciliated_Epithelials     0.4128317 0.4612808 0.4723347 0.4290862 0.4578204
# saved as Figure S7A.

# using markers from our paper
markers=ftmarkers$gene
markers=markers[which(markers %in% genes)]
  rho=cor(all[markers,],method="spearman")
rho[11:20,1:10]


rho[1:10,1:10]


marker1=ftmarker1$gene
marker1=marker1[which(marker1 %in% genes)]
marker2=ftmarker2$gene
marker2=marker2[which(marker2 %in% genes)]
  rho=cor(all[marker1,],method="spearman")
rho[11:20,1:4]

  rho=cor(all[marker2,],method="spearman")
rho[11:20,5:10]
