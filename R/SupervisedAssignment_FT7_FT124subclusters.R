# 11.8.2021 by Qianyi
# supervised assignemnt of FT7 subcluster cells to epithelial subclusters of 3 large healthy FTs: FT1, FT2 and FT4


for the disease sample FT7:
supervised assignment of FT7 cells based on the known centroid v4-x’_x4 from samples 1-2-4.
-> 17 clusters for FT7
supervised assign the epithelial cells of FT7 to subclusters of epithelial subclusters of FT124 
using HVG from post-CCA centroids of v4 epithelial subclusters 


R


library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3","2518-AJ-1","2518-AJ-2","2518-AJ-3","2518-AJ-4","2788-NU-1","2788-NU-2","4360-YS-1")
dataset1=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1","1427-NU-1","1427-NU-2","1427-NU-3","1407-NU-1","1407-NU-2","1714-YS-1","1714-YS-2","1714-YS-3","1847-NU-1","1847-NU-2","1847-NU-3","2518-AJ-1_CACTACGA-ATCAGTCT","2518-AJ-2_CACGGTGA-TGTGACGA","2518-AJ-3_ATGGCTTG-CACAACAT","2518-AJ-4_CCTTCTAG-TCGTTGTA","2788-NU-1_TATCAGCC-AGGACGAA","2788-NU-2_TGGTCCCA-ACGCCAGA","Sample_4360-YS-1_GCTACAAA-AGGGCACG")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium1","Fimbria2","Ampulla2","Isthmus2","FimAmp3","Isthmus3","Fimbria4","Ampulla4","Isthmus4","Uterus2","FT5","Ovary1","FT7","Ovary2","Myometrium2","Endometrium2","Myometrium3","Endometrium3","FT7")
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

all="FallopianTube124"


subsets=list(1,2:6,7:10)
subsetsname=c("1","2-6","7-10")



### load FT7 object 
subjectorgans 
load(file=paste0(subjectorgans[13],".Robj"))
dgeall
# extract subset cells of FT7 assigned to FT124 using HVG of FT124 post-CCA centroids
classify0=read.table("plot/FT7_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
classify0[1:2,]

i=1  # did each iteration manually
i=2

cc=subsets[[i]]
ccname=subsetsname[i]
dgefile=paste0("/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/10xFallopialTube/plot/C",ccname,"_")

cells.use=rownames(classify0)[which(classify0$assign %in% cc)]
ft7data=dgeall@assays$RNA@data[,cells.use]
length(cells.use) # 
table(classify0[cells.use,19])
  1
143

  2   3   4   6 
  1 194  11   1

### load FT124 CCA-3subjects-allgenes subcluster post-CCA centroids
# object of FT124 CCA-3subjects-allgenes 
load(file=paste0(all,"_C",ccname,".Robj"))
dge
table(Idents(dge)) # 5 subclusters
# calculate cluster centroids
avg=AverageExpression(dge)
centroid=log(avg$RNA+1)
write.table(centroid,paste0("plot/CCA3allgenes_ft124C",ccname,"_Centroid.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,paste0("plot/CCA3allgenes_ft124C",ccname,"_postCCA_Centroid.txt"),row.names=T,col.names=T,quote=F,sep="\t")
# load post-CCAallgenes centroids 
centroid=read.table(paste0("plot/CCA3allgenes_ft124C",ccname,"_postCCA_Centroid.txt"))

### supervised assignment of FT7 cells into FT124 C7-10 subcluster post-CCA centroids
#3. using HVG in FT124 post-CCAallgenes C7-10 subcluster centroids
g.mean <- apply(centroid,1,mean)
g.var <- apply(centroid,1,var)
g.var <- g.var[g.mean!=0]
g.mean <- g.mean[g.mean!=0]
hvg2 <- names(g.mean)[g.mean > 0.1 & g.var/g.mean > 0.02]
print(length(hvg2))  # 4,057
if(i==2){
hvg2 <- names(g.mean)[g.mean > 0.03 & g.var/g.mean > 0.01]
print(length(hvg2)) # 3,745
}

length(which(hvg2 %in% rownames(ft7data))) 
genes=hvg2[which(hvg2 %in% rownames(ft7data))]
length(genes) # [1]4,057 3,742

d1=ft7data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using HVG in FT124 C7-10 subcluster post-CCA centroids
#save(rho,file="FT7_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.Robj")


dd=rho[-c(1:ncol(centroid)),1:ncol(centroid)]
dd[1:2,]

###### Assign each cell of FT7 to each of the subclusters of FT124
maxcorthres=0
#maxcorthres=0.9
ncell=ncol(d1)
n=ncol(dd)
classify=matrix(0,nrow=ncell,ncol=n+2)
rownames(classify)=rownames(dd)
colnames(classify)=c(colnames(dd),"maxcor","assign")
classify[,1:n]=dd
for(i in 1:ncell){
  classify[i,n+1]=max(classify[i,1:n],na.rm=T)
  if(classify[i,n+1]<=maxcorthres){
    classify[i,n+2]=n+1
  } else  {
    classify[i,n+2]=which(classify[i,1:n]==classify[i,n+1]) 
  }  
}
summary(classify[,n+1])
summary(classify[,n+2])
table(classify[,n+2])

  1   2   3 
129  11   3 

  1   2   3   5   6 
  1   1 196   8   1
  
classify3=classify
write.table(classify,paste0("plot/FT7_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"),row.names=T,col.names=T,quote=F,sep="\t")

# compare with global assignment for FT7
classify3=read.table(paste0("plot/FT7_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id3=classify3[,ncol(classify3)]
names(id3)=rownames(classify3)
table(id3)

id0=classify0[,19]
names(id0)=rownames(classify0)

table(id0[names(id3)],id3)


### assigned cluster ID for FT7
dgeall$assignedFT124 <- id0

dgeall$assignedFT124C1 <- id3

dgeall$assignedFT124C2to6 <- id3

save(dgeall,file=paste0(subjectorgans[13],".Robj"))


