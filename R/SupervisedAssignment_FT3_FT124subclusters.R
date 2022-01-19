# 3.8.2021 by Qianyi
# supervised assignemnt of FT3 subcluster cells to epithelial subclusters of 3 large healthy FTs: FT1, FT2 and FT4


for the disease sample FT3:
supervised assign the subcluster cells of FT3 to subclusters of FT124
using HVG from post-CCA centroids of v4 FT124 subclusters 



library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

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


subsets=list(2:6,7:10,1:6,1)
subsetsname=c("2-6","7-10","1-6","1")


### load FT3 object of merged 2 segments
subjectorgans 
load(file=paste0(subjectorgans[4],".Robj"))
dgeall
# 27971 features across 15571 samples within 1 assay
# extract subset cells of FT3 assigned to FT124 using HVG of FT124 post-CCA centroids
classify0=read.table("plot/FT3_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
classify0[1:2,]

i=1 # did each iteration manually
i=2
i=4
cc=subsets[[i]]
ccname=subsetsname[i]
dgefile=paste0("/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/10xFallopialTube/plot/C",ccname,"_")

cells.use=rownames(classify0)[which(classify0$assign %in% cc)]
ft3data=dgeall@assays$RNA@data[,cells.use]
length(cells.use) # 6,081
table(classify0[cells.use,19])
   7    8    9   10 
2665 2727  443  246

  1 
559 


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

### supervised assignment of FT3 cells into FT124 C7-10 subcluster post-CCA centroids
#3. using HVG in FT124 post-CCAallgenes C7-10 subcluster centroids
g.mean <- apply(centroid,1,mean)
g.var <- apply(centroid,1,var)
g.var <- g.var[g.mean!=0]
g.mean <- g.mean[g.mean!=0]
if(i==1){
hvg2 <- names(g.mean)[g.mean > 0.03 & g.var/g.mean > 0.01]
print(length(hvg2)) # 3,745
}
if(i>1){
hvg2 <- names(g.mean)[g.mean > 0.1 & g.var/g.mean > 0.02]
print(length(hvg2))  # 3,737
}

length(which(hvg2 %in% rownames(ft3data))) 
genes=hvg2[which(hvg2 %in% rownames(ft3data))]
length(genes) # [1]3,735 4,057

d1=ft3data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") 

dd=rho[-c(1:ncol(centroid)),1:ncol(centroid)]
dd[1:2,]

###### Assign each cell of FT3 to each of the subclusters of FT124
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

   1    2    3    4    5 
2606 2788  387   48  252 

  1   3   4 
444  52  63 

classify3=classify
write.table(classify,paste0("plot/FT3_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"),row.names=T,col.names=T,quote=F,sep="\t")

# compare with global assignment for FT3
classify3=read.table(paste0("plot/FT3_C",ccname,"CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt"))
id3=classify3[,ncol(classify3)]
names(id3)=rownames(classify3)
table(id3)

id0=classify0[,19]
names(id0)=rownames(classify0)

table(id0[names(id3)],id3)

### assigned cluster ID for FT3
dgeall$assignedFT124 <- id0

dgeall$assignedFT124C2to6 <- id3

dgeall$assignedFT124C7to10 <- id3

dgeall$assignedFT124C1 <- id3

save(dgeall,file=paste0(subjectorgans[4],".Robj"))


