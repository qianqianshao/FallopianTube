# 8.16.2021 by Qianyi
# supervised assignemnt of FT6 to global 12 cell types of 3 large healthy FTs: FT1, FT2 and FT4

FT6: 1 whole fallopian tube from Sample3 (2518-AJ)  - included FT6, ovary, 2 Uterus parts from 1 caucasian woman, sequenced in 1 batch
For the 4th large healthy sample FT6:
supervised assignment of FT6 cells based on the known centroid v4-x’_x4 from FT1-2-4.
generate the v4-based supervised assignment for FT6 cells, and cross-tabulation with the CCA4-FT1246-integrated result for this sample. 
1. using all overlapped genes;
2. using HVG from FT6 cells;
3. using HVG from 17 cluster centroids of v4.

R


library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

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
all="FallopianTube124"


### load FT6 object 
subjectorgans 
load(file=paste0(subjectorgans[9],".Robj"))
dgeall
# 26320 features across 6659 samples within 1 assay 
FT6data=dgeall@assays$RNA@data
hvg6=VariableFeatures(dgeall)

### load FT124 CCA-3subjects-allgenes post-CCA centroids
# object of FT124 CCA-3subjects-allgenes 
load(file=paste0(all,"_CCA-3subjects-allgenes.Robj"))
dgeall
Idents(dgeall) <- dgeall$ft124CCA3all
table(Idents(dgeall)) # 17 clusters
# calculate cluster centroids
avg=AverageExpression(dgeall)
centroid=log(avg$RNA+1)
write.table(centroid,"plot/CCA3allgenes_ft124_Centroid.txt",row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,"plot/CCA3allgenes_ft124_postCCA_Centroid.txt",row.names=T,col.names=T,quote=F,sep="\t")
# load post-CCAallgenes centroids 
centroid=read.table("plot/CCA3allgenes_ft124_postCCA_Centroid.txt")

### supervised assignment of FT6 cells into FT124 centroids
#1. using all overlapped genes
gene1=rownames(FT6data)
gene2=rownames(centroid)
genes=unique(intersect(gene1,gene2))
length(gene1) # [1] 26320
length(gene2) # [1] 31281
length(genes) # [1] 26007

d1=FT6data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using all genes
save(rho,file="FT6_allgenes_rankcor_FT124postCCAallgenescentroids.Robj")

#2. using HVG in FT6
length(which(hvg6 %in% rownames(centroid))) # 2000
genes=hvg6[which(hvg6 %in% rownames(centroid))]
length(genes) # [1] 2000

d1=FT6data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using HVG in FT6
save(rho,file="FT6_HVG_rankcor_FT124postCCAallgenescentroids.Robj")

#3. using HVG in FT124 post-CCAallgenes cluster centroids
g.mean <- apply(centroid,1,mean)
g.var <- apply(centroid,1,var)
g.var <- g.var[g.mean!=0]
g.mean <- g.mean[g.mean!=0]
hvg2 <- names(g.mean)[g.mean > 0.2 & g.var/g.mean > 0.05]
length(hvg2) # 3,083

length(which(hvg2 %in% rownames(FT6data))) # 3083
genes=hvg2[which(hvg2 %in% rownames(FT6data))]
length(genes) # [1]3083

d1=FT6data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using HVG in FT124 cluster centroids
save(rho,file="FT6_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.Robj")


dd=rho[-c(1:17),1:17]
dd[1:2,]

###### Assign each cell of FT6 to each of the 17 clusters of FT124
maxcorthres=0
#maxcorthres=0.9
ncell=ncol(d1)
n=ncol(dd)
classify=matrix(0,nrow=ncell,ncol=n+2)
rownames(classify)=rownames(dd)
colnames(classify)=c(colnames(dd),"maxcor","assign")
classify[,1:n]=dd
for(i in 1:ncell){
  classify[i,18]=max(classify[i,1:n],na.rm=T)
  if(classify[i,18]<=maxcorthres){
    classify[i,19]=18
  } else  {
    classify[i,19]=which(classify[i,1:n]==classify[i,18]) 
  }  
}
summary(classify[,18])
summary(classify[,19])
table(classify[,19])

write.table(classify,"plot/FT6_allgenes_rankcor_FT124postCCAallgenescentroids.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(classify,"plot/FT6_HVG_rankcor_FT124postCCAallgenescentroids.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(classify,"plot/FT6_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt",row.names=T,col.names=T,quote=F,sep="\t")

# compare with global clusters of directly-merged FT1-5 (v1)
classify1=read.table("plot/FT6_allgenes_rankcor_FT124postCCAallgenescentroids.txt")
classify2=read.table("plot/FT6_HVG_rankcor_FT124postCCAallgenescentroids.txt")
id1=classify1[,19]
id2=classify2[,19]
names(id1)=rownames(classify1)
names(id2)=rownames(classify2)
table(id1)
table(id2)
id1=factor(id1,levels=1:17)
id2=factor(id2,levels=1:17)
table(id1,id2)



classify3=read.table("plot/FT6_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
id3=classify3[,19]
names(id3)=rownames(classify3)
id3=factor(id3,levels=1:17)
table(id3)
table(id1,id3)
table(id2,id3)


v1=read.table("plot/FallopianTube1246_CCA4-allgenes_ident.txt")
v1id=v1[,2]
names(v1id)=v1[,1]

table(v1id[names(id1)],id1)
table(v1id[names(id2)],id2)
table(v1id[names(id3)],id3)




