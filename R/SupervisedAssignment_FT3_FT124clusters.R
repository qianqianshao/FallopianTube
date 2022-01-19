# 1.18.2021 by Qianyi
# supervised assignemnt of FT3 to global 12 cell types of 3 large healthy FTs: FT1, FT2 and FT4

for the diseased sample FT3:
supervised assignment of FT3 cells based on the known centroid v4-x’_x4 from integrated FT1-2-4.
generate the v4-based supervised assignment for FT3 cells, and cross-tabulation with the v1 result for this sample. 
1. using all overlapped genes;
2. using HVG from FT3 cells;
3. using HVG from 17 cluster centroids of v4 - used this 
-> 17 clusters for FT3


R


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


### load FT3 object of merged 2 segments
subjectorgans 
load(file=paste0(subjectorgans[4],".Robj"))
dgeall
# 27971 features across 15571 samples within 1 assay
ft3data=dgeall@assays$RNA@data
hvg3=VariableFeatures(dgeall)

### load FT124 CCA-3subjects-allgenes post-CCA centroids
# object of FT124 CCA-3subjects-allgenes 
load(file=paste0(all,"_CCA-3subjects-allgenes.Robj"))
dgeall
table(Idents(dgeall)) # 17 clusters
# calculate cluster centroids
avg=AverageExpression(dgeall)
centroid=log(avg$RNA+1)
write.table(centroid,"plot/CCA3allgenes_ft124_Centroid.txt",row.names=T,col.names=T,quote=F,sep="\t")
centroid=log(avg$integrated+1)
write.table(centroid,"plot/CCA3allgenes_ft124_postCCA_Centroid.txt",row.names=T,col.names=T,quote=F,sep="\t")
# load post-CCAallgenes centroids 
centroid=read.table("plot/CCA3allgenes_ft124_postCCA_Centroid.txt")

### supervised assignment of FT3 cells into FT124 centroids
#1. using all overlapped genes
gene1=rownames(ft3data)
gene2=rownames(centroid)
genes=unique(intersect(gene1,gene2))
length(gene1) # [1] 27971
length(gene2) # [1] 31281
length(genes) # [1] 27422

d1=ft3data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using all genes
save(rho,file="FT3_allgenes_rankcor_FT124postCCAallgenescentroids.Robj")

#2. using HVG in FT3
length(which(hvg3 %in% rownames(centroid))) # 1997
genes=hvg3[which(hvg3 %in% rownames(centroid))]
length(genes) # [1] 1997

d1=ft3data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using HVG in FT3
save(rho,file="FT3_HVG_rankcor_FT124postCCAallgenescentroids.Robj")

#3. using HVG in FT124 post-CCAallgenes cluster centroids
g.mean <- apply(centroid,1,mean)
g.var <- apply(centroid,1,var)
g.var <- g.var[g.mean!=0]
g.mean <- g.mean[g.mean!=0]
hvg2 <- names(g.mean)[g.mean > 0.2 & g.var/g.mean > 0.05]
length(hvg2) # 3,083

length(which(hvg2 %in% rownames(ft3data))) # 3083
genes=hvg2[which(hvg2 %in% rownames(ft3data))]
length(genes) # [1]3083

d1=ft3data[genes,]
d2=centroid[genes,]
d1=as.matrix(d1)
d2=as.matrix(d2)
data=cbind(d2,d1)
data[1:2,1:20]

rho=cor(data,method="sp") # rank correlation using HVG in FT124 cluster centroids
save(rho,file="FT3_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.Robj")


dd=rho[-c(1:17),1:17]
dd[1:2,]

###### Assign each cell of FT3 to each of the 17 clusters of FT124
# if correlation with all 8 patterns <= 0.8, assign to the 9th pattern "all"
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

write.table(classify,"plot/FT3_allgenes_rankcor_FT124postCCAallgenescentroids.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(classify,"plot/FT3_HVG_rankcor_FT124postCCAallgenescentroids.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(classify,"plot/FT3_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt",row.names=T,col.names=T,quote=F,sep="\t")

# compare with global clusters of directly-merged FT1-5 (v1)
classify1=read.table("plot/FT3_allgenes_rankcor_FT124postCCAallgenescentroids.txt")
classify2=read.table("plot/FT3_HVG_rankcor_FT124postCCAallgenescentroids.txt")
id1=classify1[,19]
id2=classify2[,19]
names(id1)=rownames(classify1)
names(id2)=rownames(classify2)
table(id1)
table(id2)
table(id1,id2)
id1
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
 575   27 2025  407   19   48 2639 2785  358  230  475  458   74  286 1854  582 
  17 
2729 
id2
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
 584   82 2006  355   29   52 2625 2619  399  247  493  461   82  307 1809  603 
  17 
2818 
id1     1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
  1   575    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  2     0   27    0    0    0    0    0    0    0    0    0    0    0    0    0
  3     3   21 1989    2    6    2    0    0    0    0    0    0    0    0    0
  4     2   34   17  335    1    6    1    0    0    0    0    0    0    0    0
  5     0    0    0    0   19    0    0    0    0    0    0    0    0    0    0
  6     1    0    0    0    3   43    0    0    0    0    1    0    0    0    0
  7     0    0    0    5    0    0 2552   21    2    6    7    3   12    0    0
  8     2    0    0   11    0    0   72 2596   53    6    3    4    8    0    0
  9     0    0    0    1    0    0    0    2  342    6    0    1    3    0    0
  10    0    0    0    0    0    0    0    0    0  229    1    0    0    0    0
  11    0    0    0    0    0    0    0    0    0    0  475    0    0    0    0
  12    0    0    0    0    0    1    0    0    0    0    2  453    0    0    0
  13    0    0    0    1    0    0    0    0    2    0    0    0   59    0    0
  14    0    0    0    0    0    0    0    0    0    0    0    0    0  274    2
  15    0    0    0    0    0    0    0    0    0    0    0    0    0   33 1806
  16    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  17    1    0    0    0    0    0    0    0    0    0    4    0    0    0    1
    id2
id1    16   17
  1     0    0
  2     0    0
  3     0    2
  4     0   11
  5     0    0
  6     0    0
  7     7   24
  8    11   19
  9     1    2
  10    0    0
  11    0    0
  12    1    1
  13    4    8
  14    0   10
  15    0   15
  16  579    3
  17    0 2723


classify3=read.table("plot/FT3_FT124CentroidsHVG_rankcor_FT124postCCAallgenescentroids.txt")
id3=classify3[,19]
names(id3)=rownames(classify3)
table(id3)
table(id1,id3)
table(id2,id3)

v1=read.table("plot/FallopianTube1-5_ident.txt")
v1id=v1[,2]
names(v1id)=v1[,1]

table(v1id[names(id1)],id1)
table(v1id[names(id2)],id2)
table(v1id[names(id3)],id3)




dge=subset(dgeall,assignedFT124 != 13)
print(c(mean(dge@meta.data$nFeature_RNA),mean(dge@meta.data$nCount_RNA),mean(dge@meta.data$percent.mt)))
#[1]  2892.981644 12099.098953     4.131885



### visualize the fraction of cells in global clusters for healthy vs disease
myBrewerPalette=c(brewer.pal(12,"Paired")[1],"grey",brewer.pal(12,"Paired")[c(2:7,9:12)])
cols=c("#F8766D", "#078992") 

ncell=read.table("plot/Global_nCell.txt")[,2:1]
frac=prop.table(as.matrix(ncell),2)

library(reshape2)
data=melt(frac)
names(data)=c("ID","Status","Percentage")
data$label=round(data$Percentage*100,1)

png("plot/Global_frac_barplot.png",res=300,height=1200,width=1000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
barplot(frac,col=myBrewerPalette)
dev.off()

library(scales)
png("plot/Global_frac_barplot2.png",res=300,height=1200,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
ggplot(data, aes(x = Status, y = Percentage, fill = ID, label=label)) +
  geom_bar(stat = "identity") +
  geom_text(size=3,position = position_stack(vjust = 0.5))+
  scale_y_continuous(label=percent_format()) +
  scale_fill_manual(values=myBrewerPalette)+
  theme_bw()
dev.off()

png("plot/Global_frac_barplot33.png",res=300,height=1000,width=2000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
ggplot(data, aes(x = ID, y = Percentage, fill = Status, label=label)) +
  geom_col(position="dodge") +
  geom_text(size=3,position = position_dodge(0.9),
            vjust = 1.2,hjust = 0.5)+
  scale_y_continuous(label=percent_format()) +
  scale_fill_manual(values=rep(cols,6))+
  theme_bw()
dev.off()


png("plot/Global_frac_barplot3.png",res=300,height=1000,width=2000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
ggplot(data, aes(x = ID, y = Percentage, fill = Status, label=label)) +
  geom_col(position="dodge") +
  geom_text(size=3,position = position_dodge(0.9),
            vjust = 1.2,hjust = 0.5)+
  scale_y_continuous(label=percent_format()) +
  scale_fill_manual(values=rep(myBrewerPalette,each=2))+
  theme_bw()
dev.off()
