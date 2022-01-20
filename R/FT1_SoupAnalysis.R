# Soup analysis for each of the 3 segments of FT1 on 7.20.2020 by Qianyi
# 1) using SoupX
# 2) our modified Soup analysis

R

library(SoupX)

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
myBrewerPalette=c(brewer.pal(7,"Set2"),brewer.pal(12,"Paired")[c(10,12)])
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

dgefile=getwd()
setwd(dgefile)
dataset=c("747-YS-1","747-YS-2","747-YS-3","1073-NU-1")
n=length(dataset)
names=c("Fimbria1","Ampulla1","Isthmus1","Myometrium2")
part=c("Fimbria","Ampulla","Isthmus","Myometrium")
organ=c("FallopianTube","FallopianTube","FallopianTube","Uterus")
run=c("NovaA-223","NovaA-223","NovaA-223","NovaA-237")
menopause=c("pre","pre","pre","post")
subject=c("Human1","Human1","Human1","Human2")
datainfo=data.frame(run,dataset,names,menopause,subject,organ,part)
subjectorgans=unique(paste(subject,organ,sep="-"))
organs=c("FallopianTube","Uterus")




### load directly-merged 3 segments of fallopian tube1
indiv=1
  load(file=paste0(subjectorgans[1],".Robj"))  
dge=dgeall



### 1. using SoupX
# estimate soup for our filtered cells with >500 genes and <10% MT and our clusters - use this
library(SoupX)
sclist=list()
### 1 dataset
i=1
i=2
i=3
for(i in 1:3){
dataDir=paste0("/nfs/turbo/umms-hammou/10x/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/")
tod = Seurat::Read10X(file.path(dataDir,'raw_feature_bc_matrix'))
colnames(tod)=paste(names[i],colnames(tod),sep="_")
tocall=GetAssayData(dge,slot="counts")
which(colnames(tocall) != rownames(dge@meta.data))
toc1=tocall[,which(dge$Part==part[i])]
toc0gene=rownames(tod)[which(!(rownames(tod) %in% rownames(toc1)))]
toc0=matrix(0,nrow=length(toc0gene),ncol=ncol(toc1))
rownames(toc0)=toc0gene
colnames(toc0)=colnames(toc1)
toc=rbind(toc1,toc0)
clusters=dge$ordered[which(dge$Part==part[i])]
id=data.frame(clusters=clusters)
sc = SoupChannel(tod,toc,metaData=id)
pdf(paste0("plot/Soup_",part[i],"_autoEstCont.pdf"),width=6,height=4)
sc = autoEstCont(sc)
dev.off()
pdf(paste0("plot/Soup_",part[i],"_autoEstCont2.pdf"),width=6,height=4)
sc = autoEstCont(sc,contaminationRange=c(0.001,0.8)) # allow smaller min 
dev.off()
out = adjustCounts(sc)

save(sc,file=paste0("Soup_",names[i],".Robj"))
sclist[[i]]=sc
}

3117 genes passed tf-idf cut-off and 711 soup quantile filter.  Taking the top 100.
Using 1271 independent estimates of rho.
Estimated global rho of 0.01
Expanding counts from 15 clusters to 1861 cells.

3951 genes passed tf-idf cut-off and 934 soup quantile filter.  Taking the top 100.
Using 1048 independent estimates of rho.
Estimated global rho of 0.01
Expanding counts from 15 clusters to 4944 cells.

2354 genes passed tf-idf cut-off and 419 soup quantile filter.  Taking the top 100.
Using 969 independent estimates of rho.
Estimated global rho of 0.01
Expanding counts from 15 clusters to 3722 cells.


library(ggplot2)
umap=Embeddings(dge,"umap")
umap=umap[colnames(sc$toc),]
sc=setDR(sc,umap)
dd=sc$metaData
mids = aggregate(cbind(UMAP_1,UMAP_2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(UMAP_1,UMAP_2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle(part[i]) + 
    guides(colour = guide_legend(override.aes = list(size = 1)))+
    scale_colour_manual(values=myBrewerPalette)
plot(gg)

dd$VIM = sc$toc['VIM',]
gg = ggplot(dd,aes(UMAP_1,UMAP_2)) +
  geom_point(aes(colour=VIM>0))
plot(gg)

# assuming each droplet contains nothing but background contamination
# visualize the ratio of observed counts for a markers to this expectation value
gg = plotMarkerMap(sc, "VIM")
plot(gg)
gg = plotMarkerMap(sc, "KRT7")
plot(gg)

### Estimating the contamination fraction 
# estimating the level of background contamination (represented as rho) 
# idea: identify genes that should be not expressed by some cells in our data and the expression that we observe for these genes in these cells must be due to contamination.
plotChangeMap(sc, out, "VIM")
plotChangeMap(sc, out, "KRT7")


pdf(paste0("plot/Soup_",part[i],"_markerDistribution.pdf"),width=5,height=4)
plotMarkerDistribution(sc)
dev.off()

pdf(paste0("plot/Soup_",part[i],"_changemap_knownmarkers.pdf"),width=5,height=4)
plotChangeMap(sc, out, gene)
dev.off()

### Investigating changes in expression
# check the fraction of cells that were non-zero now set to zero after correction
library(Matrix)
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

# check genes for which there is a quantative difference after correction
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
# mostly mitochondrial genes 




### 2. our modified soup analysis
### 1) estimated soup for our filtered cells with >500 genes and <10% MT and our clusters
sclist=list()
### 1 dataset
i=1
i=2
i=3
load(file=paste0("Soup_",names[i],".Robj"))

dataDir=paste0("/nfs/turbo/umms-hammou/10x/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/")
tod = Seurat::Read10X(file.path(dataDir,'raw_feature_bc_matrix'))
colnames(tod)=paste(names[i],colnames(tod),sep="_")
dataDir=paste0("/nfs/turbo/umms-hammou/10x/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/")
sc11 <- load10X(dataDir)
toc11 = sc11$toc
tocall=GetAssayData(dge,slot="counts")
which(colnames(tocall) != rownames(dge@meta.data))
toc1=tocall[,which(dge$Part==part[i])]
toc0gene=rownames(tod)[which(!(rownames(tod) %in% rownames(toc1)))]
toc0=matrix(0,nrow=length(toc0gene),ncol=ncol(toc1))
rownames(toc0)=toc0gene
colnames(toc0)=colnames(toc1)
toc=rbind(toc1,toc0)
clusters=dge$new1[which(dge$Part==part[i])]
id=data.frame(clusters=clusters)

### Visualize nUMI distribution for all cell barcodes with at least 1 UMI
nDropUMIs=sort(sc$nDropUMIs,decreasing=TRUE)
nDropUMIs<- nDropUMIs[which(nDropUMIs>0)] # cell barcodes with >=1UMI
cell11=paste(names[i],colnames(toc11),sep="_") # Cells selected from nUMI distribution plot
cell1=colnames(toc1)   # Filtered cells with >500 detected genes and <10% MT transcripts
print(c(length(nDropUMIs),length(cell11),length(cell1)) 
# 753049   2338   1861
print(c(max(nDropUMIs[!(names(nDropUMIs) %in% cell11)]),min(nDropUMIs[cell11]),min(nDropUMIs[cell1])))
bc=1:length(nDropUMIs)
names(bc)=names(nDropUMIs)
jpeg(paste0("plot/nDropUMIs_",part[i],".jpeg"),width=1500,height=1200,res=300)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(bc,nDropUMIs,log="xy",cex=0.5,pch=16,col=rgb(0,0,0,0.5),xlab="Barcodes",ylab="nUMI")
points(bc[cell11],nDropUMIs[cell11],cex=0.5,pch=16,col=rgb(0,1,1,0.7))
points(bc[cell1],nDropUMIs[cell1],cex=0.5,pch=16,col=rgb(1,0,0,0.9))
legend("topright",legend=c("Filtered Cells","Cells","Background"),col=c("red","cyan","black"),pch=16)
dev.off()
c=cumsum(nDropUMIs)
y=c/max(c)
jpeg(paste0("plot/nDropUMIs_",part[i],"_CumulativeDistribution.jpeg"),width=1500,height=1200,res=300)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(bc,y,cex=0.5,pch=16,col=rgb(0,0,0,0.5),xlab="Barcodes with Decreasing Order",ylab="Fraction of Cumulative nUMI")
points(bc[cell11],y[cell11],cex=0.5,pch=16,col=rgb(0,1,1,0.7))
points(bc[cell1],y[cell1],cex=0.5,pch=16,col=rgb(1,0,0,0.9))
legend("bottomright",legend=c("Filtered Cells","Cells","Background"),col=c("red","cyan","black"),pch=16)
dev.off()


soupRange=c(0,100);keepDroplets=FALSE # default
estimateSoup = function(sc,soupRange=c(0,100),keepDroplets=FALSE){
  if(!is(sc,'SoupChannel'))
    stop("sc must be a SoupChannel object.")
  #Estimate the soup 
  w = which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
  sc$soupProfile = data.frame(row.names=rownames(sc$tod),
                              est = rowSums(sc$tod[,w,drop=FALSE])/sum(sc$tod[,w]),
                              counts = rowSums(sc$tod[,w,drop=FALSE]))
  #Saves a lot of space if we can drop the droplets now we're done with them
  if(!keepDroplets)
    sc$tod=NULL
  return(sc)
}

### decided to divide cells and soup by nUMI<100 and select 2000-3000 barcodes with <100 UMIs as soup
nDropUMIs=sort(sc$nDropUMIs,decreasing=TRUE)
nDropUMIs<- nDropUMIs[which(nDropUMIs>0)] # cell barcodes with >=1UMI
cells=names(nDropUMIs)[which(nDropUMIs>=100)]
length(cells)  # 2875
soup=names(nDropUMIs)[which(nDropUMIs<100)][1:3000]
summary(nDropUMIs[soup]) # 23 - 99 UMIs per cell barcode
table(nDropUMIs[soup])
# choose barcodes with >23 and <100 UMIs as soup
which(nDropUMIs<=23)[1]-1 # 5786
which(nDropUMIs<=23)[1]-1-length(cells) # 2911
soup=names(nDropUMIs)[which(nDropUMIs<100)][1:2911]
summary(nDropUMIs[soup]) # 23 - 99 UMIs per cell barcode
table(nDropUMIs[soup])
length(soup)   # 2911
print(c(length(nDropUMIs),length(cells),length(soup)) )
# 753049   2875   2911
bc=1:length(nDropUMIs)
names(bc)=names(nDropUMIs)
jpeg(paste0("plot/nDropUMIs_",part[i],"_3kSoupbc.jpeg"),width=1500,height=1200,res=300)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(bc,nDropUMIs,log="xy",cex=0.5,pch=16,col=rgb(0,0,0,0.5),xlab="Barcodes",ylab="nUMI")
abline(v=length(cells),lty=2)
dev.off()
c=cumsum(nDropUMIs)
y=c/max(c)
jpeg(paste0("plot/nDropUMIs_",part[i],"_3kSoupbc_CumulativeDistribution.jpeg"),width=1500,height=1200,res=300)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(bc,y,cex=0.5,pch=16,col=rgb(0,0,0,0.5),xlim=c(0,50000),xlab="Barcodes with Decreasing Order",ylab="Fraction of Cumulative nUMI")
abline(v=length(cells),lty=2)
dev.off()

summary(nDropUMIs[soup]) # 24 - 99 UMIs per cell barcode
table(nDropUMIs[soup])
# based on this, mofiy soup range to be 23-100 (>23 and <100UMIs) instead of 0-100
soupRange=c(23,100)

# 1) Setting up SoupChannel using ~3k soup barcodes
tod1=tod[,c(cells,soup)]
dim(tod1) # 33538 5786
sc = SoupChannel(tod1,toc,metaData=id)
soupRange=c(23,100);keepDroplets=FALSE
#estimateSoup = function(sc,soupRange=c(0,100),keepDroplets=FALSE){
  #Estimate the soup 
  w = which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
  summary(w) # 2876-5786
  length(w)  # 2911
  sc$soupProfile = data.frame(row.names=rownames(tod1),
                              est = rowSums(tod1[,w,drop=FALSE])/sum(tod1[,w]),
                              counts = rowSums(tod1[,w,drop=FALSE]))
# est: total expression of each gene in soup barcodes divided by total expression of all genes in soup barcodes
# counts: total expression of each gene in soup barcodes

### 2) estimate contamination fraction
pdf(paste0("plot/Soup_",part[i],"_3kSoup_autoEstCont2.pdf"),width=6,height=4)
sc = autoEstCont(sc,contaminationRange=c(0.001,0.8)) # allow smaller min 
dev.off()
3777 genes passed tf-idf cut-off and 770 soup quantile filter.  Taking the top 100.
Using 1743 independent estimates of rho.

save(sc,file=paste0("Soup_",names[i],"_3kSoup.Robj"))

#out = adjustCounts(sc)

# autoEstCont
sc=sc;topMarkers=NULL;tfidfMin=1.0;soupQuantile=0.90;maxMarkers=100;
contaminationRange=c(0.001,0.8);rhoMaxFDR=0.2;priorRho=0.05;priorRhoStdDev=0.10;doPlot=TRUE;forceAccept=FALSE;verbose=TRUE
  if(!'clusters' %in% colnames(sc$metaData))
    stop("Clustering information must be supplied, run setClusters first.")
  #First collapse by cluster
  s = split(rownames(sc$metaData),sc$metaData$clusters)
  # total expression of each gene for each cluster
  tmp = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
  ssc = sc 
  ssc$toc = tmp
  dim(ssc$toc) # [1] 33538    19
  ssc$metaData = data.frame(nUMIs = colSums(tmp),row.names=colnames(tmp))
  ###################
  # Get best markers
  #Get the top N soup Genes
  soupProf = ssc$soupProfile[order(ssc$soupProfile$est,decreasing=TRUE),]
  soupMin = quantile(soupProf$est,soupQuantile)
  soupProf[1:5,]
#               est counts
#MT-ND3  0.02800058   3107
#MT-ATP6 0.02515275   2791
#MT-CO2  0.01767272   1961
#MT-CO1  0.01690669   1876
#MT-ND4  0.01519439   1686
  soupMin # get top 10% of soup genes with est > soupMin
#         90% 
#3.604838e-05
  #Find or load markers.
  if(is.null(topMarkers)){
    #Refine this to the best markers we can manage
    mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
# quickMarkers
# gene expression is binarised in each cell so each cell is either considered to express or not each gene.  
# i.e., we replace the counts with toc > zeroCut
# The frequency with which a gene is expressed within the target group is compared to the global frequency to calculate the tf-idf score. 
# tf-idf: Term Frequency - Inverse Document Frequency 
#1-4_1 1-4_3 1-4_4 1-4_5 1-4_6  10_1    11    12    13    14    15     5     6 
# 5679  5845  8229 10456   103   184   477   410   200  1052   507   841    77 
#    7     8     9 
#  214   579  1030 
    #And only the most specific entry for each gene
    mrks = mrks[order(mrks$gene,-mrks$tfidf),]
    mrks = mrks[!duplicated(mrks$gene),]
    #Order by tfidif maxness
    mrks = mrks[order(-mrks$tfidf),]
    #Apply tf-idf cut-off
    mrks = mrks[mrks$tfidf > tfidfMin,]
table(mrks$cluster)
#1-4_1 1-4_3 1-4_4 1-4_5 1-4_6  10_1    11    12    13    14    15     5     6 
# 1461   706    10   251    75   166    16    26    93   242   315     1    74 
#    7     8     9 
#   77   123   141 
summary(mrks$tfidf)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.088   1.255   1.488   1.664   5.131
  }else{
    mrks = topMarkers
  }
  #Filter to include genes that exist in soup - top 10% genes with high fraction of expression in soup
  tgts = rownames(soupProf)[soupProf$est>soupMin]
  length(tgts) #[1] 2978

  #And get the ones that pass our tfidf cut-off
  filtPass = mrks[mrks$gene %in% tgts,]
  dim(filtPass) # [1] 770  10
  table(filtPass$cluster)
#1-4_1 1-4_3 1-4_4 1-4_5 1-4_6  10_1    11    12    13    14    15     6     7 
#  305   183     3   107    37    25     1     2     7    21    27    13     9 
#    8     9 
#   11    19
  maxMarkers # 100
  table(filtPass$cluster[1:100])
#1-4_1  10_1    14    15     6     7     9 
#   79     7     3     4     1     3     3 
  tgts = head(filtPass$gene,n=maxMarkers)

  if(verbose)
    message(sprintf("%d genes passed tf-idf cut-off and %d soup quantile filter.  Taking the top %d.",nrow(mrks),nrow(filtPass),length(tgts)))
  #mrks = mrks[mrks$gene %in% tgts,]
  #tgts = head(mrks$gene,nMarks)
  if(length(tgts)==0){
    stop("No plausible marker genes found.  Reduce tfidfMin or soupQuantile")
  }
  if(length(tgts)<10){
    warning("Fewer than 10 marker genes found.  Consider reducing tfidfMin or soupQuantile")
  }
  ############################
  # Get estimates in clusters
  #Get which ones we'd use and where with canonical method
  tmp = as.list(tgts)
  names(tmp) = tgts
  ute = estimateNonExpressingCells(sc,tmp,maximumContamination=max(contaminationRange),FDR=rhoMaxFDR)
# estimateNonExpressingCells
# To decide if a cell is genuinely expressing a set of genes, a Poisson test is used. 
# This tests whether the observed expression is greater than {maximumContamination} times the expected number of counts for a set of genes, if the cell were assumed to be derived wholly from the background.  

  m = rownames(sc$metaData)[match(rownames(ssc$metaData),sc$metaData$clusters)]
  ute = t(ute[m,,drop=FALSE])
  colnames(ute) = rownames(ssc$metaData)
  #Now calculate the observed and expected counts for each cluster for 
  expCnts = outer(ssc$soupProfile$est,ssc$metaData$nUMIs)
  rownames(expCnts) = rownames(ssc$soupProfile)
  colnames(expCnts) = rownames(ssc$metaData)
  expCnts = expCnts[tgts,,drop=FALSE]
  #And the observed ones
  obsCnts = ssc$toc[tgts,,drop=FALSE]
  #We're done, but record some extra data for fun and profit
  #Filter out the shite
  #Get the p-value for this being less than 1
  pp = ppois(obsCnts,expCnts*max(contaminationRange),lower.tail=TRUE)
  qq = p.adjust(pp,method='BH')
  qq = matrix(qq,nrow=nrow(pp),ncol=ncol(pp),dimnames=dimnames(pp))
  #Get the cluster level ratio
  rhos = obsCnts/expCnts
  #Index in range
  rhoIdx = t(apply(rhos,1,function(e) order(order(e))))
  #Make a data.frame with everything
  dd = data.frame(gene = rep(rownames(ute),ncol(ute)),
                  passNonExp = as.vector(ute),
                  rhoEst = as.vector(rhos),
                  rhoIdx = as.vector(rhoIdx),
                  obsCnt = as.vector(obsCnts),
                  expCnt = as.vector(expCnts),
                  isExpressedFDR = as.vector(qq)
                  )
  dd$geneIdx = match(dd$gene,mrks$gene)
  dd$tfidf = mrks$tfidf[dd$geneIdx]
  dd$soupIdx = match(dd$gene,rownames(soupProf))
  dd$soupExp = soupProf$est[dd$soupIdx]
  dd$useEst = #dd$obsCnt >= minCnts & 
    #dd$isExpressedFDR < rhoMaxFDR & 
    #dd$rhoIdx <= min(clustPerGene,floor(ncol(rhoIdx)*maxClustFrac)) & 
    dd$passNonExp
  #The logic of piling up desity around the true value gets wonky if the number of estimates is low
  if(sum(dd$useEst)<10)
    warning("Fewer than 10 independent estimates, rho estimation is likely to be unstable.  Consider reducing tfidfMin or increasing SoupMin.")
  if(verbose)
    message(sprintf("Using %d independent estimates of rho.",sum(dd$useEst)))
  #Now aggregate the posterior probabilities for the ones we're including
  p.L = function(x,alpha){if(x==0){0}else{qgamma(alpha,x)}}
  p.U = function(x,alpha){qgamma(1-alpha,x+1)}
  alpha=0.95
  alpha=(1-alpha)/2
  dd$rhoHigh=sapply(seq(nrow(dd)),function(e) p.U(dd$obsCnt[e],alpha)/dd$expCnt[e])
  dd$rhoLow=sapply(seq(nrow(dd)),function(e) p.L(dd$obsCnt[e],alpha)/dd$expCnt[e])
  rhoProbes=seq(0,1,.001)
  #Using 95% confidence intervals
  #tmp = sapply(rhoProbes,function(e) {w=which(dd$useEst & dd$tfidf<1.5);sum(e>=dd$rhoLow[w] & e<=dd$rhoHigh[w])/length(w)})
  #Do a posterior estimation instead.  Use gamma prior defined by mode (priorRho) and standard deviation (priorRhoStdDev), which yields a posterior distribution for gamma of the form dgamma(rho,obsCnt+k,scale=theta/(1+theta*expCnts)). Where k and theta are the parameters for prior distribution derived using the above constraints.
  v2 = (priorRhoStdDev/priorRho)**2
  k = 1 +v2**-2/2*(1+sqrt(1+4*v2))
  theta = priorRho/(k-1)
  tmp = sapply(rhoProbes,function(e) {
                 tmp = dd[dd$useEst,]
                 mean(dgamma(e,k+tmp$obsCnt,scale=theta/(1+theta*tmp$expCnt)))
                  })
  #Calculate prior curve
  xx=dgamma(rhoProbes,k,scale=theta)
  #Get estimates
  w = which(rhoProbes>=contaminationRange[1] & rhoProbes<=contaminationRange[2])
  rhoEst = (rhoProbes[w])[which.max(tmp[w])]
  rhoFWHM = range((rhoProbes[w])[which(tmp[w]>=(max(tmp[w])/2))])
  contEst = rhoEst
  if(verbose)
    message(sprintf("Estimated global rho of %.2f",rhoEst))
  ##I think the best way to do this is based on the density.
  #tmp = density(dd$rhoEst[dd$useEst],...)
  #contEst = tmp$x[which.max(tmp$y)]
  if(doPlot){
    plot(rhoProbes,tmp,'l',
         xlim=c(0,1),
         ylim=c(0,max(c(xx,tmp))),
         frame.plot=FALSE,
         xlab='Contamination Fraction',
         ylab='Probability Density')
    #Add prior
    lines(rhoProbes,xx,lty=2)
    abline(v=rhoProbes[which.max(tmp)],col='red')
    legend(x='topright',
           legend=c(sprintf('prior rho %g(+/-%g)',priorRho,priorRhoStdDev),
                    sprintf('post rho %g(%g,%g)',rhoEst,rhoFWHM[1],rhoFWHM[2]),
                    'rho max'),
           lty=c(2,1,1),
           col=c('black','black','red'),
           bty='n')
    #plot(0,
    #     xlim=c(0,1),
    #     ylim=c(0,max(tmp$y)),
    #     type='n',
    #     frame.plot=FALSE,
    #     xlab='Contamination Fraction',
    #     ylab='Density'
    #     )
    #lines(tmp$x,tmp$y)
    #abline(v=contEst,col='red')
  }
  sc$fit = list(dd=dd,
                priorRho=priorRho,
                priorRhoStdDev=priorRhoStdDev,
                posterior = tmp,
                rhoEst = rhoEst,
                rhoFWHM = rhoFWHM
                )
  #Set the contamination fraction
  sc = setContaminationFraction(sc,contEst,forceAccept=forceAccept)
  return(sc)
}



# Estimate contamination rate for each cluster separately 
for(i in 1:ncol(ute)){
  dd = data.frame(gene = rownames(ute),
                  passNonExp = ute[,i],
                  rhoEst = rhos[,i],
                  rhoIdx = rhoIdx[,i],
                  obsCnt = obsCnts[,i],
                  expCnt = expCnts[,i],
                  isExpressedFDR = qq[,i]
                  )
  dd$geneIdx = match(dd$gene,mrks$gene)
  dd$tfidf = mrks$tfidf[dd$geneIdx]
  dd$soupIdx = match(dd$gene,rownames(soupProf))
  dd$soupExp = soupProf$est[dd$soupIdx]
  dd$useEst = #dd$obsCnt >= minCnts & 
    #dd$isExpressedFDR < rhoMaxFDR & 
    #dd$rhoIdx <= min(clustPerGene,floor(ncol(rhoIdx)*maxClustFrac)) & 
    dd$passNonExp
  #The logic of piling up desity around the true value gets wonky if the number of estimates is low
  if(sum(dd$useEst)<10)
    warning("Fewer than 10 independent estimates, rho estimation is likely to be unstable.  Consider reducing tfidfMin or increasing SoupMin.")
  if(verbose)
    message(sprintf("Using %d independent estimates of rho.",sum(dd$useEst)))
  #Now aggregate the posterior probabilities for the ones we're including
  p.L = function(x,alpha){if(x==0){0}else{qgamma(alpha,x)}}
  p.U = function(x,alpha){qgamma(1-alpha,x+1)}
  alpha=0.95
  alpha=(1-alpha)/2
  dd$rhoHigh=sapply(seq(nrow(dd)),function(e) p.U(dd$obsCnt[e],alpha)/dd$expCnt[e])
  dd$rhoLow=sapply(seq(nrow(dd)),function(e) p.L(dd$obsCnt[e],alpha)/dd$expCnt[e])
  rhoProbes=seq(0,1,.001)
  #Using 95% confidence intervals
  #tmp = sapply(rhoProbes,function(e) {w=which(dd$useEst & dd$tfidf<1.5);sum(e>=dd$rhoLow[w] & e<=dd$rhoHigh[w])/length(w)})
  #Do a posterior estimation instead.  Use gamma prior defined by mode (priorRho) and standard deviation (priorRhoStdDev), which yields a posterior distribution for gamma of the form dgamma(rho,obsCnt+k,scale=theta/(1+theta*expCnts)). Where k and theta are the parameters for prior distribution derived using the above constraints.
  v2 = (priorRhoStdDev/priorRho)**2
  k = 1 +v2**-2/2*(1+sqrt(1+4*v2))
  theta = priorRho/(k-1)
  tmp = sapply(rhoProbes,function(e) {
                 tmp = dd[dd$useEst,]
                 mean(dgamma(e,k+tmp$obsCnt,scale=theta/(1+theta*tmp$expCnt)))
                  })
  #Calculate prior curve
  xx=dgamma(rhoProbes,k,scale=theta)
  #Get estimates
  w = which(rhoProbes>=contaminationRange[1] & rhoProbes<=contaminationRange[2])
  rhoEst = (rhoProbes[w])[which.max(tmp[w])]
  rhoFWHM = range((rhoProbes[w])[which(tmp[w]>=(max(tmp[w])/2))])
  contEst = rhoEst
  print(rhoEst)
}

Using 19 independent estimates of rho.
[1] 0.001
Using 93 independent estimates of rho.
[1] 0.002
Using 97 independent estimates of rho.
[1] 0.001
Using 97 independent estimates of rho.
[1] 0.001
Using 95 independent estimates of rho.
[1] 0.001
Using 98 independent estimates of rho.
[1] 0.001
Using 94 independent estimates of rho.
[1] 0.003
Using 97 independent estimates of rho.
[1] 0.004
Using 96 independent estimates of rho.
[1] 0.002
Using 96 independent estimates of rho.
[1] 0.001
Using 94 independent estimates of rho.
[1] 0.001
Using 93 independent estimates of rho.
[1] 0.001
Using 99 independent estimates of rho.
[1] 0.005
Using 97 independent estimates of rho.
[1] 0.003
Using 99 independent estimates of rho.
[1] 0.001
Using 97 independent estimates of rho.
[1] 0.001
Using 99 independent estimates of rho.
[1] 0.001
Using 91 independent estimates of rho.
[1] 0.002
Using 92 independent estimates of rho.
[1] 0.001


plotSoupCorrelation = function(sc){
  if(!is(sc,'SoupChannel'))
    stop("sc not a valid SoupChannel object.")
  #Calculate the cell profile
  cellProfile = rowSums(sc$toc)
  cellProfile = (cellProfile/sum(cellProfile))
  df = data.frame(cellProfile,soupProfile=sc$soupProfile$est)
  gg = ggplot(df,aes(log10(cellProfile),log10(soupProfile))) +
    geom_point(alpha=1/3) +
    geom_abline(intercept=0,slope=1) +
    ylab('log10(Soup Expression)')+
    xlab('log10(Aggregate cell Expression)')
  return(gg)
}

tiff("plot/SoupX_Fimbria_3kSoupbc_CorSoupCluster.tiff",res=300,height=800,width=800)
plotSoupCorrelation(sc)
dev.off()



# rank correlation between total expression in each of 19 cluster and soup 
dim(ssc$soupProfile)
#[1] 33538     2
dim(ssc$toc)
#[1] 33538    19
all=cbind(ssc$toc,ssc$soupProfile[,2])
# using all detected genes
rall=cor(all,method="sp") 
# using top 10% soup genes that are markers for 19 clusters (N=770)
rsoupgenemarker=cor(all[filtPass$gene,],method="sp")
# using top 100 genes belonging to the above list (N=100)
rtgts=cor(all[tgts,],method="sp")
print(cbind(rall[,20],rsoupgenemarker[,20],rtgts[,20]))
1-4_1 0.012634069 -3.098486e-02 -0.100537423
1-4_2 0.011929361 -2.777706e-02 -0.163538307
1-4_3 0.016947333  1.339369e-02 -0.150478417
1-4_4 0.017505999  1.229255e-02 -0.134409523
1-4_5 0.019971793  2.402282e-02 -0.115299010
1-4_6 0.014117913  7.970833e-03 -0.065116340
5     0.017169005 -1.110640e-02  0.006563424
6     0.010309935  1.289082e-02 -0.047654969
7     0.011670425  2.145846e-02 -0.144460067
8     0.013605573  2.531106e-02  0.062042288
9     0.015600300 -3.277332e-02 -0.010558272
10_1  0.013753098 -9.995585e-05  0.112054811
10_2  0.004899765 -3.829064e-02 -0.133965179
10_3  0.011871510 -2.737206e-02  0.013241452
11    0.015527447 -8.436088e-02  0.092206931
12    0.014144066 -7.664134e-02 -0.100256251
13    0.009993825 -5.254525e-02 -0.057890738
14    0.018849799 -4.634259e-02 -0.072900296
15    0.019786544 -7.618746e-03  0.038381733
soup  1.000000000  1.000000e+00  1.000000000


# rank correlation between 19 cluster centroids and soup centroid
  clustercentroids = do.call(cbind,lapply(s,function(e) apply(sc$toc[,e,drop=FALSE],1,mean)))
  soupcentroid=ssc$soupProfile[,2]/2911
all=cbind(clustercentroids,soupcentroid)
# using all detected genes
rall=cor(all,method="sp") 
# using top 10% soup genes that are markers for 19 clusters (N=770)
rsoupgenemarker=cor(all[filtPass$gene,],method="sp")
# using top 100 genes belonging to the above list (N=100)
rtgts=cor(all[tgts,],method="sp")
print(cbind(rall[,20],rsoupgenemarker[,20],rtgts[,20]))


# rank correlation between all cells in each of the 19 clusters and the soup and the 19 cluster centroids and soup centroid
library(SoupX)
i=1
load(file=paste0("Soup_",names[i],"_3kSoup.Robj"))

dataDir=paste0("/nfs/turbo/umms-hammou/10x/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/")
tod = Seurat::Read10X(file.path(dataDir,'raw_feature_bc_matrix'))
colnames(tod)=paste(names[i],colnames(tod),sep="_")
nDropUMIs=sort(sc$nDropUMIs,decreasing=TRUE)
soup=names(nDropUMIs)[which(nDropUMIs<100)][1:2911]
length(soup)

### normalized all cells 
  clusterinfo=sc$metaData[order(sc$metaData$clusters),]
  clusters=clusterinfo$clusters
  names(clusters)=rownames(clusterinfo)
  clustercells=sc$toc[,names(clusters)]
  which(clustercells !=  tod[rownames(sc$toc),names(clusters)]) # integer(0)
  soupcells=tod[rownames(sc$toc),soup]
  which(rownames(clustercells) != rownames(soupcells)) # integer(0)
allcells=cbind(clustercells,soupcells)
  soupclusters=rep("Soup",length(soup))
  names(soupclusters)=soup
allclusters=c(as.character(clusters),soupclusters)
names(allclusters)=c(names(clusters),soup)
allclusters=factor(allclusters,levels=c(levels(clusters),"Soup"))
table(allclusters)
dge1 <- CreateSeuratObject(counts = allcells, project = "Fimbria3ksoup", min.cells = 0, min.features = 0)
### Normalize data
dge1<-NormalizeData(dge1)
dge1$clusters <- allclusters
Idents(dge1) <- allclusters
### centroids
centroids1=AverageExpression(dge1)
centroids=log(centroids1$RNA+1)
  which(rownames(centroids) != rownames(allcells)) # integer(0)
### combining log-transformed normalized expression for all cells and centroids
all=cbind(centroids,GetAssayData(dge1))
save(all,file="plot/Fimbria3ksoup_normalizedexp_allcellscentroids.Robj")
# using HVG for cluster centroids
dge2<-subset(dge1,clusters != "Soup")
dge2<-FindVariableFeatures(dge2)
hvg<-VariableFeatures(dge2)
length(hvg) # 2000
allhvg=as.matrix(all[hvg,])
rhvg=cor(as.matrix(allhvg),method="sp") 
data=rhvg[-c(1:20),1:20]
dim(data) # 4772 20
data.std=(data-apply(data,1,mean))/apply(data,1,sd)
write.table(data,"plot/rho_allcells_centroids_clustersHVG.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(data.std,"plot/rho_allcells_centroids_clustersHVG_std.txt",row.names=T,col.names=T,quote=F,sep="\t")
# using top 10% soup genes that are markers for 19 clusters (N=770)
allmarker=as.matrix(all[filtPass$gene,])
dim(allmarker) # 770 4792
rsoupgenemarker=cor(allmarker,method="sp")
data=rsoupgenemarker[-c(1:20),1:20]
dim(data) # 4772 20
data.std=(data-apply(data,1,mean))/apply(data,1,sd)
write.table(data,"plot/rho_allcells_centroids_770soupmarkers.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(data.std,"plot/rho_allcells_centroids_770soupmarkers_std.txt",row.names=T,col.names=T,quote=F,sep="\t")

library(RColorBrewer)
tmp=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
myBrewerPalette=c(gg_color_hue(6),tmp[5:9],gg_color_hue(3),tmp[11:15])

data.use=data.std
levels=colnames(data.std)
colsep.use=cumsum(table(levels)[levels])
col.lab=rep("",length(levels))
col.lab=levels
row.lab=NULL
ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

redblue100<-rgb(read.table(paste0(home,"data_DGE/redblue100.txt"),sep='\t',row.names=1,header=T))
col.use=redblue100

library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

jpeg(file=paste0("plot/rho_HVG_allcellsVsCentroids_std.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0("plot/rho_HVG_allcellsVsCentroids_std2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()

jpeg(file=paste0("plot/rho_SoupMarkers_allcellsVsCentroids_std.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0("plot/rho_SoupMarkers_allcellsVsCentroids_std2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
