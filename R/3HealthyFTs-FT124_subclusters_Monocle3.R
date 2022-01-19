# 7.2.2021 Monocle3 for FT124 subclusters by Qianyi



R

library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

library(SeuratDisk)
library(SeuratWrappers)
library(monocle3)


gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette1=gg_color_hue(5)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(7,"Set2"),brewer.pal(12,"Paired")[c(10,12)])

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

### using original 17 clusters
subsets=list(2:6,7:10,1:6)
subsetsname=c("2-6","7-10","1-6")
# 2:6  secretory subset   -> 6 subclusters
# 7:10 old stromal subset -> 5 subclusters - very similar as global clusters
# 1:6  ciliated cluster 1 + 6 secretory subclusters

### after renaming the 12 clusters to 1:12
subsets=list(1:8,3:8,1,c(1:2,7),c(1:2,6:8))
subsetsname=c("1-8","3-8",1,"127","12678")
# 1:8   epithelial + new stromal, removing immune cells
# 3:8   new stromal including myofibroblast and endothelial
# 1     ciliated only -> 4 ciliated subclusters
# 127   epithelial   + endothelial cluster 7
# 12678 epithelial + pericyte 6 + endothelial 7&8


subsets=list(c(1,"2_2"))
subsetsname=c("1-2_2")


subsets=list(1:6,c(1,"2_2"),c(1,2,7))
subsetsname=c("1-6","1-2_2","127")
resSeurat=c("ft124CCA3all2","ordered","ft124CCA3all22")
myBrewerPalette=list(c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")),
c(gg_color_hue(4),brewer.pal(7,"Set2")[2]),
c(brewer.pal(12,"Paired")[1],brewer.pal(7,"Set2")[1:6],brewer.pal(12,"Paired")[6])
  )

# in the future
### load subclustering for CCAallgenes-corrected 3 healthy FTs - 17 re-ordered clusters
### CCAallgenes -> global cluster -> extract subset -> use integration of the subset to select HVG, use global corrected gene subsetsnameression matrix as integrated normalized data, rescale for the subset -> PCA and subcluster
dgelist=list()
for(indiv in 1:length(subsets)){
cc=subsets[[indiv]]
ccname=subsetsname[indiv]
dgefile=paste0("/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/10xFallopialTube/plot/C",ccname,"_")
load(file=paste0(all,"_C",ccname,".Robj"))
dgelist[[indiv]]=dge
}
dge
indiv=1


### load Monocle3 of indiviudal subsetsnameeriment
cdslist=list()
for(indiv in 1:length(subsetsname)){
  load(file=paste0(subsetsname[indiv],"_monocle3.Robj"))
  cdslist[[indiv]]=cds
}


library(monocle3)
library(ggplot2)
library(dplyr)



for(indiv in 1:length(subsetsname)){

indiv=1
indiv=2
indiv=3
indiv=4

cds=cdslist[[indiv]]
levels=as.character(sort(unique(pData(cds)[resSeurat[indiv]]))[,1])



dge=dgelist[[indiv]]

DefaultAssay(dge) <- "integrated"

cds <- as.cell_data_set(dge)
reducedDims(cds)
reducedDims(cds)$UMAP <- reducedDims(cds)$umap
colnames(reducedDims(cds)$UMAP) <- NULL
#reducedDims(cds)$PCA <- reducedDims(cds)$pca
if(indiv==2){
reducedDims(cds)$UMAP <- reducedDims(cds)$pca
colnames(reducedDims(cds)$UMAP) <- NULL
}

levels=as.character(sort(unique(pData(cds)[resSeurat[indiv]]))[,1])


#Cluster your cells
## Step 3: Reduce the dimensions using UMAP
pdf(paste0("plot/C",subsetsname[indiv],"_umap.pdf"),height=2.5,width=3.5)
# visualize clusters from Seurat
plot_cells(cds, color_cells_by=resSeurat[indiv],
           label_cell_groups = FALSE,
           show_trajectory_graph=FALSE,
           graph_label_size=3) + scale_color_manual(values=myBrewerPalette[[indiv]])
dev.off()
knownmarkers=c("Lyz2","Vim","Pdgfra","Tcf21","Ly6a","Hsd3b1","Cyp17a1","Cyp11a1","Star","Insl3","Thra","Thrb","Acta2","Myh11","Clu","Sox9","Kit","Mki67","Nr2f2","Zbtb16","Stra8","Prm1")
jpeg(paste0("plot/C",subsetsname[indiv],"_umap_knownmarkers.jpeg"),res=300,height=1600,width=2000)
plot_cells(cds,
           genes=knownmarkers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

res=c(NULL or 1e-4?,8e-4)
res=c(4e-6,3e-6)
## Step 4: Cluster the cells
if(res[indiv]==0){
	cds = cluster_cells(cds)  # default resolution=NULL
} else{
	cds = cluster_cells(cds, resolution=res[indiv])
}
pdf(paste0("plot/C",subsetsname[indiv],"_cluster.pdf"),height=2.5,width=2.5)
# visualize clusters from Monocle3
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster",
           show_trajectory_graph=FALSE,
           group_label_size=5)
# visualize partitions - more separated groups of clusters
plot_cells(cds, color_cells_by="partition", group_cells_by="partition",
           show_trajectory_graph=FALSE,
           group_label_size=5) 
dev.off()


#Order cells in pseudotime along a trajectory
## Step 5: Learn a graph
cds <- learn_graph(cds)
pdf(paste0("plot/C",subsetsname[indiv],"_graph.pdf"),height=2.5,width=2.5)
plot_cells(cds,
           color_cells_by = resSeurat[indiv],
           label_groups_by_cluster=FALSE,
           group_label_size=0,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + scale_color_manual(values=myBrewerPalette[[indiv]])
plot_cells(cds,
           color_cells_by = resSeurat[indiv],
           label_groups_by_cluster=FALSE,
           group_label_size=3,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=2) + scale_color_manual(values=myBrewerPalette[[indiv]])
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           group_label_size=5,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           group_label_size=3,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=2)
dev.off()

## Step 6: Order cells
### manually choose the root nodes that indicates time point 0 
cds <- order_cells(cds)
pdf(paste0("plot/C",subsetsname[indiv],"_pseudotime_manual.pdf"),height=2.5,width=3.5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3)
dev.off()

save(cds,file=paste0(subsetsname[indiv],"_monocle3.Robj"))

### programmatically choose the root nodes that indicates time 0 
# a helper function to identify the root principal points:
if(indiv<3){
	start=levels(colData(cds)[, "time2"])[2]
} else {
	start=min(colData(cds)[, "time2"])
}

if(indiv==1){
  start="2-6_2"
}
if(indiv==2){
  start="2_2"
}

get_earliest_principal_node <- function(cds, time_bin=start){
  cell_ids <- which(colData(cds)[, resSeurat[indiv]] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pdf(paste0("plot/C",subsetsname[indiv],"_pseudotime_program.pdf"),height=2.5,width=3.5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3)
dev.off()

save(cds,file=paste0(subsetsname[indiv],"_monocle3.Robj"))
save(cds,file=paste0(subsetsname[indiv],"_monocle3_1partition.Robj"))


### plot genes in pseudotime
knownmarkers=c("Pdgfra","Tcf21","Ly6a","Hsd3b1","Cyp17a1","Cyp11a1","Star","Insl3","Thra","Thrb","Acta2","Myh11","Clu","Kit","Mki67","Nr2f2")
knownmarkers_cds <- cds[rowData(cds)$gene_short_name %in% knownmarkers,
                       ]
jpeg(paste0("plot/C",subsetsname[indiv],"_pseudotime_knownmarkers.jpeg"),res=300,height=1200,width=2500)
plot_genes_in_pseudotime(knownmarkers_cds,
                         color_cells_by=resSeurat[indiv],
                         min_subsetsnamer=0.5,ncol=4) + scale_color_manual(values=myBrewerPalette[[indiv]])
dev.off()

### visualize more markers - 3/13/2020
knownmarkers=read.table("3.12.20plotmarkers.txt",stringsAsFactors=F)[,1]
length(knownmarkers) # 37
knownmarkers[which(!(knownmarkers %in% rowData(cds)$gene_short_name))]
knownmarkers_cds <- cds[rowData(cds)$gene_short_name %in% knownmarkers,
                       ]
jpeg(paste0("plot/C",subsetsname[indiv],"_pseudotime_more_markers.jpeg"),res=300,height=1800,width=3750)
plot_genes_in_pseudotime(knownmarkers_cds,
                         color_cells_by=resSeurat[indiv],
                         min_subsetsnamer=0.5,ncol=6) + scale_color_manual(values=myBrewerPalette[[indiv]])
dev.off()
jpeg(paste0("plot/C",subsetsname[indiv],"_pseudotime_more_markers2.jpeg"),res=300,height=1500,width=3125)
plot_genes_in_pseudotime(knownmarkers_cds,cell_size=0.5,
                         color_cells_by=resSeurat[indiv],
                         min_subsetsnamer=0.5,ncol=6) + scale_color_manual(values=myBrewerPalette[[indiv]])
dev.off()




### Working with 3D trajectories
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d) # manual
#cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")



}

