# 3.12.2021 scVelo on 3 healthy FTs by Qianyi
# run scVelo
# Related to Figure 5


module add python/3.7.4
python

import os 
import pandas as pd
import numpy as np
os.chdir("/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/10xFallopialTube/")
import scvelo as scv
adata = scv.read("FT124.h5ad")       # global clusters
adata = scv.read("FT124_C2-6.h5ad")  # secretory subclusters
pca_cord = pd.read_csv("plot/pcaC2-6.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata = scv.read("FT124_C7-10.h5ad") # stromal zoomed-in view (kept global clusters)
pca_cord = pd.read_csv("plot/pcaC7-10.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata = scv.read("FT124_C1-6.h5ad") # stromal zoomed-in view (kept global clusters)
pca_cord = pd.read_csv("plot/pcaC1-6.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata = scv.read("FT124_C1-8.h5ad") # stromal zoomed-in view (kept global clusters)
pca_cord = pd.read_csv("plot/pcaC1-8.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata = scv.read("FT124_C1-2_2.h5ad") # stromal zoomed-in view (kept global clusters)
pca_cord = pd.read_csv("plot/pcaC1-2_2.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata = scv.read("FT124_C127.h5ad") # stromal zoomed-in view (kept global clusters)
pca_cord = pd.read_csv("plot/pcaC127.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata = scv.read("FT124_C12678.h5ad") # stromal zoomed-in view (kept global clusters)
pca_cord = pd.read_csv("plot/pcaC12678.csv").rename(columns = {'Unnamed: 0':'Cell ID'})
adata
## AnnData object with n_obs × n_vars = 6667 × 24421
##     obs: 'orig.ident', 'nCount_spliced', 'nFeature_spliced', 'nCount_unspliced', 'nFeature_unspliced', 'nCount_ambiguous', 'nFeature_ambiguous', 'nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT', 'SCT_snn_res.0.8', 'seurat_clusters'
##     var: 'features', 'ambiguous_features', 'spliced_features', 'unspliced_features'
##     obsm: 'X_umap'
##     layers: 'ambiguous', 'spliced', 'unspliced'

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
adata.obsm['X_pca']
adata.obsm['X_pca']=pd.DataFrame(adata.obs.index).rename(columns = {0:'Cell ID'}).merge(pca_cord,on="Cell ID").iloc[:,1:].values


adata.obs['seurat_clusters']
# Global
cell_clusters = pd.read_csv("plot/clusters.csv", dtype=str)[['CellID','Seurat_clusters']]
adata.obs['clusters']=pd.DataFrame(adata.obs.index).rename(columns = {0:'CellID'}).merge(cell_clusters,on="CellID").set_index('CellID').iloc[:,0]
adata.uns['clusters_colors'] = np.asarray(['#A6CEE3','grey','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#CAB2D6','#6A3D9A','#FFFF99','#B15928']) # global C2-6 together, No C13
adata.uns['clusters_colors'] = np.asarray(['#A6CEE3','#FB9A99','#E31A1C','#FDBF6F','#CAB2D6','#6A3D9A','#FFFF99','#B15928','grey','#1F78B4','#B2DF8A','#33A02C']) # global C2-6 together, No C13
# Secretory subset
### if it is a category:
adata.obs['clusters']=adata.obs['seurat_clusters'].cat.codes + 1
### if it is int
adata.obs['clusters']=adata.obs['seurat_clusters'] + 1
adata.uns['clusters_colors'] = np.asarray(['#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F']) # secretory subclusters C2-6
# Stromal subset C7-10
adata.obs['clusters']=adata.obs['seurat_clusters'] + 1
adata.uns['clusters_colors'] = np.asarray(['#1F78B4','#B2DF8A','#33A02C','#FB9A99']) 
# Epithelial subset C1-6
cell_clusters = pd.read_csv("plot/clustersC1-6.csv", dtype=str)[['CellID','Seurat_clusters']]
adata.obs['clusters']=pd.DataFrame(adata.obs.index).rename(columns = {0:'CellID'}).merge(cell_clusters,on="CellID").set_index('CellID').iloc[:,0]
adata.uns['clusters_colors'] = np.asarray(['#A6CEE3','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F']) 
# Epithelial+new stromal subset C1-8
cell_clusters = pd.read_csv("plot/clustersC1-8.csv", dtype=str)[['CellID','Seurat_clusters']]
adata.obs['clusters']=pd.DataFrame(adata.obs.index).rename(columns = {0:'CellID'}).merge(cell_clusters,on="CellID").set_index('CellID').iloc[:,0]
adata.uns['clusters_colors'] = np.asarray(['#A6CEE3','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#CAB2D6','#6A3D9A','#FFFF99','#B15928']) 
# Ciliated subclusters + secretory 2_2
cell_clusters = pd.read_csv("plot/clustersC1-2_2.csv", dtype=str)[['CellID','Seurat_clusters']]
adata.obs['clusters']=pd.DataFrame(adata.obs.index).rename(columns = {0:'CellID'}).merge(cell_clusters,on="CellID").set_index('CellID').iloc[:,0]
adata.uns['clusters_colors'] = np.asarray(['#F8766D','#7CAE00','#00BFC4','#C77CFF','#FC8D62']) 
# Epithelial subset C1, 2 + endothelial C7
cell_clusters = pd.read_csv("plot/clustersC127.csv", dtype=str)[['CellID','Seurat_clusters']]
adata.obs['clusters']=pd.DataFrame(adata.obs.index).rename(columns = {0:'CellID'}).merge(cell_clusters,on="CellID").set_index('CellID').iloc[:,0]
adata.uns['clusters_colors'] = np.asarray(['#A6CEE3','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F', '#E31A1C']) 
# Epithelial subset C1, 2 + endothelial C678
cell_clusters = pd.read_csv("plot/clustersC12678.csv", dtype=str)[['CellID','Seurat_clusters']]
adata.obs['clusters']=pd.DataFrame(adata.obs.index).rename(columns = {0:'CellID'}).merge(cell_clusters,on="CellID").set_index('CellID').iloc[:,0]
adata.uns['clusters_colors'] = np.asarray(['#A6CEE3','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#FB9A99','#E31A1C','#FDBF6F']) 

adata.obs['clusters'].value_counts()

# default, very fast - used this
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# dynamic mode, slow
scv.tl.recover_dynamics(adata)  
scv.tl.velocity(adata,mode='dynamical') 
scv.tl.velocity_graph(adata)

# Global
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="Global.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="GlobalPCA.png")
# Secretory subset
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="fast1C2-6.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC2-6.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="fastC2-6PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC2-6tSNE.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC2-6.pdf")
# Stromal subset
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="fast1C7-10.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC7-10.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="fastC7-10PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC7-10.pdf")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC7-10tSNE.png")
# Epithelial subset
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="fast1C1-6.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="C1-6.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="C1-6PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC1-6tSNE.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", save="fastC1-6.1.pdf")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC1-6.pdf")
# Epithelial+newStromal subset (removed immune cells)
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="fastC1-8.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="C1-8.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="fastC1-8PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC1-8tSNE.png")
# Ciliated subclusters + secretory 2_2
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="fastC1-2_2.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="C1-2_2.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="fastC1-2_2PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="pca",legend_loc='right margin', color="clusters", dpi=300,save="C1-2_2PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC1-2_2tSNE.png")
# Epithelial subset C1, 2 + endothelial C7
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="C127_2.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC127_2.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="fastC127_2PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="pca",legend_loc='right margin', color="clusters", dpi=300,save="C127_2PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC127_2tSNE.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne",legend_loc='right margin', color="clusters", dpi=300,save="C127_2tSNE.png")
# Epithelial subset C1, 2 + endothelial C678
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc='right margin',color="clusters", dpi=300,save="C12678_2.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters", dpi=300,save="fastC12678_2.png")
scv.pl.velocity_embedding_stream(adata, basis="pca", color="clusters", dpi=300,save="fastC12678_2PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="pca",legend_loc='right margin', color="clusters", dpi=300,save="C12678_2PCA.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="clusters", dpi=300,save="fastC12678_2tSNE.png")
scv.pl.velocity_embedding_stream(adata, basis="tsne",legend_loc='right margin', color="clusters", dpi=300,save="C12678_2tSNE.png")
# saved as Figure 5.


scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata,mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color="clusters")


scv.pl.velocity_embedding(adata, basis="umap", color="clusters", arrow_length=3, arrow_size=1, dpi=120,save="Global_cells.png")


scv.tl.recover_dynamics(adata)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", dpi=300,save="Global_latent.png")


top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color="clusters", n_convolve=100, dpi=300,save="Global_gene.png")
