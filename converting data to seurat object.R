library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("restingCells_CD4only_HVGs_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc <- LoadH5Seurat("restingCells_CD4only_HVGs_processed.h5seurat")
pbmc

pbmc <- NormalizeData(pbmc) #normalising the data according to seurat vignete

lookatmatrix <- pbmc

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #identifying key genes

# Scale the data
pbmc <- ScaleData(pbmc)
# The expression values of each gene in each cell are now scaled to have zero mean and unit variance

# visualize the top 10 variable features
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0, )
plot1 + plot2

#perfoming linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) 
VizDimLoadings(pbmc, dims = 1:5, reduction = "pca") # visualize gene along the first five principal components and how they contribute +ve/-ve to PC.  

# show the expression distribution of the top 10 variable features across different cell populations or groups.
VlnPlot( pbmc, features = top10, cols = NULL, pt.size = NULL,
  idents = NULL,
  sort = FALSE,
  group.by = NULL,
  split.by = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  slot = "data",
  split.plot = FALSE,
  stack = FALSE,
  combine = TRUE,
  fill.by = "feature",
  flip = FALSE,
  raster = NULL
)

# visualize cells in reduced dimensions (PCA space)
DimPlot(pbmc, reduction = "pca", raster=FALSE) #this gave the scatterlooking plot.
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) #Examine and visualize PCA results another way: lsit


# heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. 
DimHeatmap(pbmc, dims = 1:10, nfeatures = 30, cells = 500, reduction = "pca", balanced = TRUE) 
# Allows for nice visualization of sources of heterogeneity in the dataset.


#trying to run jackstraw plot and elbow plot
pbmc <- JackStraw(pbmc, num.replicate = 100) 
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20, cols = NULL, reduction = "pca", xmax = 0.1, ymax = 0.3) 

# Elbow plot used to visualize the explained variance by each principal component and helps in determining the optimal number of dimensions to retain for downstream analysis.
ElbowPlot(pbmc)

#cluster the cells 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#Run non-linear dimensional reduction (UMAP)
pbmc <- RunUMAP(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "umap")

cluster2.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0)
head(cluster2.markers, n = 10)

RunTSNE( pbmc, reduction = "pca", cells = NULL,
  dims = 1:5,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
)

RunTSNE(pbmc, reduction = "pca")









