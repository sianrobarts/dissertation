library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("restingCells_CD4only_HVGs_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc <- LoadH5Seurat("restingCells_CD4only_HVGs_processed.h5seurat")
pbmc

pbmc <- NormalizeData(pbmc) #normalising the data according to seurat vignete

lookatmatrix <- pbmc

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #identifying key genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#perfoming linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) 
VizDimLoadings(pbmc, dims = 1:5, reduction = "pca") #this gave the list of gene plot
DimPlot(pbmc, reduction = "pca", raster=F) #this gave the scatterlooking plot.
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) #Examine and visualize PCA results another way: lsit

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #was unable to run this...

#trying to run jackstraw plot
pbmc <- JackStraw(pbmc, num.replicate = 100) #run later babe. Update, did run later now doesnt work idk why
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(
  pbmc,
  dims = 1:20,
  cols = NULL,
  reduction = "pca",
  xmax = 0.1,
  ymax = 0.3
) 

ElbowPlot(pbmc)

#cluster the cells 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "umap")

cluster2.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0)
head(cluster2.markers, n = 5)











