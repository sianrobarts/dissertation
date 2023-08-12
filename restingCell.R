library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("restingCells_CD4only_HVGs_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
resting <- LoadH5Seurat("restingCells_CD4only_HVGs_processed.h5seurat")
resting
