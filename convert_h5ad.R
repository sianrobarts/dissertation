install.packages("renv") # installing {renv}
renv::init() # ask {renv} to create an environment 

# need to install the following R commands 
pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)

bioc_pkgs <- c(
    "SingleCellExperiment",
    "scater",
    "multtest"
)
# If you are using an {renv} environment
renv::install(pkgs)


# Install Bioconductor packages
BiocManager::install(bioc_pkgs, update = FALSE)


# Install the following Python packages
py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)
renv::snapshot()

# The {reticulate} approach 
library(Seurat)
library(SeuratData)
library(SeuratDisk)

head(adata$obs)
