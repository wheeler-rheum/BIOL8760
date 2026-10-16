# Stage 4. Construct Uniform Manifold Approximation and Projection (UMAP) and identify per-cluster marker genes
library(dplyr)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')

## Stage 4.i.	Construct UMAP for data reduction based on first 14 dimensions
# Construct UMAP using the first 14 principal components
RA1_seurat <- RunUMAP(RA1_seurat, dims = 1:14)

# Visualize the UMAP, colored by cluster
DimPlot(RA1_seurat, reduction = "umap", label = TRUE)

# Repeat for RA2
RA2_seurat <- RunUMAP(RA2_seurat, dims = 1:14)
DimPlot(RA2_seurat, reduction = "umap", label = TRUE)

# Repeat for RA4
RA4_seurat <- RunUMAP(RA4_seurat, dims = 1:14)
DimPlot(RA4_seurat, reduction = "umap", label = TRUE)

# Repeat for RA5
RA5_seurat <- RunUMAP(RA5_seurat, dims = 1:14)
DimPlot(RA5_seurat, reduction = "umap", label = TRUE)

# Repeat for RA6
RA6_seurat <- RunUMAP(RA6_seurat, dims = 1:14)
DimPlot(RA6_seurat, reduction = "umap", label = TRUE)

## Stage 4.ii.	Identify marker genes for every cluster
# Identify marker genes for every cluster
RA1_seurat <- JoinLayers(RA1_seurat)
RA1_markers <- FindAllMarkers(RA1_seurat)

# Repeat for RA2
RA2_seurat <- JoinLayers(RA2_seurat)
RA2_markers <- FindAllMarkers(RA2_seurat)

# Repeat for RA4
RA4_seurat <- JoinLayers(RA4_seurat)
RA4_markers <- FindAllMarkers(RA4_seurat)

# Repeat for RA5
RA5_seurat <- JoinLayers(RA5_seurat)
RA5_markers <- FindAllMarkers(RA5_seurat)

# Repeat for RA6
RA6_seurat <- JoinLayers(RA6_seurat)
RA6_markers <- FindAllMarkers(RA6_seurat)
