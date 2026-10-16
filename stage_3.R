## Stage 3: Construct Shared Nearest Neighbor (SNN) graph based on first 14 PCs with Euclidean distance and cluster results

# Step 3.i.	Construct SNN Graph
RA1_seurat <- FindNeighbors(RA1_seurat, dims = 1:14)
RA2_seurat <- FindNeighbors(RA2_seurat, dims = 1:14)
RA4_seurat <- FindNeighbors(RA4_seurat, dims = 1:14)
RA5_seurat <- FindNeighbors(RA5_seurat, dims = 1:14)
RA6_seurat <- FindNeighbors(RA6_seurat, dims = 1:14)

# Step 3.ii.	Cluster on SNN graph using Louvain algorithm and resolution of 0.5
RA1_seurat <- FindClusters(
  RA1_seurat,
  resolution = 0.5,
  algorithm = 1
)
table(Idents(RA1_seurat))

RA2_seurat <- FindClusters(
  RA2_seurat,
  resolution = 0.5,
  algorithm = 1
)
table(Idents(RA2_seurat))

RA4_seurat <- FindClusters(
  RA4_seurat,
  resolution = 0.5,
  algorithm = 1
)
table(Idents(RA4_seurat))

RA5_seurat <- FindClusters(
  RA5_seurat,
  resolution = 0.5,
  algorithm = 1
)
table(Idents(RA5_seurat))

RA6_seurat <- FindClusters(
  RA6_seurat,
  resolution = 0.5,
  algorithm = 1
)
table(Idents(RA6_seurat))
