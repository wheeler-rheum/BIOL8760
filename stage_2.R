### Stage 2: Selection of 2000 most variable features and scaling to construct principal component analysis (PCA)

## Step 2.i. Select the 2000 most variable features
RA1_seurat <- FindVariableFeatures(
  RA1_seurat,
  selection.method = "vst",
  nfeatures = 2000
)
RA2_seurat <- FindVariableFeatures(
  RA2_seurat,
  selection.method = "vst",
  nfeatures = 2000
)
RA4_seurat <- FindVariableFeatures(
  RA4_seurat,
  selection.method = "vst",
  nfeatures = 2000
)
RA5_seurat <- FindVariableFeatures(
  RA5_seurat,
  selection.method = "vst",
  nfeatures = 2000
)
RA6_seurat <- FindVariableFeatures(
  RA6_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

## Step 2.ii.	Scale the data
RA1_seurat <- ScaleData(RA1_seurat)
RA2_seurat <- ScaleData(RA2_seurat)
RA4_seurat <- ScaleData(RA4_seurat)
RA5_seurat <- ScaleData(RA5_seurat)
RA6_seurat <- ScaleData(RA6_seurat)

## Step 2.iii.	Construct PCA
RA1_seurat <- RunPCA(RA1_seurat, features = VariableFeatures(object = RA1_seurat))
RA2_seurat <- RunPCA(RA2_seurat, features = VariableFeatures(object = RA2_seurat))
RA4_seurat <- RunPCA(RA4_seurat, features = VariableFeatures(object = RA4_seurat))
RA5_seurat <- RunPCA(RA5_seurat, features = VariableFeatures(object = RA5_seurat))
RA6_seurat <- RunPCA(RA6_seurat, features = VariableFeatures(object = RA6_seurat))
