# Stage 5. Annotate clusters based on expression of top marker genes with Monaco immune reference dataset as reference

install.packages("curl")
BiocManager::install(c("SingleR", "celldex"))
library(curl)
library(SingleR)
library(celldex)
library(SummarizedExperiment)
library(ggplot2)

# Download the Monaco immune reference from celldex
monaco_ref <- MonacoImmuneData()

# Inspect the available labels (fine-grained and main cell types)
table(monaco_ref$label.main)

# Extract the expression matrix from the Seurat object
# JoinLayers() was run previously, so there is a single counts/data layer
RA1_expr_matrix <- GetAssayData(RA1_seurat, assay = "RNA", layer = "data")
RA2_expr_matrix <- GetAssayData(RA2_seurat, assay = "RNA", layer = "data")
RA4_expr_matrix <- GetAssayData(RA4_seurat, assay = "RNA", layer = "data")
RA5_expr_matrix <- GetAssayData(RA5_seurat, assay = "RNA", layer = "data")
RA6_expr_matrix <- GetAssayData(RA6_seurat, assay = "RNA", layer = "data")

# Run SingleR to annotate each cell against the Monaco reference
# Use label.main for broad cell types; switch to "label.fine" for subtypes
RA1_pred <- SingleR(
  test = RA1_expr_matrix,
  ref = monaco_ref,
  labels = monaco_ref$label.main
)
RA2_pred <- SingleR(
  test = RA2_expr_matrix,
  ref = monaco_ref,
  labels = monaco_ref$label.main
)
RA4_pred <- SingleR(
  test = RA4_expr_matrix,
  ref = monaco_ref,
  labels = monaco_ref$label.main
)
RA5_pred <- SingleR(
  test = RA5_expr_matrix,
  ref = monaco_ref,
  labels = monaco_ref$label.main
)
RA6_pred <- SingleR(
  test = RA6_expr_matrix,
  ref = monaco_ref,
  labels = monaco_ref$label.main
)

# Check how many cells received each label and the score distribution
table(RA1_pred$labels)
summary(RA1_pred$pruned.labels == "NA" | is.na(RA1_pred$pruned.labels))
table(RA2_pred$labels)
summary(RA2_pred$pruned.labels == "NA" | is.na(RA2_pred$pruned.labels))
table(RA4_pred$labels)
summary(RA4_pred$pruned.labels == "NA" | is.na(RA4_pred$pruned.labels))
table(RA5_pred$labels)
summary(RA5_pred$pruned.labels == "NA" | is.na(RA5_pred$pruned.labels))
table(RA6_pred$labels)
summary(RA6_pred$pruned.labels == "NA" | is.na(RA6_pred$pruned.labels))

# Add the SingleR labels to the Seurat metadata
RA1_seurat$SingleR_label <- RA1_pred$labels
RA1_seurat$SingleR_pruned <- RA1_pred$pruned.labels
RA2_seurat$SingleR_label <- RA2_pred$labels
RA2_seurat$SingleR_pruned <- RA2_pred$pruned.labels
RA4_seurat$SingleR_label <- RA4_pred$labels
RA4_seurat$SingleR_pruned <- RA4_pred$pruned.labels
RA5_seurat$SingleR_label <- RA5_pred$labels
RA5_seurat$SingleR_pruned <- RA5_pred$pruned.labels
RA6_seurat$SingleR_label <- RA6_pred$labels
RA6_seurat$SingleR_pruned <- RA6_pred$pruned.labels

# 5. Summarize annotations per cluster
# Each cluster gets its most frequent label
RA1_cluster_annotation <- data.frame(
  cluster = RA1_seurat$seurat_clusters,
  label = RA1_pred$pruned.labels
) %>%
  group_by(cluster, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = n, n = 1)

RA2_cluster_annotation <- data.frame(
  cluster = RA2_seurat$seurat_clusters,
  label = RA2_pred$pruned.labels
) %>%
  group_by(cluster, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = n, n = 1)

RA4_cluster_annotation <- data.frame(
  cluster = RA4_seurat$seurat_clusters,
  label = RA4_pred$pruned.labels
) %>%
  group_by(cluster, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = n, n = 1)

RA5_cluster_annotation <- data.frame(
  cluster = RA5_seurat$seurat_clusters,
  label = RA5_pred$pruned.labels
) %>%
  group_by(cluster, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = n, n = 1)

RA6_cluster_annotation <- data.frame(
  cluster = RA6_seurat$seurat_clusters,
  label = RA6_pred$pruned.labels
) %>%
  group_by(cluster, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = n, n = 1)

print(RA1_cluster_annotation)
print(RA2_cluster_annotation)
print(RA4_cluster_annotation)
print(RA5_cluster_annotation)
print(RA6_cluster_annotation)

# Visualize the annotations on the UMAP
DimPlot(RA1_seurat, reduction = "umap", group.by = "SingleR_pruned", label = TRUE) +
  ggtitle("SingleR annotations (Monaco reference)")
DimPlot(RA2_seurat, reduction = "umap", group.by = "SingleR_pruned", label = TRUE) +
  ggtitle("SingleR annotations (Monaco reference)")
DimPlot(RA4_seurat, reduction = "umap", group.by = "SingleR_pruned", label = TRUE) +
  ggtitle("SingleR annotations (Monaco reference)")
DimPlot(RA5_seurat, reduction = "umap", group.by = "SingleR_pruned", label = TRUE) +
  ggtitle("SingleR annotations (Monaco reference)")
DimPlot(RA6_seurat, reduction = "umap", group.by = "SingleR_pruned", label = TRUE) +
  ggtitle("SingleR annotations (Monaco reference)")
