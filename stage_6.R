
# Stage 6.	Subsetting into cell lineages
## 6.i.	Subset into T cells, B cells, and myeloid cells
t_cells   <- subset(RA1_seurat, subset = SingleR_pruned %in% c("CD4+ T cells", "CD8+ T cells", "T cells"))
b_cells   <- subset(RA1_seurat, subset = SingleR_pruned %in% c("B cells"))
myeloid   <- subset(RA1_seurat, subset = SingleR_pruned %in% c("Monocytes", "Non classical monocytes", "Dendritic cells"))

# 2.ii.	Integration of each subset across batches to account for variance due to batch effect using anchor-based CCA intregration
integrate_subset <- function(obj, batch_col = "batch", dims = 1:30) {
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[batch_col]][, 1])

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)

  anchors <- FindIntegrationAnchors(
    object.list = NULL,
    reduction = "cca",
    dims = dims
  )
  integrated <- IntegrateData(anchorset = anchors, dims = dims)

  DefaultAssay(integrated) <- "integrated"
  integrated
}

RA1_t_int <- integrate_subset(t_cells, batch_col = "orig.ident")
RA1_b_int <- integrate_subset(b_cells, batch_col = "orig.ident")
RA1_myeloid_int <- integrate_subset(myeloid, batch_col = "orig.ident")

# Save the integrated objects
saveRDS(t_int, "t_cells_integrated.rds")
saveRDS(b_int, "b_cells_integrated.rds")
saveRDS(myeloid_int, "myeloid_integrated.rds")
