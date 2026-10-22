# Install Seurat and hdf5r
install.packages("Seurat")
install.packages("hdf5r")
library("Seurat")
library("hdf5r")


## Step 1.i. Data from each individual will be imported as a Seurat object

# Read in the hd5 file
counts <- Read10X_h5("~/galaxy_inputs/1_RA1_raw_feature_bc_matrix.h5.h5", use.names = TRUE, unique.features = TRUE)

# Need to limit to the Gene Expression features (there are also antibody counts)
rna_counts <- counts[["Gene Expression"]]

# Create the Seurat object
RA1_seurat <- CreateSeuratObject(
  counts = rna_counts,
  project = "RA1",
  min.cells = 3,
  min.features = 200
)

# Create the correct timepoint definitions
RA1_key <- read.delim("~/galaxy_input/7_RA1_QC_categories.txt.tabular") 
RA1_key$timepoint <- sub(".*_(Baseline|3w_MTX|6w_MTX|12w_MTX)$", "\\1", RA1_key$sample)
RA1_key$timepoint <- c(Baseline = "BL", "3w_MTX" = "wk3",
                   "6w_MTX" = "wk6", "12w_MTX" = "wk12")[RA1_key$timepoint]

# Keep only singlets with a timepoint that pass QC
# (if pass_QC was read as text, use key$pass_QC == "True")
RA1_key <- RA1_key[!is.na(RA1_key$timepoint) & RA1_key$pass_QC == "True", ]
RA1_key$timepoint <- factor(RA1_key$timepoint, levels = c("BL", "wk3", "wk6", "wk12"))

table(RA1_key$timepoint, useNA = "ifany")
table(RA1_key$hto_demux, RA1_key$timepoint)   # hto_demux should agree with timepoint


# Seurat cell names from the merge are "<orig.ident>_<barcode>"
RA1_seurat$cell_id <- paste(RA1_seurat$orig.ident, sub("^[A-Z0-9]+_", "", colnames(RA1_seurat)), sep = "_")

# Create the timepoint and cell ID columns
RA1_seurat$timepoint <- RA1_key$timepoint[match(RA1_seurat$cell_id, RA1_key$cell_id)]
RA1_seurat$hto_call  <- RA1_key$hto_demux[match(RA1_seurat$cell_id, RA1_key$cell_id)]

# Check: most cells should have a timepoint; NAs are cells removed by HTO QC
table(RA1_seurat$orig.ident, RA1_seurat$timepoint, useNA = "ifany")

# Repeat above steps for RA2
counts <- Read10X_h5("~/galaxy_inputs/4_RA2_raw_feature_bc_matrix.h5.h5", use.names = TRUE, unique.features = TRUE)
rna_counts <- counts[["Gene Expression"]]
RA2_seurat <- CreateSeuratObject(
  counts = rna_counts,
  project = "RA2",
  min.cells = 3,
  min.features = 200
)


# Repeat above steps for RA4
counts <- Read10X_h5("~/galaxy_inputs/3_RA4_raw_feature_bc_matrix.h5.h5", use.names = TRUE, unique.features = TRUE)
rna_counts <- counts[["Gene Expression"]]
RA4_seurat <- CreateSeuratObject(
  counts = rna_counts,
  project = "RA4",
  min.cells = 3,
  min.features = 200
)

# Repeat above steps for RA5
counts <- Read10X_h5("~/galaxy_inputs/2_RA5_raw_feature_bc_matrix.h5.h5", use.names = TRUE, unique.features = TRUE)
rna_counts <- counts[["Gene Expression"]]
RA5_seurat <- CreateSeuratObject(
  counts = rna_counts,
  project = "RA5",
  min.cells = 3,
  min.features = 200
)

# Repeat above steps for RA6
counts <- Read10X_h5("~/galaxy_inputs/1_RA6_raw_feature_bc_matrix.h5.h5", use.names = TRUE, unique.features = TRUE)
rna_counts <- counts[["Gene Expression"]]
RA6_seurat <- CreateSeuratObject(
  counts = rna_counts,
  project = "RA6",
  min.cells = 3,
  min.features = 200
)


## Step 1.ii.	Filtering will be performed to remove cells with <300 features or >10% mitochondrial reads
# Compute percent mitochondrial reads per cell
RA1_seurat[["percent.mt"]] <- PercentageFeatureSet(RA1_seurat, pattern = "^MT-")

# Remove cells with <300 features or >10% mitochondrial reads
RA1_seurat <- subset(
  RA1_seurat,
  subset = nFeature_RNA >= 300 & percent.mt <= 10
)

# Check how many cells remain
RA1_seurat

# Repeat for RA2
RA2_seurat[["percent.mt"]] <- PercentageFeatureSet(RA2_seurat, pattern = "^MT-")
RA2_seurat <- subset(
  RA2_seurat,
  subset = nFeature_RNA >= 300 & percent.mt <= 10
)
RA2_seurat

# Repeat for RA4
RA4_seurat[["percent.mt"]] <- PercentageFeatureSet(RA4_seurat, pattern = "^MT-")
RA4_seurat <- subset(
  RA4_seurat,
  subset = nFeature_RNA >= 300 & percent.mt <= 10
)
RA4_seurat

# Repeat for RA5
RA5_seurat[["percent.mt"]] <- PercentageFeatureSet(RA5_seurat, pattern = "^MT-")
RA5_seurat <- subset(
  RA5_seurat,
  subset = nFeature_RNA >= 300 & percent.mt <= 10
)
RA5_seurat

# Repeat for RA6
RA6_seurat[["percent.mt"]] <- PercentageFeatureSet(RA6_seurat, pattern = "^MT-")
RA6_seurat <- subset(
  RA6_seurat,
  subset = nFeature_RNA >= 300 & percent.mt <= 10
)
RA6_seurat

## Step 1.iii.	Normalization with a scale factor of 10,000
RA1_seurat <- NormalizeData(
  RA1_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
RA2_seurat <- NormalizeData(
  RA2_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
RA4_seurat <- NormalizeData(
  RA4_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
RA5_seurat <- NormalizeData(
  RA5_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
RA6_seurat <- NormalizeData(
  RA6_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
