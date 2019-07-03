# Load packages
library(loomR)
library(Seurat)

# Load seurat object (query)
epid_seu <- readRDS("data/epid_seurat_filtered_normalized_annotated.rds")

# Connect to loom object
brain_loom <- connect(filename = "data/TM_brain_ref.loom")
brain_loom 

# Convert to Seurat
brain_seu <- as.Seurat(brain_loom)

# Pre-processing
brain_seu <- NormalizeData(
  object = brain_seu, 
  normalization.method = "LogNormalize",
  scale.factor = 10000, 
  display.progress = TRUE
)
brain_seu <- FindVariableFeatures(epid_seu, selection.method = "vst", nfeatures = 2000)
brain_seu <- ScaleData(brain_seu)
brain_seu <- RunPCA(brain_seu, features = VariableFeatures(object = epid_seu))

# Project data
anchors <- FindTransferAnchors(reference = brain_seu, query = epid_seu)
predictions <- TransferData(anchorset = anchors, refdata = brain_seu$Class)
epid_seu <- AddMetaData(epid_seu, metadata = predictions)

# Save Seurat object