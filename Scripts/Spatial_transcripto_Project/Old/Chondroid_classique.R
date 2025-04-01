library(Seurat)
library(ggplot2)
library(dplyr)
library(hdf5r)
library(ggfortify)
library(harmony)
library(sctransform)

# Définir le chemin des données
dir_path <- "/mnt/datadisk/Jordan/Data/"

## MpBC9
## -----
MpBC9 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC9"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC9")
annotation_MpBC9_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC9/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC9"))
MpBC9@meta.data <- cbind(MpBC9@meta.data, annotation_MpBC9_with_unlabels)

## MpBC13
## ------
MpBC13 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC13"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC13")
annotation_MpBC13_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC13/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC13"))
MpBC13@meta.data <- cbind(MpBC13@meta.data, annotation_MpBC13_with_unlabels)

## MpBC15
## ------
MpBC15 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC15"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC15")
annotation_MpBC15_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC15/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC15"))
MpBC15@meta.data <- cbind(MpBC15@meta.data, annotation_MpBC15_with_unlabels)

## MpBC16
## ------
MpBC16 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC16"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC16")
annotation_MpBC16_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC16/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC16"))
MpBC16@meta.data <- cbind(MpBC16@meta.data, annotation_MpBC16_with_unlabels)

"MpBC9" -> MpBC9@meta.data$orig.ident
"MpBC13" -> MpBC13@meta.data$orig.ident
"MpBC15" -> MpBC15@meta.data$orig.ident
"MpBC16" -> MpBC16@meta.data$orig.ident

# Filtrage des cellules avec nCount_Spatial > 0
MpBC9 <- subset(MpBC9, nCount_Spatial > 0)
MpBC13 <- subset(MpBC13, nCount_Spatial > 0)
MpBC15 <- subset(MpBC15, nCount_Spatial > 0)
MpBC16 <- subset(MpBC16, nCount_Spatial > 0)

# Normalisation avec SCTransform pour chaque patient
MpBC9 <- SCTransform(MpBC9,
                     vst.flavor = "v2",
                     assay = "Spatial",
                     do.scale = FALSE,
                     do.center = TRUE,
                     verbose = FALSE)
MpBC13 <- SCTransform(MpBC13,
                      vst.flavor = "v2",
                      assay = "Spatial",
                      do.scale = FALSE,
                      do.center = TRUE,
                      verbose = FALSE)
MpBC15 <- SCTransform(MpBC15,
                      vst.flavor = "v2",
                      assay = "Spatial",
                      do.scale = FALSE,
                      do.center = TRUE,
                      verbose = FALSE)
MpBC16 <- SCTransform(MpBC16,
                      vst.flavor = "v2",
                      assay = "Spatial",
                      do.scale = FALSE,
                      do.center = TRUE,
                      verbose = FALSE)

merged_seurat <- merge(MpBC9, 
                       y = list(MpBC13, MpBC15, MpBC16), 
                       add.cell.ids = c("MpBC9", "MpBC13", "MpBC15", "MpBC16"))

# Convertir les données dans orig.ident en facteur
merged_seurat$orig.ident <- factor(merged_seurat$orig.ident, levels = c("MpBC9", "MpBC13", "MpBC15", "MpBC16"))

# Filtrer les NAs dans Annotations_new
merged_seurat <- merged_seurat %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")

VariableFeatures(merged_seurat) <- SelectIntegrationFeatures(object.list = list(MpBC9, MpBC13, MpBC15, MpBC16), 
                                                             assay = "SCT",
                                                             nfeatures = 100)

# merged_seurat <- FindVariableFeatures(merged_seurat, 
#                                       selection.method = "vst", 
#                                       nfeatures = 3000)

merged_seurat <- RunPCA(merged_seurat, 
                        npcs = 30, 
                        verbose = FALSE)

ElbowPlot(merged_seurat, reduction = "pca")


merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident")

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30, reduction = "harmony")

merged_seurat <- FindClusters(merged_seurat, resolution = 0.5, algorithm = 1)

merged_seurat <- RunUMAP(merged_seurat, 
                         dims = 1:30,
                         reduction = "harmony",
                         n.neighbors = 50, 
                         min.dist = 0.0001,
                         spread = 2, 
                         n.epochs = 500)


DimPlot(merged_seurat,
              reduction = "umap",
              group.by = "Annotations_new",
              pt.size = 2)



# merged_seurat <- RunTSNE(merged_seurat, dims = 1:30, reduction = "harmony")
# 
# DimPlot(merged_seurat,
#               reduction = "tsne",
#               group.by = "Annotations_new",
#               pt.size = 2)
