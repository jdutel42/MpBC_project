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

MpBC9[["RNA"]] <- MpBC9[["Spatial"]]
MpBC13[["RNA"]] <- MpBC13[["Spatial"]]
MpBC15[["RNA"]] <- MpBC15[["Spatial"]]
MpBC16[["RNA"]] <- MpBC16[["Spatial"]]

# Sélection des features communes
features <- SelectIntegrationFeatures(object.list = list(MpBC9, MpBC13, MpBC15, MpBC16), nfeatures = 3000)

# Sélectionne les gènes d'ancrage communs à tous les objets
common_genes <- Reduce(intersect, lapply(list(MpBC9, MpBC13, MpBC15, MpBC16), rownames))
features <- common_genes

# Préparation des objets pour l'intégration
prepsct <- PrepSCTIntegration(object.list = list(MpBC9, MpBC13, MpBC15, MpBC16), anchor.features = features)

# Trouver les ancres d'intégration
anchors <- FindIntegrationAnchors(object.list = prepsct, normalization.method = "SCT", anchor.features = features, dims = 1:30)

# Intégration des datasets
integrated_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

# Suite du pipeline : PCA, UMAP, clustering
integrated_seurat <- RunPCA(integrated_seurat)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30)
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)

# Visualisation
DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident")

