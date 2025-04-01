library(Seurat)
library(ggplot2)
library(dplyr)
library(hdf5r)
library(ggfortify)
library(harmony)
library(sctransform)

# Définir le chemin des données
dir_path <- "/mnt/datadisk/Jordan/Data/"


## MpBC1
## -----
### On load les données 10X Sptiale
MpBC1 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC1"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC1")
### On load les annotations (avec les spots non annotés)
annotation_MpBC1_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC1/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC1"))    # Add MpBC1 to Barcode column
### On ajoute les annotations aux meta.data
MpBC1@meta.data <- cbind(MpBC1@meta.data, annotation_MpBC1_with_unlabels)

# ## MpBC2
# ## -----
# MpBC2 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC2"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC2")
# annotation_MpBC2_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC2/Annotations_new_with_unlabels.csv")) %>%
#   mutate(Barcode = paste0(Barcode, "-MpBC2"))
# MpBC2@meta.data <- cbind(MpBC2@meta.data, annotation_MpBC2_with_unlabels)

## MpBC3
## -----
MpBC3 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC3"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC3")
annotation_MpBC3_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC3/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC3"))
MpBC3@meta.data <- cbind(MpBC3@meta.data, annotation_MpBC3_with_unlabels)
# 
# ## MpBC4
# ## -----
# MpBC4 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC4"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC4")
# annotation_MpBC4_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC4/Annotations_new_with_unlabels.csv")) %>%
#   mutate(Barcode = paste0(Barcode, "-MpBC4"))
# MpBC4@meta.data <- cbind(MpBC4@meta.data, annotation_MpBC4_with_unlabels)


"MpBC1" -> MpBC1@meta.data$orig.ident
# "MpBC2" -> MpBC2@meta.data$orig.ident
"MpBC3" -> MpBC3@meta.data$orig.ident
# "MpBC4" -> MpBC4@meta.data$orig.ident

# Filtrage des cellules avec nCount_Spatial > 0
MpBC1 <- subset(MpBC1, nCount_Spatial > 0)
# MpBC2 <- subset(MpBC2, nCount_Spatial > 0)
MpBC3 <- subset(MpBC3, nCount_Spatial > 0)
# MpBC4 <- subset(MpBC4, nCount_Spatial > 0)

MpBC1 <- MpBC1 %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")
# MpBC2 <- MpBC2 %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")
MpBC3 <- MpBC3 %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")
# MpBC4 <- MpBC4 %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")


# Fusion des objets Seurat
seurat_object <- merge(MpBC1, y = list(MpBC3), add.cell.ids = c("MpBC1", "MpBC3"), project = "MpBC")

# Normalisation et identification des gènes variables
# seurat_object <- NormalizeData(seurat_object)

seurat_object <- SCTransform(seurat_object,
                             vst.flavor = "v2",
                             assay = "Spatial",
                             do.scale = FALSE,
                             do.center = TRUE,
                             verbose = TRUE,
                             n_genes = 1000)

# Création des ancres pour l'intégration
seurat_list <- SplitObject(seurat_object, 
                           split.by = "orig.ident")

# seurat_object <- FindVariableFeatures(seurat_object, 
#                                       selection.method = "vst", 
#                                       nfeatures = 3000)

features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000)

seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                  assay = "SCT",
                                  anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                  dims = 1:30, 
                                  reduction = "cca", 
                                  normalization.method = "SCT",
                                  anchor.features = features)

seurat_object <- IntegrateData(anchorset = anchors,
                               dims = 1:30,
                               normalization.method = "SCT")

# Scaling et PCA
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, npcs = 30)

# Construction du graphe des voisins et clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Correction du batch
seurat_object <- RunHarmony(seurat_object,
                            group.by.vars = "orig.ident",
                            theta = 2,
                            lambda = 1,
                            sigma = 0.1)

# Exécution de l'UMAP
seurat_object <- RunUMAP(seurat_object, 
                         assay = "integrated",
                         reduction = "harmony",
                         dims = 1:20,
                         n_neighbors = 100,
                         min.dist = 0.001,
                         spread = 1.5,
                         metric = "cosine")

# Visualisation
DimPlot(seurat_object, 
        reduction = "umap", 
        group.by = "Annotations_new", 
        label = TRUE, 
        repel = TRUE, 
        pt.size = 1) +
  theme_minimal()
