library(Seurat)

MpBC9 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC9"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC1")
### On load les annotations (avec les spots non annotés)
annotation_MpBC9_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC9/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC9"))    # Add MpBC1 to Barcode column
### On ajoute les annotations aux meta.data
MpBC9@meta.data <- cbind(MpBC9@meta.data, annotation_MpBC9_with_unlabels)


"MpBC16" -> MpBC16@meta.data$orig.ident
MpBC16 <- subset(MpBC16, nCount_Spatial > 0)

MpBC16 <- MpBC16 %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")


# Normalisation des données
seurat_object <- NormalizeData(MpBC16)

# Sélection des gènes variables
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Scaling des données
seurat_object <- ScaleData(seurat_object)

# Calcul des composantes principales
seurat_object <- RunPCA(seurat_object, npcs = 30)

# Construction du graphe des voisins
seurat_object <- FindNeighbors(seurat_object, dims = 1:20)

# Clustering (peut être ajusté selon le niveau de granularité souhaité)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Exécution de l'UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:20)

DimPlot(seurat_object, reduction = "umap", group.by = "Annotations_new", label = TRUE, repel = TRUE) +
  theme_minimal()
