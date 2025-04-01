# install.packages("Seurat")
# install.packages("SeuratData")
# install.packages("spatstat")
# install.packages("ggplot2")

library(Seurat)
library(ggplot2)

# Liste des échantillons
samples <- c("Visium_MpBC_batch1/Mp_BC1", "Visium_MpBC_batch1/Mp_BC2", "Visium_MpBC_batch1/Mp_BC3", "Visium_MpBC_batch1/Mp_BC4", "Visium_MpBC_batch1/Mp_BC5", "Visium_MpBC_batch1/Mp_BC6", "Visium_MpBC_batch1/Mp_BC7", "Visium_MpBC_batch1/Mp_BC8", "Visium_MpBC_batch2/MpBC9", "Visium_MpBC_batch2/MpBC10", "Visium_MpBC_batch2/MpBC13", "Visium_MpBC_batch2/MpBC14", "Visium_MpBC_batch2/MpBC15", "Visium_MpBC_batch2/MpBC16")

# Charger chaque échantillon séparément
seurat_list <- lapply(samples, function(sample) {
  dir_path <- paste0("/mnt/datadisk/Jordan/Data/", sample)
  seurat_obj <- Load10X_Spatial(data.dir = dir_path, filename = "filtered_feature_bc_matrix.h5")
  #seurat_obj$sample <- sample  # Ajouter une annotation du nom d’échantillon
  return(seurat_obj)
})

# Vérifier les objets chargés
seurat_list

# # Charger les données de transcriptomique spatiale
# visium_data <- Load10X_Spatial(data.dir = "/mnt/datadisk/Jordan/Data/Visium_MpBC_batch1/Mp_BC1", filename = "filtered_feature_bc_matrix.h5")
# # Afficher un résumé des données
# visium_data

# Normalisation des données
visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE)

# Réduction de dimension (PCA & UMAP)
visium_data <- RunPCA(visium_data)
visium_data <- RunUMAP(visium_data, dims = 1:30)

# Clustering des spots
visium_data <- FindNeighbors(visium_data, dims = 1:30)
visium_data <- FindClusters(visium_data, resolution = 0.5)

# Visualisation des clusters
p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)
p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE)
p1 + p2

# Les annotations ne sont pas les bonnes