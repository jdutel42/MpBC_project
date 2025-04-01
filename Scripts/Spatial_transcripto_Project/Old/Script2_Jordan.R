dir_path <- "/mnt/datadisk/Jordan/Data/"


# Charger les données

## MpBC1
## -----
### On load les données 10X Sptiale
MpBC1 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC1"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC1")
### On load les annotations (avec les spots non annotés)
annotation_MpBC1_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC1/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC1"))    # Add MpBC1 to Barcode column
### On ajoute les annotations aux meta.data
MpBC1@meta.data <- cbind(MpBC1@meta.data, annotation_MpBC1_with_unlabels)

## MpBC2
## -----
MpBC2 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC2"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC2")
annotation_MpBC2_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC2/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC2"))
MpBC2@meta.data <- cbind(MpBC2@meta.data, annotation_MpBC2_with_unlabels)

## MpBC3
## -----
MpBC3 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC3"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC3")
annotation_MpBC3_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC3/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC3"))
MpBC3@meta.data <- cbind(MpBC3@meta.data, annotation_MpBC3_with_unlabels)

## MpBC4
## -----
MpBC4 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC4"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC4")
annotation_MpBC4_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC4/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC4"))
MpBC4@meta.data <- cbind(MpBC4@meta.data, annotation_MpBC4_with_unlabels)

## MpBC5
## -----
MpBC5 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC5"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC5")
annotation_MpBC5_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC5/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC5"))
MpBC5@meta.data <- cbind(MpBC5@meta.data, annotation_MpBC5_with_unlabels)

## MpBC6
## -----
MpBC6 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC6"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC6")
annotation_MpBC6_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC6/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC6"))
MpBC6@meta.data <- cbind(MpBC6@meta.data, annotation_MpBC6_with_unlabels)

## MpBC7
## -----
MpBC7 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC7"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC7")
annotation_MpBC7_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC7/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC7"))
MpBC7@meta.data <- cbind(MpBC7@meta.data, annotation_MpBC7_with_unlabels)

## MpBC8
## -----
MpBC8 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch1/Mp_BC8"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC8")
annotation_MpBC8_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch1/Mp_BC8/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC8"))
MpBC8@meta.data <- cbind(MpBC8@meta.data, annotation_MpBC8_with_unlabels)

## MpBC9
## -----
MpBC9 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC9"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC9")
annotation_MpBC9_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC9/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC9"))
MpBC9@meta.data <- cbind(MpBC9@meta.data, annotation_MpBC9_with_unlabels)

## MpBC10
## ------
MpBC10 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC10"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC10")
annotation_MpBC10_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC10/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC10"))
MpBC10@meta.data <- cbind(MpBC10@meta.data, annotation_MpBC10_with_unlabels)

## MpBC11
## ------
MpBC11 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC11"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC11")
annotation_MpBC11_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC11/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC11"))
MpBC11@meta.data <- cbind(MpBC11@meta.data, annotation_MpBC11_with_unlabels)

## MpBC13
## ------
MpBC13 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC13"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC13")
annotation_MpBC13_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC13/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC13"))
MpBC13@meta.data <- cbind(MpBC13@meta.data, annotation_MpBC13_with_unlabels)

## MpBC14
## ------
MpBC14 <- Load10X_Spatial(data.dir = paste0(dir_path, "Visium_MpBC_batch2/MpBC14"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "Slice_MpBC14")
annotation_MpBC14_with_unlabels <- read.csv(paste0(dir_path, "Visium_MpBC_batch2/MpBC14/Annotations_new_with_unlabels.csv")) %>%
  mutate(Barcode = paste0(Barcode, "-MpBC14"))
MpBC14@meta.data <- cbind(MpBC14@meta.data, annotation_MpBC14_with_unlabels)

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

# MpBC10 et MpBC16 n'ont pas au total 4992 barcodes/cellules/spots (annotées ou non) alors que les autres samples si ??? ==> Étrange, problème de séquençage ??? ==> En fait c'est nomal car certains certaines tumeurs petites et certains spots n'ont pas été séquencé car il n'y avait rien à séquencer

"MpBC1" -> MpBC1@meta.data$orig.ident
"MpBC2" -> MpBC2@meta.data$orig.ident
"MpBC3" -> MpBC3@meta.data$orig.ident
"MpBC4" -> MpBC4@meta.data$orig.ident
"MpBC5" -> MpBC5@meta.data$orig.ident
"MpBC6" -> MpBC6@meta.data$orig.ident
"MpBC7" -> MpBC7@meta.data$orig.ident
"MpBC8" -> MpBC8@meta.data$orig.ident
"MpBC9" -> MpBC9@meta.data$orig.ident
"MpBC10" -> MpBC10@meta.data$orig.ident
"MpBC11" -> MpBC11@meta.data$orig.ident
"MpBC13" -> MpBC13@meta.data$orig.ident
"MpBC14" -> MpBC14@meta.data$orig.ident
"MpBC15" -> MpBC15@meta.data$orig.ident
"MpBC16" -> MpBC16@meta.data$orig.ident

# Filtrage des cellules avec nCount_Spatial > 0
MpBC1 <- subset(MpBC1, nCount_Spatial > 0)
MpBC2 <- subset(MpBC2, nCount_Spatial > 0)
MpBC3 <- subset(MpBC3, nCount_Spatial > 0)
MpBC4 <- subset(MpBC4, nCount_Spatial > 0)
MpBC5 <- subset(MpBC5, nCount_Spatial > 0)
MpBC6 <- subset(MpBC6, nCount_Spatial > 0)
MpBC7 <- subset(MpBC7, nCount_Spatial > 0)
MpBC8 <- subset(MpBC8, nCount_Spatial > 0)
MpBC9 <- subset(MpBC9, nCount_Spatial > 0)
MpBC10 <- subset(MpBC10, nCount_Spatial > 0)
MpBC11 <- subset(MpBC11, nCount_Spatial > 0)
MpBC13 <- subset(MpBC13, nCount_Spatial > 0)
MpBC14 <- subset(MpBC14, nCount_Spatial > 0)
MpBC15 <- subset(MpBC15, nCount_Spatial > 0)
MpBC16 <- subset(MpBC16, nCount_Spatial > 0)

# ------------------------------------------------------------------------------

# On a nos objets Seurat individuels, il faut les combiner un en seul objet


# Normalisation avec SCTransform pour chaque patient
sample_list <- list(MpBC9, MpBC13, MpBC15, MpBC16)  # Liste des Seurat objets

sample_list <- lapply(sample_list, function(x) {
  SCTransform(x, 
              vst.flavor = "v2", 
              assay = "Spatial", 
              verbose = FALSE,
              new.assay.name = "SCT")
})

# Trouver les ancres pour aligner les échantillons
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 100)

# Il faut que les échantillons aient des features communs
common_features <- Reduce(intersect, lapply(sample_list, function(x) rownames(x@assays$SCT@scale.data)))
features <- intersect(features, common_features)


sample_list <- PrepSCTIntegration(object.list = sample_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = sample_list, normalization.method = "SCT", anchor.features = features)

# Intégration des données
merged_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")



# Filtrer les NAs dans Annotations_new
merged_seurat <- merged_seurat %>% subset(Annotations_new != "" & Annotations_new != "Cluster 1" & Annotations_new != "Cluster 3")

Idents(merged_seurat) <- "Annotations_new"

degs <- FindAllMarkers(merged_seurat, assay = "integrated", only.pos = TRUE)  


top_genes <- degs %>% group_by(cluster) %>% top_n(500, avg_log2FC) %>% pull(gene)


VariableFeatures(merged_seurat) <- top_genes


# Réduction de dimension PCA
merged_seurat <- RunPCA(merged_seurat,
                        assay = "integrated",
                        npcs = 30,
                        features = top_genes)

# Correction du batch
merged_seurat <- RunHarmony(merged_seurat,
                            group.by.vars = "orig.ident",
                            theta = 2,
                            lambda = 1,
                            sigma = 0.1)

merged_seurat <- RunUMAP(merged_seurat,
                         assay = "integrated",
                         reduction = "harmony",
                         dims = 1:10,
                         n.neighbors = 50, 
                         min.dist = 0.0001,
                         spread = 2, 
                         n.epochs = 500,
                         metric = "manhattan")
DimPlot(merged_seurat,
        pt.size = 2,
        reduction = "umap",
        group.by = "Annotations_new",
        shape.by = "orig.ident")
