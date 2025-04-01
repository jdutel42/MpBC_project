```{r}
# samples <- c("MpBC1", "MpBC2", "MpBC3", "MpBC4", "MpBC5", "MpBC6",
#              "MpBC7", "MpBC8", "MpBC9", "MpBC10", "MpBC11", "MpBC13",
#              "MpBC14", "MpBC15", "MpBC16")
# 
# base_path <- "/mnt/datadisk/Jordan/Data/Visium/"
# 
# for (sample in samples) {
#   file_path <- paste0(base_path, sample, "/filtered_feature_bc_matrix.h5")
#   count_matrix <- Read10X_h5(file_path)
# 
#   obj <- get(sample)
# 
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[1]] <- count_matrix@Dimnames[[1]]
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- count_matrix@Dimnames[[2]]
#   
#   # Rajouter MpBCX_ pour chaque barcodes
#   new_barcode_vec <- c()
#   ## On stocke tous les barcodes de la matrice de comptage qui ne sont pas présents dans les barcodes filtrés dans metadata (via les annotations) 
#   # filtered_out_barcode <- setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]], obj@meta.data$Barcode)
#   for (old_barcode in obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]]) {
#     new_barcode <- paste0(sample, "_", old_barcode)
#     new_barcode_vec <- c(new_barcode_vec, new_barcode)
#   }
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- new_barcode_vec
#   
#   # ON filtre les spots avec une annotations non désirée
#   obj <- subset(obj,
#                 subset = (nCount_Spatial > 0 &
#                           Annotations_old != "" & Annotations_old != "Adipose" &
#                           Annotations_old != "Artifacts" & Annotations_old != "Immune cells" &
#                           Annotations_old != "Intermediate tumour cells" & Annotations_old != "Mixed cells" &
#                           Annotations_old != "Necrosis" & Annotations_old != "Necrotic and apoptotic tissue" &
#                           Annotations_old != "NST surrounded by spindle" & Annotations_old != "Scar-like fibrous stroma" &
#                           Annotations_new != "Pleiomorphic tumor"))
#   
# 
#   new_filtered_out_barcode <- setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]], obj@meta.data$Barcode)
#   
#   # Réassigner l'objet mis à jour
#   assign(sample, obj)
#   
#   if (length(new_filtered_out_barcode) == 0) {
#     print(paste("✅ Données mises à jour pour", sample))
#   } else {
#     print(paste("❌ Erreur pour", sample))
#   }
# 
# }
# 
# samples <- list(MpBC1, MpBC2, MpBC3, MpBC4, MpBC5, MpBC6, MpBC7, MpBC8, MpBC9, MpBC10, MpBC11, MpBC13, MpBC14, MpBC15, MpBC16)
# names(samples) <- c("MpBC1", "MpBC2", "MpBC3", "MpBC4", "MpBC5", "MpBC6", "MpBC7", "MpBC8", "MpBC9", "MpBC10", "MpBC11", "MpBC13", "MpBC14", "MpBC15", "MpBC16")
# samples_names <- names(samples)
# 

```

```{r}
# # samples <- c("MpBC1", "MpBC2", "MpBC3", "MpBC4", "MpBC5", "MpBC6", 
# #              "MpBC7", "MpBC8", "MpBC9", "MpBC10", "MpBC11", "MpBC13", 
# #              "MpBC14", "MpBC15", "MpBC16")
# samples <- c("MpBC1")
# 
# base_path <- "/mnt/datadisk/Jordan/Data/Visium/"
# 
# for (sample in samples) {
#   file_path <- paste0(base_path, sample, "/filtered_feature_bc_matrix.h5")
#   count_matrix <- Read10X_h5(file_path)
# 
#   obj <- get(sample)
# 
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[1]] <- count_matrix@Dimnames[[1]]
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- count_matrix@Dimnames[[2]]
#   
#   # Rajouter MpBCX_ pour chaque barcodes
#   new_barcode_vec <- c()
#   
#   ## On stocke tous les barcodes de la matrice de comptage qui ne sont pas présents dans les barcodes filtrés dans metadata (via les annotations) 
#   filtered_out_barcode <- setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]], obj@meta.data$Barcode)
# 
#   for (old_barcode in obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]]) {
#     new_barcode <- paste0(sample, "_", old_barcode)
#     new_barcode_vec <- c(new_barcode_vec, new_barcode)
#   }
# 
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- new_barcode_vec
# 
#   # for (new_barcode in obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]]) {
#   #   if (new_barcode %in% setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]])) {
#   #     # ON retire cette colonne de la matrice de comptage
#   #     obj@assays[["Spatial"]]@layers[["counts"]] <- obj@assays[["Spatial"]]@layers[["counts"]][, -which(colnames(obj@assays[["Spatial"]]@layers[["counts"]]) == new_barcode)]
#   #   } else {
#   #     next
#   #   }
#   # 
#   # }
#   # 
#   # # Renommer tous les barcodes en ajoutant le préfixe "sample_"
#   # new_barcode_vec <- paste0(sample, "_", obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]])
#   # obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- new_barcode_vec
#   # 
#   # # Filtrer la matrice de comptage en supprimant les colonnes indésirables en une seule étape
#   # keep_columns <- !(new_barcode_vec %in% filtered_out_barcode)
#   # obj@assays[["Spatial"]]@layers[["counts"]] <- obj@assays[["Spatial"]]@layers[["counts"]][, keep_columns]
#   # 
#   # print("OK")
#   # 
#   # # Vérification
#   # if (all(!colnames(obj@assays[["Spatial"]]@layers[["counts"]]) %in% filtered_out_barcode_prefixed)) {
#   #   print(paste("✅ Colonnes correctement filtrées pour", sample))
#   # } else {
#   #   print(paste("❌ Erreur dans le filtrage des colonnes pour", sample))
#   # }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   
#   new_filtered_out_barcode <- setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]], obj@meta.data$Barcode)
#   
#   # Réassigner l'objet mis à jour
#   assign(sample, obj)
#   
#   if (length(new_filtered_out_barcode) == 0) {
#     print(paste("✅ Données mises à jour pour", sample))
#   } else {
#     print(paste("❌ Erreur pour", sample))
#   }
# 
# }
# 
# samples <- list(MpBC1, MpBC2, MpBC3, MpBC4, MpBC5, MpBC6, MpBC7, MpBC8, MpBC9, MpBC10, MpBC11, MpBC13, MpBC14, MpBC15, MpBC16)
# names(samples) <- c("MpBC1", "MpBC2", "MpBC3", "MpBC4", "MpBC5", "MpBC6", "MpBC7", "MpBC8", "MpBC9", "MpBC10", "MpBC11", "MpBC13", "MpBC14", "MpBC15", "MpBC16")
# samples_names <- names(samples)
# 

```

```{r}
# for (sample in samples) {
#   file_path <- paste0(base_path, sample, "/filtered_feature_bc_matrix.h5")
#   count_matrix <- Read10X_h5(file_path)
# 
#   obj <- get(sample)
# 
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[1]] <- count_matrix@Dimnames[[1]]
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- count_matrix@Dimnames[[2]]
# 
#   # Identifier les barcodes à filtrer (ceux qui ne sont pas dans le metadata)
#   filtered_out_barcode <- setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]], obj@meta.data$Barcode)
# 
#   # Appliquer le préfixe "sample_" aux barcodes de la matrice
#   new_barcode_vec <- paste0(sample, "_", obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]])
#   obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]] <- new_barcode_vec
# 
#   # Appliquer le même préfixe aux barcodes filtrés pour qu'ils soient comparables
#   filtered_out_barcode_prefixed <- paste0(sample, "_", filtered_out_barcode)
# 
#   # Filtrer la matrice pour retirer les colonnes indésirables
#   keep_columns <- !(new_barcode_vec %in% filtered_out_barcode_prefixed)
# 
#   # Vérifier la cohérence des dimensions
#   if (length(keep_columns) == ncol(obj@assays[["Spatial"]]@layers[["counts"]])) {
#     obj@assays[["Spatial"]]@layers[["counts"]] <- obj@assays[["Spatial"]]@layers[["counts"]][, keep_columns]
#   } else {
#     stop(paste("❌ Problème de dimension pour", sample))
#   }
# 
#   # Vérification
#   if (all(!colnames(obj@assays[["Spatial"]]@layers[["counts"]]) %in% filtered_out_barcode_prefixed)) {
#     print(paste("✅ Colonnes correctement filtrées pour", sample))
#   } else {
#     print(paste("❌ Erreur dans le filtrage des colonnes pour", sample))
#   }
# 
#   # Vérifier que tous les barcodes restants sont bien dans le metadata
#   new_filtered_out_barcode <- setdiff(obj@assays[["Spatial"]]@layers[["counts"]]@Dimnames[[2]], obj@meta.data$Barcode)
# 
#   # Réassigner l'objet mis à jour
#   assign(sample, obj)
# 
#   if (length(new_filtered_out_barcode) == 0) {
#     print(paste("✅ Données mises à jour pour", sample))
#   } else {
#     print(paste("❌ Erreur pour", sample))
#   }
# }

```


```{r}
# ## INES
# p <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16)
# 
# for ( i in p ){
#   
#   setwd(sprintf("/mnt/datadisk/Jordan/Data/Visium/MpBC%s", i))
# 
#   # Visium files --------------------------------------------------------------
#   
#   assign(sprintf("matrix_%s", i),
#          Load10X_Spatial("./", "/filtered_feature_bc_matrix.h5"))  
# 
#   # Import seurat object and convert it to a matrix object for inferCNV
#   assign(sprintf("matrix_%s", i),
#          importSrat(eval(str2expression(sprintf("matrix_%s",i))),
#                     slot = "counts", 
#                     assay = "Spatial", 
#                     log2tpm_tr = TRUE))
#   
#   
#   # from counts to transcripts per million -------------------------------------
#   
#   assign(sprintf("matrix_%s",i), 
#          umi_to_log2tpm(eval(str2expression(sprintf("matrix_%s",i)))))
#   
#   
#   # Annotation files ----------------------------------------------------------
# 
#   ref <- read.csv(sprintf("Annotations_new_with_unlabels.csv", i),
#                   header = TRUE, 
#                   sep = ",", 
#                   stringsAsFactors = TRUE)
# 
#   # FILTRER LES ANNOTATION !!!
#   
#   
#   # change column name to "Annotations" ----------------------------------------
#   
#   colnames(ref)[2] <- "Annotations"
#   
#   # Retrieve cell types names --------------------------------------------------
#   
#   assign(sprintf("cells%s", i), unique( ref$Annotations))
#   
#   assign(sprintf("ref%s", i), ref)
# 
# }


```



# Create Refgroup
```{r}
# On va créer un vecteur contenant tous les ID des cellules/spots de la matrice de compte (dans l'obet Seurat) considérés comme réference = "Normal épithélium", "Normal mesenchyme", "Blood"...

# On merge les MpBC en un seul objet seurat plus pratique pour ne garder que les cellules reférences
MpBC_merged <- merge(MpBC1, 
                     y = c(MpBC2, MpBC3, MpBC4, MpBC5, MpBC6, MpBC7, MpBC8, 
                           MpBC9, MpBC10, MpBC11, MpBC13, MpBC14, MpBC15, MpBC16), 
                     project = "MpBC_Visium") %>% 
  JoinLayers()


# Convertir les données dans orig.ident en facteur
MpBC_merged$Patient <- factor(MpBC_merged$Patient, 
                              levels = c("MpBC1", "MpBC2", "MpBC3", "MpBC4", "MpBC5", "MpBC6", "MpBC7", "MpBC8", 
                                         "MpBC9", "MpBC10", "MpBC11", "MpBC13", "MpBC14", "MpBC15", "MpBC16"))
MpBC_merged$Seq_batch <- factor(MpBC_merged$Seq_batch, 
                                levels = c("batch1", "batch2"))
MpBC_merged$Visium_slide <- factor(MpBC_merged$Visium_slide, 
                                   levels = c("slide1", "slide2", "slide3", "slide4", 
                                              "slide5", "slide6", "slide7", "slide8"))

print("✅ MpBC_merged créé")

# # ON récupère le barcodes des spots annotés comme "Blood" et "Normal epithelium" et "Normal mesenchyme" pour tous les patients
# Ref_cells_barcodes <- WhichCells(MpBC_merged, 
#                                  expression = (Annotations_old %in% c("Blood")) | 
#                                               (Annotations_new %in% c("Normal epithelium", "Normal mesenchyme")))
# 
# # On extrait la sous matrice de comptage correspondant aux cellules de référence
# Ref_cells_matrix <- MpBC_merged[["Spatial"]]@layers$counts[, Ref_cells_barcodes]

# On crée un objet seurt ne contenant que les cellules de référence de tous les patents (ainsi on a directement notre matrice de reférence)
Ref_MpBC <- subset(MpBC_merged,
                   subset = 
                     (Annotations_old %in% c("Blood")) |
                     (Annotations_new %in% c("Normal epithelium", "Normal mesenchyme")) )

# Sélectionner les cellules normales et récupérer leurs rownames, on stocke les ID des spots dans un vecteur pour plus tard
ref_cells <- Cells(Ref_MpBC)

# On créé des nouveaux nouveaux objets Seurat pour chaque échantillon AVEC les comptages des cellules de référence de tous les patients
for (sample_name in samples_names) {
  
  sample <- get(sample_name)
  
  # On enlève les ref_cells déjà present dans chaque échantillon (sinon il y a des duplications de cellule)
  sample <- subset(sample, 
                   subset = Annotations_old != "Blood" )
  sample <- subset(sample,
                   subset = Annotations_new != "Normal epithelium" )
  sample <- subset(sample,
                   subset = Annotations_new != "Normal mesenchyme")
  
  # On ajoute TOUTES les cellules de reference pour chaque échantillons avec un merge des matrices de comptage
  assign(sprintf("%s_with_ref", sample_name),
         merge(sample, Ref_MpBC) %>% 
           JoinLayers())
  
  print(paste("✅ Création de l'objet Seurat avec les cellules de référence pour", sample_name))
  
}



# # Créer le vecteur ref_obs pour infercnvplus
# for (sample_name in samples_names) {
#   
#   sample <- get(sprintf("%s_with_ref", sample_name))
#   
#   # assign(sprintf("ref_obs_%s", "MpBC1"), 
#   #        WhichCells(MpBC1_with_ref, 
#   #          cells = rownames(MpBC1_with_ref@meta.data[
#   #            !(MpBC1_with_ref@meta.data$Annotations_old %in% c("Blood")) & 
#   #            !(MpBC1_with_ref@meta.data$Annotations_new %in% c("Normal epithelium", "Normal mesenchyme")), 
#   #          ])) )
#   
#   assign(sprintf("ref_obs_%s", sample_name), 
#          rownames(sample@meta.data) %>% setNames(ifelse(. %in% ref_cells, "reference", "observed")) )
#   
#   #        ref_obs <- all_data@meta.data %>%
#   # rownames() %>%
#   # setNames(ifelse(. %in% ref_cells, "reference", "observed"))
# 
#   
#   print(paste("✅ Création du vecteur ref_obs pour", sample_name))
# }
#        
#        
#        
# ref_obs <- all_data@meta.data %>%
# rownames() %>%
# setNames(ifelse(. %in% ref_cells, "reference", "observed"))

```



```{r}

# Charger les librairies nécessaires
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gtools)
library(circlize)  # Pour l'annotation des couleurs

# Charger le dataframe (remplace ceci par ton import)
# df <- read.csv("chemin/vers/tes_donnees.csv")

# Filtrer uniquement les types cellulaires d'intérêt
df_filtered <- df_list[[1]] %>%
  filter(Annotations_new %in% c("Epithelial tumor", "Mesenchymal tumor"))

# Trier les cytobandes par chromosome (converti chromosome en facteur ordonné)
df_filtered <- df_filtered %>%
  mutate(chromosomes = factor(chromosomes, levels = mixedsort(unique(chromosomes))))

# Transformer en matrice avec les cytobandes en colonnes et les cellules en lignes
df_wide <- df_filtered %>%
  select(Barcode, chromosomes, cnv_cyt, CNA_score, Annotations_new) %>%
  pivot_wider(names_from = cnv_cyt, values_from = CNA_score) 

# Convertir en matrice numérique (ignorer les colonnes non numériques)
mat_CNA <- as.matrix(df_wide[, -c(1, 3)])  # Suppression de Barcode et Annotations_new

# Remplacement des NA par 0 (ou autre stratégie, selon besoin)
mat_CNA[is.na(mat_CNA)] <- 0 

mat_CNA <- apply(mat_CNA, 2, as.numeric)

# Calcul du Z-score par colonne (cytoband)
mat_CNA_zscore <- t(scale(t(mat_CNA)))

mat_CNA_zscore[is.na(mat_CNA_zscore)] <- 0 

mat_CNA_zscore <- apply(mat_CNA_zscore, 2, as.numeric)

# Création des annotations pour la heatmap
row_ha <- rowAnnotation(
  Type = df_wide$Annotations_new,
  col = list(Type = c("Epithelial tumor" = "red", "Mesenchymal tumor" = "blue"))
)

# Dessiner la heatmap
Heatmap(
  mat_CNA_zscore,
  name = "Z-score CNA",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,  # Optionnel : Clusteriser les cellules
  cluster_columns = FALSE,  # Ne pas clusteriser les cytobands (respecter l’ordre)
  right_annotation = row_ha
)

```