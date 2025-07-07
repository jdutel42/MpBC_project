import scanpy as sc
import pandas as pd
import os

# Chemin du dossier contenant tous les sous-dossiers d’échantillons
root_dir = "/mnt/datadisk/Jordan/Data/Visium"

# Liste de tous les noms de dossiers (chaque dossier = un échantillon)
# sample_dirs = [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]
sample_dirs = [
    "MpBC1",
    "MpBC2",
    "MpBC3",
    "MpBC4",
    "MpBC5",
    "MpBC6",
    "MpBC7",
    "MpBC8",
    "MpBC9",
    "MpBC10",
    "MpBC11",
    "MpBC13",
    "MpBC14",
    "MpBC15",
    "MpBC16",
]

# Dictionnaire pour stocker tous les objets AnnData par type cellulaire
cell_type_dict = {}

for sample in sample_dirs:
    print(f"→ Traitement de : {sample}")
    
    sample_path = os.path.join(root_dir, sample)
    h5_file = os.path.join(sample_path, "filtered_feature_bc_matrix.h5")
    annot_file = os.path.join(sample_path, "Annotations_old_with_unannotated.csv")
    
    if not os.path.exists(h5_file) or not os.path.exists(annot_file):
        print(f"  Fichiers manquants pour {sample}, on saute.")
        continue
    
    # Charger la matrice Visium
    adata = sc.read_10x_h5(h5_file)

    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)

    # Préfixer les barcodes pour éviter conflits
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    
    # Charger les annotations
    annot = pd.read_csv(annot_file, dtype={"Barcode": str})
    annot['Barcode'] = annot['Barcode'].apply(lambda bc: f"{sample}_{bc}")
    
    # Merge des annotations
    adata.obs = adata.obs.merge(annot, left_index=True, right_on='Barcode', how='left')
    adata.obs.set_index('Barcode', inplace=True)  # ← très important pour que l'index reste correct
    
    # Filtrage uniquement des spots annotés
    # adata = adata[adata.obs['Annotations_old_with_unannotated'].notna()].copy()
    
    # Ajouter une colonne "sample" dans .obs
    adata.obs["sample"] = sample
    
    # Ajouter les AnnData au bon type cellulaire
    for cell_type in adata.obs['Annotations_old_with_unannotated'].unique():
        sub_adata = adata[adata.obs['Annotations_old_with_unannotated'] == cell_type].copy()
        if cell_type not in cell_type_dict:
            cell_type_dict[cell_type] = [sub_adata]
        else:
            cell_type_dict[cell_type].append(sub_adata)

# Dossier de sortie
output_dir = os.path.join(root_dir, "cell_type_combined")
os.makedirs(output_dir, exist_ok=True)

# Fusion et sauvegarde par type cellulaire

for cell_type, adata_list in cell_type_dict.items():
    combined = adata_list[0].concatenate(*adata_list[1:], batch_key="batch", batch_categories=None)
    
    safe_name = cell_type.replace(" ", "_").replace("/", "_")
    out_file = os.path.join(output_dir, f"{safe_name}_combined.h5ad")
    combined.write(out_file)
    print(f"✅ Sauvé : {out_file} ({combined.n_obs} spots)")
