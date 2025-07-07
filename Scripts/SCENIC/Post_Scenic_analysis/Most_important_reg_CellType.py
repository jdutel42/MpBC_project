import pandas as pd
import numpy as np
import os
import glob

# Script to analyze AUCell matrices for different cell types and extract the most important regulons.

cell_types = ['Adipose', 'Blood', 'Chondroid_tumour_cells', 'Classical_chondrosarcoma_cells',
              'Epithelial_tumor_cells', 'Immune_cells', 'Intermediate_tumour_cells',
              'Mesenchymal_tumor_cells', 'Mixed_cells', 'Mixoid_chondrosarcoma_cells',
              'Mixoid_matrix-enriched_spindle+_spindle', 'MMP9+_spindle',
              'MMP9-_spindle', 'Necrosis', 'Normal_epithelium', 'Normal_fibrous_tissue',
              'NST_cells', 'NST_surrounded_by_spindle', 'Squamous_cell_tumour',
              'Osteosarcomatoid_tumour', 'Pleiomorphic_tumour',
              'Scar-like_fibrous_stroma', 'Spindle_cell_tumour',
              'Spindle_surrounded_by_NST']

# === PARAMETERS ===
INPUT_FOLDER = f"/mnt/datadisk/Jordan/Results/SCENIC/Outputs_per_CellType/Annotations_old/SCENIC_pipeline_outputs/"
OUTPUT_FOLDER = "/mnt/datadisk/Jordan/Results/SCENIC/Outputs_per_CellType/Annotations_old/Post_SCENIC_outputs/Top_regulons/"
TOP_N = 10                                # top regulons to extract
BIN_THRESHOLD = 0.80                      # treshold quantile for binarization

# === CREATE OUTPUT FOLDER ===
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# === FOR EACH CELL TYPE ===
for cell_type in cell_types:
    file_path = f'{INPUT_FOLDER}/{cell_type}/AUCELL/{cell_type}_auc_matrix.csv'

    print(f"\n🔍 Traitement de : {cell_type}")
    
    # Charge the AUC matrix for the current cell type (cellules x régulons)
    auc = pd.read_csv(file_path, index_col=0)

    # Binarization of the AUC matrix based on a quantile threshold
    thresholds = auc.quantile(BIN_THRESHOLD)
    binary_auc = auc.gt(thresholds).astype(int)

    # Mean AUC for each regulon
    mean_auc = auc.mean().sort_values(ascending=False)
    top_mean_auc = mean_auc.head(TOP_N)

    # % of active cells for each regulon
    freq_active = binary_auc.sum() / binary_auc.shape[0] * 100
    top_freq_active = freq_active.sort_values(ascending=False).head(TOP_N)

    # SUmmary
    top_regulons = pd.DataFrame({
        "Regulon": top_mean_auc.index,
        "Mean AUC activity": top_mean_auc.values,
        "% Cells Actives": [freq_active.get(reg, np.nan) for reg in top_mean_auc.index]
    })


    # Save
    top_regulons.to_csv(f"{OUTPUT_FOLDER}/regulons_{cell_type}.csv", index=False)
    print(top_regulons)
