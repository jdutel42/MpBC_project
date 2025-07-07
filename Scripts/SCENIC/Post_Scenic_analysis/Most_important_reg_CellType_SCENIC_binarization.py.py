import os
import pandas as pd
from pyscenic.binarization import binarize
# from pyscenic.utils import load_regulons
import numpy as np
import warnings

# Suppress deprecation warnings (e.g., numpy float)
warnings.filterwarnings("ignore")
if not hasattr(np, 'float'):
    np.float = float

# List of cell types to process
cell_types = [
    'Adipose', 'Blood', 'Chondroid_tumour_cells', 'Classical_chondrosarcoma_cells',
    'Epithelial_tumor_cells', 'Immune_cells', 'Intermediate_tumour_cells',
    'Mesenchymal_tumor_cells', 'Mixed_cells', 'Mixoid_chondrosarcoma_cells',
    'Mixoid_matrix-enriched_spindle+_spindle', 'MMP9+_spindle',
    'MMP9-_spindle', 'Necrosis', 'Normal_epithelium', 'Normal_fibrous_tissue',
    'NST_cells', 'NST_surrounded_by_spindle', 'Squamous_cell_tumour',
    'Osteosarcomatoid_tumour', 'Pleiomorphic_tumour',
    'Scar-like_fibrous_stroma', 'Spindle_cell_tumour',
    'Spindle_surrounded_by_NST'
]

# Input/output directories (adjust as needed !!!!!!!!)
input_dir = "/data/home/jdutel/pyscenic/"
output_dir = os.path.join(input_dir, "binarized_outputs")
os.makedirs(output_dir, exist_ok=True)

print("🚀 Starting SCENIC regulon binarization process...\n")

for i, cell_type in enumerate(cell_types, 1):
    print(f"🔹 [{i}/{len(cell_types)}] Processing cell type: {cell_type}")

    regulon_file = os.path.join(input_dir, f"Outputs/{cell_type}/Regulons/{cell_type}_regulons.csv")
    input_auc_path = os.path.join(input_dir, f"Outputs/{cell_type}/AUCELL/{cell_type}_auc_matrix.csv")
    output_celltype_dir = os.path.join(output_dir, cell_type)
    os.makedirs(output_celltype_dir, exist_ok=True)

    output_bin_path = os.path.join(output_celltype_dir, f"{cell_type}_binarized.csv")
    output_topreg_path = os.path.join(output_celltype_dir, f"{cell_type}_top10_regulons.csv")
    output_allreg_path = os.path.join(output_celltype_dir, f"{cell_type}_all_regulons.csv")

    # Load input matrices
    try:
        print("  📥 Loading AUC and regulon files...")
        auc_mtx = pd.read_csv(input_auc_path, index_col=0)
        regulons = pd.read_csv(regulon_file)

    except Exception as e:
        print(f"  ❌ Failed to load files for {cell_type}: {e}")
        continue

    # Binarize AUC matrix using regulons
    try:
        print("  🧠 Binarizing regulon activity...")
        binary_mtx, auc_thresholds = binarize(auc_mtx, num_workers=20)
        binary_mtx.to_csv(output_bin_path)
        thresholds_path = os.path.join(output_celltype_dir, "thresholds_per_regulon.csv")
        auc_thresholds.to_csv(thresholds_path)
        print(f"  ✅ Binarized matrix saved at: {output_bin_path}")
    except Exception as e:
        print(f"  ❌ Binarization failed for {cell_type}: {e}")
        continue

    # Calculate top active regulons
    try:
        print("  📊 Calculating top active regulons...")
        freq_active = binary_mtx.sum(axis=0) / binary_mtx.shape[0] * 100  # percentage of active cells (meaning cells where the regulon is active (=1)) /  It's a mean, not median
        top_regulons = freq_active.sort_values(ascending=False).head(10)
        all_regulons = freq_active.sort_values(ascending=False)

        print(f"  🏆 Top 10 regulons for {cell_type}:")
        print(top_regulons)
        print("  ---------------------------------------------")

        # Save top regulons
        top_regulons.to_frame(name="Percent_cells/spots_active").to_csv(output_topreg_path)
        print(f"  💾 Top regulons saved at: {output_topreg_path}\n")
        # Save all regulons
        all_regulons.to_frame(name="Percent_cells_active").to_csv(output_allreg_path)
        print(f"  💾 All regulons saved at: {output_allreg_path}\n")
    except Exception as e:
        print(f"  ❌ Failed to analyze/save regulons for {cell_type}: {e}")

print("✅ SCENIC regulon binarization completed for all cell types.")
