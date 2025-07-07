import scanpy as sc
import numpy as np
import loompy
import argparse

# === Arguments ===
parser = argparse.ArgumentParser(description="Filter and normalize a LOOM file")
parser.add_argument("-i", "--input", required=True, help="Input loom file")
parser.add_argument("-o", "--output", required=True, help="Output filtered loom file")
args = parser.parse_args()

input_loom = args.input
output_loom = args.output

# === Lecture des données ===
adata = sc.read_loom(input_loom)
adata.var_names_make_unique()


# === Filtrage qualité minimale ===
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# === Pourcentage de mitochondries ===
if isinstance(adata.X, np.ndarray):
    mito_sum = np.sum(adata[:, adata.var_names.str.startswith("MT-")].X, axis=1)
    total_sum = np.sum(adata.X, axis=1)
else:
    mito_sum = np.sum(adata[:, adata.var_names.str.startswith("MT-")].X, axis=1).A1
    total_sum = np.sum(adata.X, axis=1).A1

adata.obs["percent_mito"] = mito_sum / total_sum

# === Normalisation + log1p ===
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

# === Sélection des gènes hautement variables ===
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var["highly_variable"]]

# === Export en loom ===
row_attrs = {
    "Gene": np.array(adata.var_names)
}
col_attrs = {
    "CellID": np.array(adata.obs_names)
}
loompy.create(output_loom, adata.X.transpose(), row_attrs, col_attrs)

print(f"Fichier loom filtré sauvegardé dans : {output_loom}")
