import scanpy as sc
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input .h5 file (10X format)")
    parser.add_argument("-o", "--output", required=True, help="Output .loom file")
    args = parser.parse_args()

    # Lire le fichier HDF5 avec le bon chemin
    adata = sc.read_10x_h5(args.input)

    # Rendre les noms de gènes uniques
    adata.var_names_make_unique()

    # Ajouter l'attribut requis pour SCENIC
    adata.var["Gene"] = adata.var_names
    adata.obs["CellID"] = adata.obs_names

    # Sauvegarde en loom
    adata.write_loom(args.output)

if __name__ == "__main__":
    main()
