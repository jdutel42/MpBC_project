import argparse
import scanpy as sc

def main():
    parser = argparse.ArgumentParser(description="Transpose a 10X HDF5 gene expression matrix.")
    parser.add_argument("--input", "-i", required=True, help="Path to the input 10X HDF5 file")
    parser.add_argument("--output", "-o", required=True, help="Path to the output transposed CSV file")
    args = parser.parse_args()

    # Lire le fichier 10X HDF5
    print(f"Lecture du fichier H5 : {args.input}")
    adata = sc.read_10x_h5(args.input)

    # Make uniques names for variables
    adata.var_names_make_unique()

    # Extraire la matrice en DataFrame
    df = adata.to_df()

    # Transposer
    df_T = df.T
    # df_T = df


    # Sauvegarder en CSV
    df_T.to_csv(args.output, sep="\t")
    print(f"✅ Fichier transposé sauvegardé : {args.output}")

if __name__ == "__main__":
    main()
