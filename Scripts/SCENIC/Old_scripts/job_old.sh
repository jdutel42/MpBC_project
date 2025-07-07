#!/bin/bash
#
#SBATCH --job-name=pyscenic
#SBATCH --output=pyscenic_out
#SBATCH --partition=research
#SBATCH --mem=40G
#
#SBATCH --ntasks=10
#SBATCH --time=10:00:00



set -e

for i in {1..16}; do
    if [ "$i" -eq 12 ]; then
        continue  # On saute Mp_BC12
    fi

    var="Mp_BC$i"

    if [ "$i" -le 8 ]; then
        batch=1
    else
        batch=2
    fi

    echo "=== Traitement de $var (batch $batch) ==="

    # === Chemins ===
    input_file_raw_MpBC_matrix="/home/martinep/data/visium/visium_mpbc_batch${batch}/${var}/filtered_feature_bc_matrix.h5"
    output_dir="/home/dutelj/SCENIC/Outputs/$var"
    raw_matrix_csv="$output_dir/Raw_matrix/Raw_${var}_matrix_transpose.csv"
    adjacencies_file="$output_dir/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv"
    regulons_file="$output_dir/Regulons/regulons_${var}.csv"
    aucell_file="$output_dir/AUCELL/auc_mtx_${var}.csv"

    # === Création des dossiers ===
    mkdir -p "$output_dir/Raw_matrix"
    mkdir -p "$output_dir/Expr_matrix_adj"
    mkdir -p "$output_dir/Regulons"
    mkdir -p "$output_dir/AUCELL"

    # === Étape 1 : Transposition ===
    echo "$var - Transposition matrice"
    source ~/softs/miniconda3/etc/profile.d/conda.sh
    conda activate basic
    python transpose.py --input "$input_file_raw_MpBC_matrix" --output "$raw_matrix_csv"

    # === Étape 2 : GRN ===
    echo "$var - STEP 1: GRN"
    singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic grn \
        "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.csv" \
        /data/allTFs_hg38.txt \
        --num_workers 6 \
        -o "/data/Outputs/$var/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv"

    # === Étape 3 : Contextualisation (ctx) ===
    echo "$var - STEP 2: Contextualisation (ctx)"
    singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic ctx \
        "/data/Outputs/$var/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv" \
        /data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname /data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.csv" \
        --mode "dask_multiprocessing" \
        --output "/data/Outputs/$var/Regulons/regulons_${var}.csv" \
        --num_workers 6

    # === Étape 4 : AUCELL ===
    echo "$var - STEP 3: AUCELL"
    singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic aucell \
        "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.csv" \
        "/data/Outputs/$var/Regulons/regulons_${var}.csv" \
        --output "/data/Outputs/$var/AUCELL/auc_mtx_${var}.csv" \
        --num_workers 6

    echo "=== $var terminé ==="
done
