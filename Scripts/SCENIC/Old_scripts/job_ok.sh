#!/bin/bash
#
#SBATCH --job-name=pyscenic
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --partition=research
#
#SBATCH --mem=60G
#SBATCH --ntasks=15
#SBATCH --time=24:00:00

# Variables
LOGFILE="logs/job_errors_all.log"
DATE=$(date "+%Y-%m-%d %H:%M:%S")
JID=$SLURM_JOB_ID
NAME=$SLURM_JOB_NAME

# Encapsule ton code dans une redirection stderr
{
  echo ""
  echo ""
  echo "===================================================================================================="
  echo "====================[ $DATE | JobID: $JID | JobName: $NAME ]===================="
  echo "===================================================================================================="







# 20min per sample

set -e

for i in {1..16}; do
    if [ "$i" -eq 12 ]; then
        continue  # On saute Mp_BC12 / MpBC12
    fi

    if [ "$i" -le 8 ]; then
        batch=1
        var="Mp_BC$i"
    else
        batch=2
        var="MpBC$i"
    fi
    echo "=== Traitement de $var (batch $batch) ==="

    # === Chemins ===
    input_file_raw_MpBC_matrix="/home/martinep/data/visium/visium_mpbc_batch${batch}/${var}/filtered_feature_bc_matrix.h5"
    output_dir="/home/dutelj/SCENIC/Outputs/$var"
    raw_matrix_csv="$output_dir/Raw_matrix/Raw_${var}_matrix_transpose.csv"
    raw_matrix_loom="$output_dir/Raw_matrix/Raw_${var}_matrix_transpose.loom"
    adjacencies_file="$output_dir/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv"
    regulons_file="$output_dir/Regulons/regulons_${var}.csv"
    aucell_file="$output_dir/AUCELL/auc_mtx_${var}.csv"

    # === Création des dossiers ===
    mkdir -p "$output_dir/Raw_matrix"
    mkdir -p "$output_dir/Expr_matrix_adj"
    mkdir -p "$output_dir/Regulons"
    mkdir -p "$output_dir/AUCELL"

    # === Étape 1 : Transposition ===
    echo "----------------------------------"
    echo "$var - Transposition matrice"
    echo "----------------------------------"
    source ~/softs/miniconda3/etc/profile.d/conda.sh
    conda activate basic
    python transpose2.py -i "$input_file_raw_MpBC_matrix" -o "$raw_matrix_loom"

    # === Étape 1.5 : Filtrage des gènes non exprimés ===
    echo "----------------------------------"
    echo "$var - Filtrage des gènes"
    echo "----------------------------------"
    filtered_loom="$output_dir/Raw_matrix/Raw_${var}_matrix_transpose.filtered.loom"
    python preprocess_to_loom.py -i "$raw_matrix_loom" -o "$filtered_loom"


    # === Étape 2 : GRN ===
    echo "---------------------"
    echo "$var - STEP 1: GRN"
    echo "---------------------"

    echo 'singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic grn \
        "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.loom" \
        /data/allTFs_hg38.txt \
        --num_workers 18 \
        -o "/data/Outputs/$var/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv"'
    singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic grn \
        "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.filtered.loom" \
        /data/allTFs_hg38.txt \
        --num_workers 18 \
        -o "/data/Outputs/$var/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv"

    # === Étape 3 : Contextualisation (ctx) ===
    echo "--------------------------------------------"
    echo "$var - STEP 2: Contextualisation (ctx)"
    echo "--------------------------------------------"

    echo 'singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic ctx \
        "/data/Outputs/$var/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv" \
        /data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname /data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.csv" \
        --mode "dask_multiprocessing" \
        --output "/data/Outputs/$var/Regulons/regulons_${var}.csv" \
        --num_workers 18'
    singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic ctx \
        "/data/Outputs/$var/Expr_matrix_adj/expr_mat_${var}.adjacencies.tsv" \
        /data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname /data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.filtered.loom" \
        --mode "dask_multiprocessing" \
        --output "/data/Outputs/$var/Regulons/regulons_${var}.csv" \
        --num_workers 18

    # === Étape 4 : AUCELL ===
    echo "-------------------------"
    echo "$var - STEP 3: AUCELL"
    echo "-------------------------"

    echo 'singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic aucell \
        "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.loom" \
        "/data/Outputs/$var/Regulons/regulons_${var}.csv" \
        --output "/data/Outputs/$var/AUCELL/auc_mtx_${var}.csv" \
        --num_workers 18'
    singularity run -B /home/dutelj/SCENIC:/data aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
      pyscenic aucell \
        "/data/Outputs/$var/Raw_matrix/Raw_${var}_matrix_transpose.filtered.loom" \
        "/data/Outputs/$var/Regulons/regulons_${var}.csv" \
        --output "/data/Outputs/$var/AUCELL/auc_mtx_${var}.csv" \
        --num_workers 18

    echo "=== Analyse SCENIC $var terminé ==="











done












echo ""
echo ""
} >> "$LOGFILE" 2>&1

