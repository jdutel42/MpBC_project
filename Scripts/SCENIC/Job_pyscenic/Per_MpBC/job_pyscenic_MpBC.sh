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
        continue  # Do not process Mp_BC12 / MpBC12
    fi

    if [ "$i" -le 8 ]; then
        batch=1
        var="Mp_BC$i"
    else
        batch=2
        var="MpBC$i"
    fi
    echo "=== Processing $var (batch $batch) ==="

    # === Host and container paths ===
    host_project_dir="/home/dutelj/SCENIC"
    container_mount_point="/data"
    host_output_dir="$host_project_dir/Outputs/$sample_id"
    container_output_dir="$container_mount_point/Outputs/$sample_id"

    # === Input files ===
    input_matrix_h5="/home/martinep/data/visium/visium_mpbc_batch${batch_id}/${sample_id}/filtered_feature_bc_matrix.h5"

    # === Output files (HOST) ===
    host_raw_loom="$host_output_dir/Raw_matrix/${sample_id}_raw_matrix.loom"
    host_filtered_loom="$host_output_dir/Raw_matrix/${sample_id}_filtered_matrix.loom"
    host_adjacencies="$host_output_dir/Expr_matrix_adj/${sample_id}_adjacencies.tsv"
    host_regulons="$host_output_dir/Regulons/${sample_id}_regulons.csv"
    host_auc_matrix="$host_output_dir/AUCELL/${sample_id}_auc_matrix.csv"

    # === Output files (CONTAINER) ===
    container_raw_loom="$container_output_dir/Raw_matrix/${sample_id}_raw_matrix.loom"
    container_filtered_loom="$container_output_dir/Raw_matrix/${sample_id}_filtered_matrix.loom"
    container_adjacencies="$container_output_dir/Expr_matrix_adj/${sample_id}_adjacencies.tsv"
    container_regulons="$container_output_dir/Regulons/${sample_id}_regulons.csv"
    container_auc_matrix="$container_output_dir/AUCELL/${sample_id}_auc_matrix.csv"

    # === Create necessary directories ===
    mkdir -p "$(dirname "$host_raw_loom")"
    mkdir -p "$(dirname "$host_adjacencies")"
    mkdir -p "$(dirname "$host_regulons")"
    mkdir -p "$(dirname "$host_auc_matrix")"

    # === Step 1: Transpose expression matrix ===
    echo "----------------------------------"
    echo "$var - Transposing expression matrix"
    echo "----------------------------------"
    source ~/softs/miniconda3/etc/profile.d/conda.sh
    conda activate basic
    python transpose2.py -i "$input_matrix_h5" -o "$host_raw_loom"

    # === Step 1.5: Filter non-expressed genes ===
    echo "----------------------------------"
    echo "$var - Filtering genes"
    echo "----------------------------------"
    python preprocess_to_loom.py -i "$host_raw_loom" -o "$host_filtered_loom"


    # === Step 2: GRN (Gene Regulatory Network inference) ===
    echo "---------------------"
    echo "$var - STEP 1: GRN"
    echo "---------------------"

    echo 'singularity run -B "$host_project_dir:$container_mount_point" aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
        pyscenic grn \
        "$container_filtered_loom" \
        "$container_mount_point/allTFs_hg38.txt" \
        --num_workers 18 \
        -o "$container_adjacencies"'
    singularity run -B "$host_project_dir:$container_mount_point" aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
        pyscenic grn \
        "$container_filtered_loom" \
        "$container_mount_point/allTFs_hg38.txt" \
        --num_workers 18 \
        -o "$container_adjacencies"
        


    # === Step 3: Regulatory modules (CTX) ===
    echo "--------------------------------------------"
    echo "$var - STEP 2: Regulatory modules (ctx)"
    echo "--------------------------------------------"

    echo 'singularity run -B "$host_project_dir:$container_mount_point" aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
        pyscenic ctx \
        "$container_adjacencies" \
        "$container_mount_point/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        "$container_mount_point/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        --annotations_fname "$container_mount_point/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
        --expression_mtx_fname "$container_filtered_loom" \
        --mode "dask_multiprocessing" \
        --output "$container_regulons" \
        --num_workers 18'
    singularity run -B "$host_project_dir:$container_mount_point" aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
        pyscenic ctx \
        "$container_adjacencies" \
        "$container_mount_point/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        "$container_mount_point/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        --annotations_fname "$container_mount_point/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
        --expression_mtx_fname "$container_filtered_loom" \
        --mode "dask_multiprocessing" \
        --output "$container_regulons" \
        --num_workers 18


    # === Step 4: AUCELL (Regulon activity matrix) ===
    echo "-------------------------"
    echo "$var - STEP 3: AUCELL"
    echo "-------------------------"

    echo 'singularity run -B "$host_project_dir:$container_mount_point" aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
        pyscenic aucell \
        "$container_filtered_loom" \
        "$container_regulons" \
        --output "$container_auc_matrix" \
        --num_workers 18'
    singularity run -B "$host_project_dir:$container_mount_point" aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
        pyscenic aucell \
        "$container_filtered_loom" \
        "$container_regulons" \
        --output "$container_auc_matrix" \
        --num_workers 18

    echo "=== SCENIC analysis $var finished ==="











done












echo ""
echo ""
} >> "$LOGFILE" 2>&1

