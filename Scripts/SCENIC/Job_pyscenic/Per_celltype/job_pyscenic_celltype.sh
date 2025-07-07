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

# Define the cell types found in annotations
# Change these names according the annotation used
cell_types=(
    Adipose Blood Chondroid_tumour_cells Classical_chondrosarcoma_cells Epithelial_tumor_cells
    Immune_cells Intermediate_tumour_cells Mesenchymal_tumor_cells Mixed_cells Mixoid_chondrosarcoma_cells
    Mixoid_matrix-enriched_spindle+_spindle MMP9+_spindle MMP9-_spindle
    Necrosis Normal_epithelium Normal_fibrous_tissue NST_cells NST_surrounded_by_spindle
    Squamous_cell_tumour Osteosarcomatoid_tumour Pleiomorphic_tumour Scar-like_fibrous_stroma
    Spindle_cell_tumour Spindle_surrounded_by_NST
)

for cell_type in "${cell_types[@]}"; do
    echo "=== Processing $cell_type ==="

    # === Host and container paths ===
    # Change these paths according the setup
    host_project_dir="/home/dutelj/SCENIC"
    path_singularity_image="$host_project_dir/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif"
    container_mount_point="/data"
    host_output_dir="$host_project_dir/Outputs/$cell_type"
    container_output_dir="$container_mount_point/Outputs/$cell_type"

    # === Input files ===
    # Change these paths according the setup
    input_matrix_h5="$host_project_dir/${cell_type}.h5ad"

    # === Output files (HOST) ===
    # Change these paths according the setup
    host_raw_loom="$host_output_dir/Raw_matrix/${cell_type}_raw_matrix.loom"
    host_filtered_loom="$host_output_dir/Raw_matrix/${cell_type}_filtered_matrix.loom"
    host_adjacencies="$host_output_dir/Expr_matrix_adj/${cell_type}_adjacencies.tsv"
    host_regulons="$host_output_dir/Regulons/${cell_type}_regulons.csv"
    host_auc_matrix="$host_output_dir/AUCELL/${cell_type}_auc_matrix.csv"

    # === Output files (CONTAINER) ===
    # Change these paths according the setup
    container_raw_loom="$container_output_dir/Raw_matrix/${cell_type}_raw_matrix.loom"
    container_filtered_loom="$container_output_dir/Raw_matrix/${cell_type}_filtered_matrix.loom"
    container_adjacencies="$container_output_dir/Expr_matrix_adj/${cell_type}_adjacencies.tsv"
    container_regulons="$container_output_dir/Regulons/${cell_type}_regulons.csv"
    container_auc_matrix="$container_output_dir/AUCELL/${cell_type}_auc_matrix.csv"

    # === Environment setup ===
    # Change these paths according the setup
    path_conda="~/softs/miniconda3/etc/profile.d/conda.sh"
    environment_name="pyscenic"

    # === Create necessary directories ===
    mkdir -p "$(dirname "$host_raw_loom")"
    mkdir -p "$(dirname "$host_adjacencies")"
    mkdir -p "$(dirname "$host_regulons")"
    mkdir -p "$(dirname "$host_auc_matrix")"

    # === Step 1: Transpose expression matrix ===
    echo "----------------------------------"
    echo "$cell_type - Transposing expression matrix"
    echo "----------------------------------"
    source $path_conda
    conda activate $environment_name
    python transpose3.py -i "$input_matrix_h5" -o "$host_raw_loom"

    # === Step 1.5: Filter non-expressed genes ===
    echo "----------------------------------"
    echo "$cell_type - Filtering genes"
    echo "----------------------------------"
    python preprocess_to_loom.py -i "$host_raw_loom" -o "$host_filtered_loom"


    # === Step 2: GRN (Gene Regulatory Network inference) ===
    echo "---------------------"
    echo "$cell_type - STEP 1: GRN"
    echo "---------------------"

    echo 'singularity run -B "$host_project_dir:$container_mount_point" $path_singularity_image \
        pyscenic grn \
        "$container_filtered_loom" \
        "$container_mount_point/allTFs_hg38.txt" \
        --num_workers 30 \
        -o "$container_adjacencies"'
    singularity run -B "$host_project_dir:$container_mount_point" $path_singularity_image \
        pyscenic grn \
        "$container_filtered_loom" \
        "$container_mount_point/allTFs_hg38.txt" \
        --num_workers 30 \
        -o "$container_adjacencies"
        


    # === Step 3: Regulatory modules (CTX) ===
    echo "--------------------------------------------"
    echo "$cell_type - STEP 2: Regulatory modules (ctx)"
    echo "--------------------------------------------"

    echo 'singularity run -B "$host_project_dir:$container_mount_point" $path_singularity_image \
        pyscenic ctx \
        "$container_adjacencies" \
        "$container_mount_point/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        "$container_mount_point/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        --annotations_fname "$container_mount_point/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
        --expression_mtx_fname "$container_filtered_loom" \
        --mode "dask_multiprocessing" \
        --output "$container_regulons" \
        --num_workers 30'
    singularity run -B "$host_project_dir:$container_mount_point" $path_singularity_image \
        pyscenic ctx \
        "$container_adjacencies" \
        "$container_mount_point/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        "$container_mount_point/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
        --annotations_fname "$container_mount_point/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
        --expression_mtx_fname "$container_filtered_loom" \
        --mode "dask_multiprocessing" \
        --output "$container_regulons" \
        --num_workers 30


    # === Step 4: AUCELL (Regulon activity matrix) ===
    echo "-------------------------"
    echo "$cell_type - STEP 3: AUCELL"
    echo "-------------------------"

    echo 'singularity run -B "$host_project_dir:$container_mount_point" $path_singularity_image \
        pyscenic aucell \
        "$container_filtered_loom" \
        "$container_regulons" \
        --output "$container_auc_matrix" \
        --num_workers 30'
    singularity run -B "$host_project_dir:$container_mount_point" $path_singularity_image \
        pyscenic aucell \
        "$container_filtered_loom" \
        "$container_regulons" \
        --output "$container_auc_matrix" \
        --num_workers 30

    echo "=== SCENIC analysis $cell_type finished ==="






done











echo ""
echo ""
} >> "$LOGFILE" 2>&1

