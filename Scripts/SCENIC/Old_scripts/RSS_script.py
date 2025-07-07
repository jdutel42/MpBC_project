from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
import math
import pandas as pd
import loompy as lp

# samples = [f'Mp_BC{i}' for i in range(1, 17)]  # Mp_BC1 to Mp_BC16, skipping Mp_BC12
samples = ['Mp_BC1', 'Mp_BC2', 'Mp_BC3', 'Mp_BC4', 'Mp_BC5',
           'Mp_BC6', 'Mp_BC7', 'Mp_BC8', 'Mp_BC9', 'Mp_BC10',
           'Mp_BC11', 'Mp_BC13', 'Mp_BC14', 'Mp_BC15', 'Mp_BC16']

output_dir = '/mnt/datadisk/Jordan/Results/SCENIC/Outputs_per_MpBC'

for sample in samples:
    if sample == 'Mp_BC12':
        print(f"Skipping sample: {sample}")
        continue

    print(f"Processing sample: {sample}")

    # Load the AUC matrix and cell annotations
    auc_mtx_path = f'{output_dir}/{sample}/AUCELL/auc_mtx_{sample}.csv'
    annot_path = f'{output_dir}/{sample}/Annotations/Annotations_old_with_unannotated.csv'

    # REad the AUC matrix and cell annotations
    auc_mtx = pd.read_csv(auc_mtx_path, index_col=0)
    cellAnnot = pd.read_csv(annot_path, index_col=0)
    if cellAnnot.shape[1] == 1:
        cellAnnot.columns = ['CellType']

    # Filter the AUC matrix and cell annotations to keep only common cells
    ## Some pots xere remove from AUCell matrix since filtering before SCENIC
    common_cells = auc_mtx.index.intersection(cellAnnot.index)
    cellAnnot = cellAnnot.loc[common_cells]

    annotated_cells = cellAnnot[cellAnnot['CellType'] != 'Unannotated'].index
    auc_mtx = auc_mtx.loc[annotated_cells,:]
    cellAnnot = cellAnnot.loc[annotated_cells]


    # Load the regulon specificity scores
    rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot['CellType'])

    # Plot the regulon specificity scores
    # Get sorted list of unique cell types (categories)
    cats = sorted(list(set(cellAnnot['CellType'])))

    # Set seaborn style for the plots
    sns.set_style("whitegrid")

    n_cats = len(cats)
    ncols = 4
    nrows = math.ceil(n_cats / ncols)

    # Create a figure with subplots for each cell type
    fig = plt.figure(figsize=(ncols * 4, nrows * 4))
    fig.suptitle("Top 5 Regulon Specificity Scores (RSS) by Cell Type", fontsize=18, y=1)

    # Loop through each cell type and create a subplot
    for c, num in zip(cats, range(1, n_cats + 1)):
        x = rss_cellType.T[c]
        ax = fig.add_subplot(nrows, ncols, num)
        plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
        ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05,
                    x.max() + (x.max() - x.min()) * 0.05)
        for t in ax.texts:
            t.set_fontsize(12)
            t.set_color("darkred")
        ax.set_title(c, fontsize=14, fontweight='bold')
        ax.set_ylabel('')
        ax.set_xlabel('')
        adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom',
                    arrowprops=dict(arrowstyle='-', color='lightgrey'), precision=0.001)

    # Adjust spacing to prevent overlap
    fig.subplots_adjust(top=0.88, bottom=0.07, left=0.06, right=0.97, hspace=0.4, wspace=0.3)

    # Save figure
    plt.savefig(f'{output_dir}/{sample}/RSS/rss_{sample}.png', dpi=300)
    plt.close()


    # Store in cats, unique value from cellAnnot['CellType']
    cats = cellAnnot['CellType'].unique().tolist()

    topreg_dict = {}

    # For each cell type, get the top 5 regulons based on RSS
    for c in cats:
        top5 = rss_cellType.T[c].sort_values(ascending=False)[:5]
        topreg_dict[c] = list(top5.index)

    # Create a DataFrame from the top regulons dictionary
    topreg_df = pd.DataFrame.from_dict(
        topreg_dict, orient='index', columns=[f'Rang {i+1}' for i in range(5)]
    )

    # Save the top regulons DataFrame to a CSV file
    topreg_df.to_csv(f'{output_dir}/{sample}/RSS/top_regulons_{sample}.csv')

    print(f"Top regulons for {sample} saved to top_regulons_{sample}.csv")