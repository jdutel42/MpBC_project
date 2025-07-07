import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
import math
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss

# Liste des échantillons
samples = ['Mp_BC1', 'Mp_BC2', 'Mp_BC3', 'Mp_BC4', 'Mp_BC5',
           'Mp_BC6', 'Mp_BC7', 'Mp_BC8', 'Mp_BC9', 'Mp_BC10',
           'Mp_BC11', 'Mp_BC13', 'Mp_BC14', 'Mp_BC15', 'Mp_BC16']

output_dir = '/mnt/datadisk/Jordan/Results/SCENIC/Outputs_per_MpBC'

# Initialiser liste pour stocker les top 5 de tous les samples
all_top5_list = []

# Initialiser liste pour stocker tous les régulons
all_list = []

for sample in samples:
    print(f"Processing sample: {sample}")

    # Load AUC matrix and annotations
    auc_mtx_path = f'{output_dir}/{sample}/AUCELL/auc_mtx_{sample}.csv'
    annot_path = f'{output_dir}/{sample}/Annotations/Annotations_old_with_unannotated.csv'

    auc_mtx = pd.read_csv(auc_mtx_path, index_col=0)
    cellAnnot = pd.read_csv(annot_path, index_col=0)
    if cellAnnot.shape[1] == 1:
        cellAnnot.columns = ['CellType']

    # Filtrer uniquement les cellules annotées
    common_cells = auc_mtx.index.intersection(cellAnnot.index)
    cellAnnot = cellAnnot.loc[common_cells]
    annotated_cells = cellAnnot[cellAnnot['CellType'] != 'Unannotated'].index
    auc_mtx = auc_mtx.loc[annotated_cells]
    cellAnnot = cellAnnot.loc[annotated_cells]

    # Calculer les RSS
    rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot['CellType'])

    # Extraire les types cellulaires uniques triés
    cats = sorted(cellAnnot['CellType'].unique())

    # ----------- PLOTS RSS ------------
    sns.set_style("whitegrid")
    n_cats = len(cats)
    ncols = 4
    nrows = math.ceil(n_cats / ncols)
    fig = plt.figure(figsize=(ncols * 4, nrows * 4))
    fig.suptitle("Top 5 Regulon Specificity Scores (RSS) by Cell Type", fontsize=18, y=1)

    for c, num in zip(cats, range(1, n_cats + 1)):
        x = rss_cellType.T[c]
        ax = fig.add_subplot(nrows, ncols, num)
        plot_rss(rss_cellType, c, top_n=5, ax=ax)
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

    fig.subplots_adjust(top=0.88, bottom=0.07, left=0.06, right=0.97, hspace=0.4, wspace=0.3)
    plt.savefig(f'{output_dir}/{sample}/RSS/rss_{sample}.png', dpi=300)
    plt.close()

    # ----------- CSV TOUS LES RÉGULONS ------------
    full_rss_df = rss_cellType.T.reset_index().melt(id_vars='index', var_name='CellType', value_name='RSS')
    full_rss_df = full_rss_df.rename(columns={'index': 'Regulon'})
    full_rss_df = full_rss_df.sort_values(by=['CellType', 'RSS'], ascending=[True, False])
    full_rss_df.to_csv(f'{output_dir}/{sample}/RSS/rss_sorted_{sample}.csv', index=False)

    # ----------- CSV TOP 5 RÉGULONS PAR CELLTYPE ------------
    top5_rows = []
    for c in cats:
        top5 = rss_cellType.T[c].sort_values(ascending=False).head(5)
        for rank, (reg, rss) in enumerate(top5.items(), start=1):
            top5_rows.append({
                'Sample': sample,
                'CellType': c,
                'Rank': rank,
                'Regulon': reg,
                'RSS': rss
            })

    # Ajouter à la liste globale
    all_top5_list.extend(top5_rows)

    # ----------- CSV TOUS LES RÉGULONS PAR CELLTYPE ------------
    all_rows = []
    for c in cats:
        for reg, rss in rss_cellType.T[c].sort_values(ascending=False).items():
            all_rows.append({
                'Sample': sample,
                'CellType': c,
                'Regulon': reg,
                'RSS': rss
            })

    # Ajouter à la liste globale
    all_list.extend(all_rows)

# ----------- FICHIER GLOBAL TOP 5 POUR TOUS LES ÉCHANTILLONS ------------
df_all_top5 = pd.DataFrame(all_top5_list)
df_all_top5 = df_all_top5.sort_values(by=['Sample', 'CellType', 'Rank'])
df_all_top5.to_csv(f'{output_dir}/all_samples_top5_RSS.csv', index=False)

# ----------- FICHIER GLOBAL POUR TOUS LES ÉCHANTILLONS ------------
df_all = pd.DataFrame(all_list)
df_all = df_all.sort_values(by=['Sample', 'CellType'])
df_all.to_csv(f'{output_dir}/all_samples_RSS.csv', index=False)

print("✔️ Tous les fichiers ont été générés avec succès.")
# ----------- FIN DU SCRIPT ------------