import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table
from textwrap import fill

# Charger les données TSV
df = pd.read_csv("~/Downloads/test.csv", sep=",")

# Sous-ensemble
df_head = df.iloc[:10, :10].copy()

# Fonction de wrap texte
def wrap_text(s, width=20):
    try:
        return fill(str(s), width=width)
    except:
        return s

# Appliquer le wrap
df_wrapped = df_head.applymap(lambda x: wrap_text(x, width=20))
df_wrapped.columns = [wrap_text(col, width=20) for col in df_wrapped.columns]

# Estimer le nombre moyen de lignes par cellule pour ajuster la hauteur
avg_lines_per_cell = df_wrapped.applymap(lambda x: str(x).count("\n") + 1).values.mean()

# Taille de figure
cell_height = 0.5  # ← augmente pour plus de hauteur
cell_width = 1.5
fig_height = 1 + cell_height * avg_lines_per_cell * df_wrapped.shape[0]
fig_width = 1 + cell_width * df_wrapped.shape[1]

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
ax.axis("off")

# Tracer le tableau
tbl = table(ax, df_wrapped, loc="center", cellLoc="center", colWidths=[0.1]*df_wrapped.shape[1])
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)

# Ajustement manuel de l’échelle
tbl.scale(1.2, 4.0)  # ← Augmente le second paramètre pour forcer plus de hauteur

# Sauvegarde
plt.tight_layout()
plt.savefig("aperçu_tsv.png", dpi=300, bbox_inches="tight")
plt.show()
