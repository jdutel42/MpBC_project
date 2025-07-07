# import pandas as pd
# import matplotlib.pyplot as plt
# from pandas.plotting import table

# # Charger le TSV
# df = pd.read_csv("./Outputs/Mp_BC1/Regulons/regulons_Mp_BC1.csv", sep=",")

# # # Sélectionner les premières lignes
# # n_lignes = 10
# # df_head = df.head(n_lignes)

# # Garder les 10 premières lignes et 5 premières colonnes
# df_head = df.iloc[:10, :5]

# # Créer une figure
# fig, ax = plt.subplots(figsize=(12, 0.5 * n_lignes))  # Ajuste la taille selon le nombre de lignes
# ax.axis("off")  # Supprime les axes

# # Ajouter le tableau
# tbl = table(ax, df_head, loc='center', cellLoc='center', colWidths=[0.2]*len(df_head.columns))

# # Personnaliser l'apparence
# tbl.auto_set_font_size(False)
# tbl.set_fontsize(10)
# tbl.scale(1.2, 1.2)  # Échelle du tableau (x, y)

# # Sauvegarder l’image
# plt.savefig("aperçu_tsv.png", bbox_inches='tight', dpi=300)
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table

# Charger le fichier CSV
df = pd.read_csv("./Outputs/Mp_BC1/RSS/top_regulons_Mp_BC1.csv", sep=",")

# Garder les 10 premières lignes et 5 premières colonnes
df_head = df.iloc[:10, :6]

# Créer une figure ajustée
fig, ax = plt.subplots(figsize=(10, 0.6 * len(df_head)))  # Largeur adaptée à 5 colonnes
ax.axis("off")

# Ajouter le tableau
tbl = table(ax, df_head, loc='center', cellLoc='center', colWidths=[0.2]*df_head.shape[1])

# Personnalisation
tbl.auto_set_font_size(False)
tbl.set_fontsize(10)
tbl.scale(1.2, 1.2)

# Sauvegarde
plt.savefig("aperçu_tsv.png", bbox_inches='tight', dpi=300)
plt.show()

