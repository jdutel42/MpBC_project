import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# # Charger les données
# df = pd.read_csv("./Outputs/all_samples_RSS.csv")

# # S'assurer que les colonnes sont bien nommées
# # Vérifie que les colonnes sont bien : "sample", "cell_type", "regulon"
# print(df.columns)

# # Créer un dossier de sortie pour les graphiques
# output_dir = "rss_frequency_plots"
# os.makedirs(output_dir, exist_ok=True)

# # Boucle par type cellulaire
# for cell_type, group in df.groupby("CellType"):
#     # Compter la fréquence des régulons pour ce type cellulaire
#     freq = group["Regulon"].value_counts().sort_values(ascending=False)

#     # Créer la figure
#     plt.figure(figsize=(10, max(4, 0.4 * len(freq))))
#     sns.barplot(x=freq.values, y=freq.index, palette="viridis")
#     plt.title(f"Fréquence des régulons pour le type cellulaire : {cell_type}")
#     plt.xlabel("Nombre d'apparitions dans le top 5")
#     plt.ylabel("Régulons (RSS)")
#     plt.tight_layout()

#     # Sauvegarder le graphique
#     plt.savefig(f"{output_dir}/{cell_type.replace(' ', '_')}_rss_frequency.png", dpi=300)
#     plt.close()

# print(f"✅ Graphiques générés dans le dossier : {output_dir}/")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("./Outputs/all_samples_RSS.csv")

# Compter les occurrences
count_df = df.groupby(["CellType", "Regulon"]).size().reset_index(name="frequency")

# Plot
plt.figure(figsize=(15, 8))
sns.scatterplot(data=count_df, x="Regulon", y="CellType", size="frequency", hue="frequency", palette="coolwarm", legend="brief", sizes=(50, 300))

plt.title("Fréquence des régulons par type cellulaire (dot plot)", fontsize=14)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("rss_dotplot.png", dpi=300)
plt.show()


