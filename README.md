# Stage M2 : Transcriptomique spatiale et évolution des tumeurs mammaires rares

Bienvenue dans ce dépôt dédié au stage de Master 2 portant sur la transcriptomique spatiale et l’évolution des tumeurs mammaires rares. Ce projet, encadré par **Dr Pierre Martinez**, se déroule au Centre de Recherche en Cancérologie de Lyon. Vous trouverez ci-dessous toutes les informations essentielles concernant le stage.
[Sujet_de_stage_M2](/Sujet_stage_M2.pdf)

---

## Contexte et Objectifs

### Contexte
Certaines tumeurs du sein présentent une transdifférenciation où les cellules épithéliales d’origine deviennent mésenchymateuses. Ces formes rares (1-3%) demeurent mal comprises et mal prises en charge en clinique. La plasticité cellulaire représente également un défi thérapeutique majeur.

### Objectifs du Stage
Le projet de stage vise à :
- **Identifier les gènes et pathways** impliqués dans la transdifférenciation.
- **Analyser les différences** génomiques et micro-environnementales entre les compartiments épithéliaux et mésenchymateux.
- **Déterminer des marqueurs spécifiques** des cellules tumorales mésenchymateuses.

Les analyses se baseront sur des données issues de coupes d’échantillons, avec notamment l’utilisation de la technologie Visium de 10X Genomics.

---

## Informations Pratiques

- **Encadrement :** Dr Pierre Martinez (bioinformaticien, Inserm) – [pierre.martinez@lyon.unicancer.fr](mailto:pierre.martinez@lyon.unicancer.fr)
- **Lieu :** Centre de Recherche en Cancérologie de Lyon, Cheney D 2e étage.
- **GitGub de l'encadrant :** [https://pierremartinez.github.io/](https://pierremartinez.github.io/)

---

## Références Bibliographiques

1. **McCart Reed et al. (2019)**  
   *Phenotypic and molecular dissection of metaplastic breast cancer and the prognostic implications.*
2. **Prat et al. (2010)**  
   *Phenotypic and molecular characterization of the claudin-low intrinsic subtype of breast cancer.*
3. **Black & McGranahan (2021)**  
   *Genetic and non-genetic clonal diversity in cancer evolution.*
4. **Coutant et al. (2023)**  
   *Spatial transcriptomics reveal pitfalls and opportunities for the detection of rare high-plasticity breast cancer subtypes.*  
   [DOI](https://doi.org/10.1016/j.labinv.2023.100258)

---

## Organisation du Dépôt - Arborescence du projet

```
📁 **Jordan**
├── 📁 Data
│   ├── Annotations
│   │   ├── Anciennes annotations # Annotation détaillé des spots Visium
│   │   ├── Marqueurs # Marqueurs moléculaires utilisés pour l'annotation
│   │   └── Nouvelles annotations # Annotation généraliste des spots Visium
│   ├── CNA # Objets R utiles pour les scripts d'InferCNV
│   ├── Cytoband 
│   │   └── cytoBand.txt # Fichier .txt des cytobandes chez l'Homme pour InferCNV
│   ├── GTF
│   │   ├── InferCNVplus # Fichier des positions des gènes humains (adapté à InferCNV)
│   │   └── Original # Fichiers d'annotations du génome humain (initiale) pour InferCNV
│   ├── Seurat_object # Objets R de Seurat utiles lors de l'analyse transcripto spatiale
│   └── Visium
│       ├── CellType_matrix # Matrice d'expression Visium par type cellulaire (SCENIC)
│       ├── MpBC1... # Matrice d'expression brute Visium de 1 à 16 
│       └── Samples_Infos # Informations de qualité séquençage des matrice d'expression brutes Visium
├── 📁 Logiciels # Repertoire des logiciels et outils utilisés
│   ├── Case_Viewer
│   ├── GitHub
│   ├── Gnomic
│   ├── InkScape
│   ├── Loupe_Browser
│   ├── RStudio
│   └── SCENIC
│       ├── aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif # Image Singularity pour utilisé SCENIC facilement
│       ├── allTFs_hg38.txt # Fichier .txt de tous les facteurs de transcription (TF) humain hg38
│       ├── hg38_*.feather # Base de données pour la recherche de motif des TF dans les promoteurs gène
│       ├── motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl # Fichier repertoriant les motif de chaque TF
│       └── pyscenic.yml # Environnement COnda pour faire fonctionner la pipeline pyscenic
├── 📁 Results
│   ├── CNA
│   │   └── InferCNVPlus # Ensemble des figures généré avec l'analyse InferCNVPlus
│   ├── SCENIC
│   │   ├── Outputs_per_CellType
│   │   │   └── Annotations_old
│   │   │       ├── Post_SCENIC_outputs
│   │   │       │   ├── Top_regulons # Liste des top regulons avec binarisation manuelle (test, fiabilité incertaine)
│   │   │       │   └── Top_regulon_binarized # Listes des top regulons binarisé selon SCENIC (fiable, automatique)
│   │   │       │       └── binarized_outputs
│   │   │       │           ├── Epithelial_tumor_cells # Ensemble des régulons pour chaque type cellulaire de l'annotation
│   │   │       │           │   ├── Epithelial_tumor_cells_all_regulons.csv # Liste entère des régulons classés
│   │   │       │           │   ├── Epithelial_tumor_cells_binarized.csv # Matrice AUCell binarisé
│   │   │       │           │   ├── Epithelial_tumor_cells_top10_regulons.csv # Liste Top10 des régulons classés
│   │   │       │           │   └── thresholds_per_regulon.csv # Seuils de binarisation utilisés par SCENIC
│   │   │       │           ...
│   │   │       └── SCENIC_pipeline_outputs
│   │   │           ├── Epithelial_tumor_cells # Ensemble des régulons pour chaque type cellulaire de l'annotation
│   │   │           │   ├── AUCell # Matrice AUCell 
│   │   │           │   ├── Expr_matrix_adj # Matrice d'adjacence (TF/target_genes)
│   │   │           │   ├── Raw_matrix # Fichiers .loom des matrices d'expression par type cellulaire filtrées et non filtrées
│   │   │           │   └── Regulons # Matrice des régulons potentiels trouvés par SCENIC
│   │   │           ...
│   │   └── Outputs_per_MpBC # Résultats de l'analyse SCENIC initiale portant sur les patients individuellement (plutot que le type cellulaire)
│   └── Seurat # Ensemble des figures généré avec l'analyse Seurat pour les données Visium
└── 📁 Scripts
    ├── CNA_Project
    │   ├── InferCNVplus 
    │   │   ├── IronHeart_v4 # Script principal et complet pour faire l'analyse des CNA par cytobande avec InferCNVPlus (peut faire appel à des données générés avec le script IronHeart_v3)
    │   │   └── IronHeart_v3 # Script secondaire pour l'analyse InferCNVPlus
    │   └── XClone # Tentatives de scripts pour faire fonctionner XClone
    ├── SCENIC
    │   ├── Job_pyscenic
    │   │   ├── Per_celltype # Sripts et job.sh pour l'analyse SCENIC des matrice d'expression par type cellulaire
    │   │   │   └── Script_convert_MpBC_matrix_to_CellTypes_matrix # Script permettant de généré les matrices d'expression par types cellulaires à partir des matrices par MpBC
    │   │   └── Per_MpBC # Sripts et job.sh pour l'analyse SCENIC des matrice d'expression par patients
    │   ├── Old_scripts # Vieux scripts pour run SCENIC
    │   ├── Post_analysis_Scenic
    │   │   ├── Most_important_reg_CellType_SCENIC_binarization.py # Binarise automatique SCENIC et génère les Top régulons à partir de la matrice AUCell
    │   │   ├── Most_important_reg_CellType.py # Binarise manuelle et génère les Top régulons à partir de la matrice AUCell
    │   │   └── RSS_per_MpBC.py # Génère les RSS plot et calcul les régulons les plus spécifiques à chaque type cellulaire (ne marche qu'avec l'analyse SCENIC par patient)
    │   └── Tutorial_pyscenic
    └── Spatial_transcripto_Project
	     ├── Actual # Pipeline d'analyse Seurat des matrice d'expression Visium
	     ├── Old # Vieux script (peu fiables)
	     └── Tests # Scripts d'essais pour faire fonctionner Seurat et analyses complémentaires
```

## Contribution et Contact

Pour toute question ou suggestion, n’hésitez pas à contacter Pierre Martinez à l’adresse pierre.martinez@lyon.unicancer.fr ou Jordan Dutel : jordan.dutel@etu.univ-lyon1.fr
