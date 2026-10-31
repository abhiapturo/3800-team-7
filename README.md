# Clustering Methods for Single Cell RNA-Sequencing Data

A reproducible comparison of three clustering methods (K-means, Leiden, and SC3s) on the same single-cell RNA seq dataset. We measure how well each method recovers known cell types, how well separated its clusters are, how many clusters it finds, how long it takes, and how much its results change with key parameters and random seeds.


## Contents
* Background
* Methods compared
* Dataset
* Repository Structure
* Setup
* How to run
* Shared file formats
* Evaluation 
* Team and responsibilities
* Timeline
* References

### Background
scRNA-seq measures gene expression in individual cells. Clustering groups cells with similar expression profiles, and these groups are then interpreted as cell types. Different clustering methods make different assumptions, so they can produce different groupings, with different stability and computational cost.

Aim
1. Compare agreement between each method's clusters and the reference cell-type labels.
2. Compare cluster separation, runtime and the number of clusters each method finds
3. Examine how results change the key parameters (number of cluster K, Leiden resolution) and random seeds.

### Methods compared

### Dataset
PBMC68k reduced, bundled with Scanpy:

```Python
import scanpy as sc
adata = sc.datasets.pbmc68k_reduced()
```

### Repository Structure
scrna-clustering-comparison/
├── README.md
├── environment.yml          # Conda environment (all dependencies, pinned)
├── config.yaml              # Shared settings: PCs, neighbours, seeds, k values, resolutions
├── .gitignore
│
├── src/
│   ├── preprocess.py        # Load data, PCA, neighbour graph, UMAP → data/pbmc_pca.h5ad
│   ├── run_leiden.py        # Leiden runs and resolution sweep
│   ├── run_kmeans.py        # K-means runs and k sweep
│   ├── run_sc3s.py          # SC3s runs and k sweep
│   ├── evaluate.py          # ARI, NMI, silhouette, runtime → results/metrics.csv
│   └── make_figures.py      # All final figures → figures/
│
├── data/                    # Generated; not committed
│   └── pbmc_pca.h5ad
│
├── results/
│   ├── clusters_leiden.csv
│   ├── clusters_kmeans.csv
│   ├── clusters_sc3s.csv
│   └── metrics.csv
│
├── figures/                 # Final figures used in the report and slides
├── notebooks/               # Exploration only; final results come from src/
└── report/                  # Report drafts and final PDF


### Setup
Requirements: Conda and git.

```bash
# 1. Clone the repository
git clone https://github.com/[username]/scrna-clustering-comparison.git
cd scrna-clustering-comparison

# 2. Create and activate the environment
conda env create -f environment.yml
conda activate scrna-clustering
```

Main dependencies
- Python 3.10+
- scanpy, anndata
- scikit-learn
- leidenalg, python-igraph
- sc3s
- pandas, numpy, matplotlib, seaborn
- pyyaml


### How to run
Run the steps in order from the repository root. All settings are read from config.yaml

```bash
# 1. Preprocess: load data, PCA, neighbour graph, UMAP
python src/preprocess.py

# 2. Run each clustering method (default runs, parameter sweeps and seed repeats)
python src/run_leiden.py
python src/run_kmeans.py
python src/run_sc3s.py

# 3. Score every run
python src/evaluate.py

# 4. Make all figures
python src/make_figures.py
```
Each step can also be run on its own once the previous steps output exists.

### Shared settings
All members use these values. Change them only after the group agrees.

```yaml
n_pcs: 30
n_neighbors: 15
seeds: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
k_values: [4, 6, 8, 10, 12, 14]
leiden_resolutions: [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
label_column: bulk_labels
```

### Shared file formats
Every method saves result in the same format, so one script can evaluate all of them.

** Cluster assignments **: `results/clusters_<method>.csv` one row per cell per run:

| Column            | Description                         |
| :---------------- | ----------------------------------: |
| cell_id           | Cell barcode from adata.obs_names   |
| Cluster           | Cluster assigned by the method      |
| method            | Kmeans, Leiden or SC3S              |
| parameter         | Value of k or Leiden resolution     |
| seed              | Random seed used                    |


** Metrics **: `results/metrics.csv` one row per cell per run:

| Column            | Description                                        |
| :---------------- | -------------------------------------------------: |
| method            | kmeans, leiden or sc3s                             |
| parameter         | value of k or leiden resolution                    |
| seed              | Random seed used                                   |
| n_clusters        | Number of clusters found                           |
| ARI               | Adjusted Rand Index vs reference labels            |
| NMI               | Normalized mutual informations vs reference labels |
| Silhouette        | Silhouette score in PCA space                      |
| runtime           | runtime in seconds                                 |



