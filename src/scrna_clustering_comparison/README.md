# Comparative Evaluation of Clustering Methods for scRNA-seq Data

CSE 3800 / BME 4800 / CSE 5800 Bioinformatics, Fall 2026 — Group #7

Compares **K-means**, **Leiden**, and **SC3s** on the same PCA representation of a
single-cell dataset, over a grid of parameters and random seeds. Every run is
evaluated with ARI, NMI, Silhouette score, runtime, and stability across seeds,
and all results are saved so plots can be made later.

## Project structure

```
scrna_clustering_comparison/
├── main.py                 # runs everything: load data -> cluster -> evaluate -> save
├── make_plots.py           # makes all plots from saved results (no re-running)
├── methods/
│   ├── kmeans_method.py    # K-means (scikit-learn)       parameter: k
│   ├── leiden_method.py    # Leiden (Scanpy)              parameter: resolution
│   └── sc3s_method.py      # SC3s (official sc3s package)  parameter: k
├── utils/
│   ├── data_loader.py      # built-in or your own dataset -> preprocessing, PCA, UMAP
│   ├── evaluation.py       # ARI, NMI, Silhouette, stability, summary table
│   └── plotting.py         # all comparison figures
└── requirements.txt
```

All three methods use the **same** PCA matrix (`data.X_pca`), so the comparison is fair.
Reference labels are used only for evaluation, never for clustering.

## Installation

```bash
python -m pip install -r requirements.txt
python -m pip install --no-deps sc3s
```

SC3s is installed with `--no-deps` because its package metadata pins `numpy<2`;
the code runs fine with numpy 2, and this avoids downgrading numpy in your environment.
In Spyder you can run the same commands in the console with `%pip` instead of
`python -m pip`, then restart the kernel.

## Running

```bash
# default: PBMC68k reduced (700 cells, built into Scanpy), all methods, 5 seeds, then plots
python main.py --plot

# quick test
python main.py --seeds 0 --k-values 8 10 --resolutions 0.5 1.0 --plot

# only some methods
python main.py --methods kmeans leiden

# other built-in datasets (downloaded on first use)
python main.py --dataset pbmc3k_processed --plot
python main.py --dataset paul15 --plot          # raw counts, preprocessed automatically

# your own data
python main.py --dataset my_data.h5ad --label-key cell_type --raw --plot
python main.py --dataset my_data.h5ad --label-key cell_type --processed
python main.py --dataset counts.csv --labels-file labels.csv --raw
python main.py --dataset path/to/10x_folder/ --labels-file labels.csv --raw

# re-make plots later (choose the metric used to pick each method's best setting)
python make_plots.py results/pbmc68k_reduced --select-by nmi
```

`--raw` runs filtering, normalization, log transform, highly variable genes, and
scaling. `--processed` skips this (use it if your data is already normalized).
A labels file is a CSV with cell IDs in the first column and labels in the second.

Main options (see `python main.py --help`):

| Option | Default | Meaning |
|---|---|---|
| `--k-values` | 4 6 8 10 12 14 | k for K-means and SC3s (the true number of cell types is added automatically; `--no-true-k` to disable) |
| `--resolutions` | 0.2 … 2.0 | Leiden resolutions |
| `--seeds` | 0 1 2 3 4 | random seeds (used for stability) |
| `--n-pcs` | 30 | number of principal components |
| `--n-neighbors` | 15 | neighbors for the Leiden graph and UMAP |
| `--sc3s-runs` | 5 | SC3s realisations per consensus |
| `--save-h5ad` | off | also save an .h5ad with every clustering in `.obs` |

## Output (`results/<dataset>/`)

| File | Content |
|---|---|
| `metrics.csv` | one row per run: method, parameter, seed, n_clusters, ARI, NMI, Silhouette, runtime |
| `summary.csv` | mean and std over seeds for each method + parameter, plus stability (mean pairwise ARI between seeds) |
| `assignments.csv` | cluster label of every cell for every run (+ reference label) |
| `embedding.csv` | UMAP coordinates and reference labels |
| `best_settings.csv` | best parameter of each method (written by the plotting step) |
| `config.json` | all settings used, for reproducibility |
| `figures/` | all plots |

Figures: parameter sensitivity, metrics vs number of clusters, best-setting comparison,
runtime, variation across seeds, UMAP comparison, and reference-vs-cluster contingency
heatmaps for each method.

## Adding another method

Create `methods/my_method.py` with `NAME`, `PARAM_NAME`, and
`run(data, param, seed=0, **_)` returning `{"labels": ..., "runtime_sec": ...}`,
then add it to `METHODS` in `methods/__init__.py`. It is picked up automatically
by `main.py`, the evaluation, and the plots.

## Related repositories

- SC3s (official code from the paper): https://github.com/hemberg-lab/sc3s
- SC3 (original R version): https://github.com/hemberg-lab/SC3
- Scanpy: https://github.com/scverse/scanpy
- Leiden algorithm (leidenalg): https://github.com/vtraag/leidenalg
- Clustering benchmark by Duò, Robinson & Soneson (2018): https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison
- Scanpy tutorials: https://github.com/scverse/scanpy-tutorials
