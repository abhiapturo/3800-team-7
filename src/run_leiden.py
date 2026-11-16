from pathlib import Path
import time

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


# Project root variable
ROOT = Path(__file__).resolve().parent.parent

with open(ROOT / "config.yaml", "r") as file:
    """Load shared experimental settings"""
    config = yaml.safe_load(file)

# Loads preprocessed dataset from ROOT
adata = sc.read_h5ad(ROOT / config["data_path"])


if "X_pca" not in adata.obsm:
    """Check if PCA features exist"""
    raise ValueError("PCA representation is missing")

if adata.obsm["X_pca"].shape[1] < config["n_pcs"]:
    """Need atleast 30 PCA components, checks for that amount"""
    raise ValueError("Not enough PCA components")

# Verify shared neighbor graph
if "neighbors" not in adata.uns:
    raise ValueError(
        "NOT FOUND!!! Run preprocess.py"
    )

# first Leiden run with initial setting (resolution = 1.0) ONLY
resolution = 1.0
seed = config["seeds"][0]

start = time.perf_counter()

sc.tl.leiden(
    adata,
    resolution=resolution,
    random_state=seed,
    key_added="leiden",
    flavor="igraph",
    n_iterations=2,
    directed=False
)

runtime = time.perf_counter() - start

# Extract cluster assignments
labels = adata.obs["leiden"].to_numpy()
n_clusters = len(np.unique(labels))

print("Leiden clustering complete")
print(f"Resolution: {resolution}")
print(f"Seed: {seed}")
print(f"Clusters found: {n_clusters}")
print(f"Runtime: {runtime:.4f} seconds")

# Save results
output_dir = ROOT / "results"
output_dir.mkdir(exist_ok=True)

results = pd.DataFrame({
    "cell_id": adata.obs_names,
    "cluster": labels
})

results.to_csv(
    output_dir / "clusters_leiden.csv",
    index=False
)

print("Saved results/clusters_leiden.csv")
