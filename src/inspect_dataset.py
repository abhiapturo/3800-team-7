"""This script inspects our dataset and returns:
    - how many cells in our dataset
    - how many represented genes
    - how many reference cell types exist
    - Any missing labels?
    - Distribution of cells among reference cell types
"""
from pathlib import Path
import anndata as ad
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent

with open(ROOT / "config.yaml", "r") as file:
    """Load the shared configuration"""
    config = yaml.safe_load(file)

adata = ad.read_h5ad(ROOT / config["data_path"])

labels = adata.obs[config["label_column"]]

print("Number of cells:", adata.n_obs)
print("Number of genes:", adata.n_vars)
print("Number of reference classes:", labels.nunique())
print("Missing labels:", labels.isna().sum())

print("\nReference label distribution:")
print(labels.value_counts(dropna=False))
