'''
Build shared input for all clustering methods

Load PBMC6K reduced, apply PCA, the neighbour graph and UMAP
with settings in config.yaml and saves result to data/pbmc_pca.h5ad
'''
from pathlib import Path

import scanpy as sc
import yaml
import importlib


config = yaml.safe_load(open("config.yaml"))
seed = config["random_state"]

# load the data (NOTE: The reduced data is already normalized and scaled)
adata = sc.datasets.pbmc68k_reduced()
label = config["label_column"]
print(f"Loaded {adata.n_obs} cells X {adata.n_vars} genes")
print(f"{adata.obs[label].nunique()} reference cell types in `{label}':")
print(adata.obs[label].value_counts().to_string())

sc.pp.pca(adata, n_comps = config["n_pcs"], random_state = seed)
sc.pp.neighbors(adata, n_neighbors = config["n_neighbors"], n_pcs =config["n_pcs"], random_state= seed)
sc.tl.umap(adata, random_state=seed)

adata.uns["preprocessing_setting"] = {
    "n_pcs": config["n_pcs"],
    "n_neighbors": config["n_neighbors"],
    "random_state":seed,
    "scanpy_version": importlib.metadata.version('scanpy'),
}

out = Path(config["data_path"])
out.parent.mkdir(parents = True, exist_ok =  True)
adata.write_h5ad(out)
print(f"Saved {out} (X_pca shape: {adata.obsm['X_pca'].shape})")