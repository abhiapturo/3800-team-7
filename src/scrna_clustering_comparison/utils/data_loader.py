"""
Loading and preparing single-cell datasets.

Every dataset ends up in the same form (a `Dataset` object) so that all three
clustering methods and the evaluation code can work with any dataset:

    data.adata    AnnData object with PCA (obsm["X_pca"]) and UMAP (obsm["X_umap"])
    data.X_pca    the PCA matrix used by ALL methods (cells x n_pcs)
    data.labels   reference cell-type labels (or None if not available)

Supported inputs
----------------
1. Built-in Scanpy datasets (see BUILTIN_DATASETS), e.g. "pbmc68k_reduced"
2. Your own file:
      - .h5ad file (AnnData)
      - .csv file (rows = cells, columns = genes)
      - a 10x Genomics folder with matrix.mtx / genes.tsv / barcodes.tsv
   Labels can be a column in adata.obs (--label-key) or a separate CSV file
   (--labels-file, first column = cell ID, second column = label).
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

log = logging.getLogger(__name__)

# name -> how to load it, which obs column holds the labels,
#         and whether the data is already normalized/scaled
BUILTIN_DATASETS = {
    # 700 PBMCs x 765 genes, already processed (our main dataset)
    "pbmc68k_reduced": dict(loader=sc.datasets.pbmc68k_reduced,
                            label_key="bulk_labels", processed=True),
    # 2,638 PBMCs, already processed (labels come from Louvain clustering)
    "pbmc3k_processed": dict(loader=sc.datasets.pbmc3k_processed,
                             label_key="louvain", processed=True),
    # 2,730 mouse bone-marrow cells, raw counts (downloads on first use)
    "paul15": dict(loader=sc.datasets.paul15,
                   label_key="paul15_clusters", processed=False),
}


@dataclass
class Dataset:
    name: str
    adata: AnnData
    label_key: Optional[str]
    n_pcs: int

    @property
    def X_pca(self) -> np.ndarray:
        """PCA representation shared by all clustering methods."""
        return self.adata.obsm["X_pca"][:, : self.n_pcs]

    @property
    def labels(self) -> Optional[np.ndarray]:
        if self.label_key is None:
            return None
        return self.adata.obs[self.label_key].astype(str).to_numpy()

    @property
    def n_cells(self) -> int:
        return self.adata.n_obs


def _read_file(path: Path) -> AnnData:
    if path.is_dir():
        return sc.read_10x_mtx(path)
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return sc.read_h5ad(path)
    if suffix == ".csv":
        return sc.read_csv(path, first_column_names=True)
    raise ValueError(f"Unsupported file type: {path}. Use .h5ad, .csv or a 10x folder.")


def load_raw(source: str,
             label_key: Optional[str] = None,
             labels_file: Optional[str] = None,
             processed: Optional[bool] = None):
    """Load a dataset without any processing.

    Returns (adata, dataset_name, label_key, processed).
    """
    if source in BUILTIN_DATASETS:
        info = BUILTIN_DATASETS[source]
        adata = info["loader"]()
        name = source
        label_key = label_key or info["label_key"]
        processed = info["processed"] if processed is None else processed
    else:
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(
                f"'{source}' is not a built-in dataset ({list(BUILTIN_DATASETS)}) "
                f"and no file/folder with that name exists.")
        adata = _read_file(path)
        name = path.stem
        processed = bool(processed)  # files are treated as raw counts unless told otherwise

    adata.var_names_make_unique()

    if labels_file is not None:
        lab = pd.read_csv(labels_file, index_col=0).iloc[:, 0]
        adata.obs["label"] = lab.reindex(adata.obs_names).astype(str)
        label_key = "label"

    if label_key is not None and label_key not in adata.obs:
        raise KeyError(f"Label column '{label_key}' not found in adata.obs. "
                       f"Available columns: {list(adata.obs.columns)}")

    return adata, name, label_key, processed


def preprocess(adata: AnnData, processed: bool, n_top_genes: int = 2000,
               min_genes: int = 200, min_cells: int = 3) -> AnnData:
    """Standard Scanpy preprocessing, skipped if the data is already processed."""
    adata = adata.copy()
    if processed:
        log.info("Data marked as already processed: skipping normalization.")
        return adata

    log.info("Preprocessing raw counts (filter, normalize, log, HVG, scale).")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars),
                                subset=True)
    sc.pp.scale(adata, max_value=10)
    return adata


def load_dataset(source: str = "pbmc68k_reduced",
                 label_key: Optional[str] = None,
                 labels_file: Optional[str] = None,
                 processed: Optional[bool] = None,
                 n_pcs: int = 30,
                 n_neighbors: int = 15,
                 n_top_genes: int = 2000,
                 seed: int = 0) -> Dataset:
    """Load, preprocess, and compute PCA + UMAP. Returns a ready-to-use Dataset."""
    adata, name, label_key, processed = load_raw(source, label_key, labels_file, processed)
    log.info("Loaded '%s': %d cells x %d genes", name, adata.n_obs, adata.n_vars)

    adata = preprocess(adata, processed, n_top_genes=n_top_genes)

    # drop cells without a reference label (they cannot be evaluated)
    if label_key is not None:
        keep = adata.obs[label_key].notna() & (adata.obs[label_key].astype(str) != "nan")
        if (~keep).sum() > 0:
            log.info("Removing %d cells without a label", int((~keep).sum()))
            adata = adata[keep].copy()

    n_pcs = int(min(n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=seed)

    # neighbor graph + UMAP are only used for visualization here
    # (Leiden builds its own graph inside leiden_method.py so its runtime is fair)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(adata, random_state=seed)

    if label_key is not None:
        n_types = adata.obs[label_key].nunique()
        log.info("Reference labels: '%s' (%d cell types)", label_key, n_types)
    else:
        log.info("No reference labels: ARI/NMI will be skipped.")

    return Dataset(name=name, adata=adata, label_key=label_key, n_pcs=n_pcs)
