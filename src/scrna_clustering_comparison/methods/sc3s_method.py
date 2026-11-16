"""
SC3s: single-cell consensus clustering with speed (Quah & Hemberg 2022).

Uses the official package (https://github.com/hemberg-lab/sc3s) on the shared
PCA representation. SC3s runs streaming K-means several times and combines
the runs into a consensus clustering.
Key parameter: k (number of clusters).
"""

import logging
import time

import pandas as pd
from anndata import AnnData

NAME = "sc3s"
PARAM_NAME = "k"


def run(data, k, seed=0, n_runs=5, d_range=None, **_):
    """Run SC3s consensus clustering.

    Parameters
    ----------
    data    : Dataset from utils.data_loader (uses data.X_pca)
    k       : number of clusters
    seed    : random seed
    n_runs  : number of K-means realisations combined in the consensus
    d_range : list of numbers of PCs to use (None = SC3s default: min(25, n_pcs))

    Returns dict with "labels" (array of str) and "runtime_sec".
    """
    import sc3s  # imported here so the other methods work even if sc3s is missing

    k = int(k)
    tmp = AnnData(obs=pd.DataFrame(index=data.adata.obs_names.copy()),
                  obsm={"X_pca": data.X_pca.copy()})

    # SC3s prints many INFO/WARNING messages; hide them during the run
    root = logging.getLogger()
    old_level = root.level
    root.setLevel(logging.ERROR)
    try:
        start = time.perf_counter()
        sc3s.tl.consensus(tmp, n_clusters=[k], d_range=d_range,
                          n_runs=n_runs, random_state=seed)
        runtime = time.perf_counter() - start
    finally:
        root.setLevel(old_level)

    labels = tmp.obs[f"sc3s_{k}"].astype(str).to_numpy()
    return {"labels": labels, "runtime_sec": runtime}
