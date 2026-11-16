"""
Leiden clustering (graph-based community detection, Traag et al. 2019).

Builds a k-nearest-neighbor graph on the shared PCA representation and runs
Scanpy's Leiden algorithm on it. Graph construction is included in the
runtime, because Leiden cannot run without it.
Key parameter: resolution (higher = more, smaller clusters).
"""

import inspect
import time

import pandas as pd
import scanpy as sc
from anndata import AnnData

NAME = "leiden"
PARAM_NAME = "resolution"


def _leiden_kwargs():
    """Pick Leiden options that work with the installed Scanpy version.

    Newer Scanpy (>= 1.10) has a `flavor` option and can use the fast igraph
    implementation. Older Scanpy does not know `flavor` and always uses the
    leidenalg package, so the option must not be passed at all.
    """
    if "flavor" not in inspect.signature(sc.tl.leiden).parameters:
        return {}  # older Scanpy: default leidenalg implementation
    try:
        import igraph  # noqa: F401
        return dict(flavor="igraph", n_iterations=2, directed=False)
    except ImportError:
        return dict(flavor="leidenalg")


def run(data, resolution, seed=0, n_neighbors=15, **_):
    """Run Leiden with the given resolution.

    Returns dict with "labels" (array of str) and "runtime_sec".
    """
    # small AnnData that only holds the PCA matrix, so the main object is untouched
    tmp = AnnData(obs=pd.DataFrame(index=data.adata.obs_names.copy()),
                  obsm={"X_pca": data.X_pca.copy()})

    start = time.perf_counter()
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(tmp, resolution=float(resolution), random_state=seed,
                 key_added="cluster", **_leiden_kwargs())
    runtime = time.perf_counter() - start

    return {"labels": tmp.obs["cluster"].astype(str).to_numpy(), "runtime_sec": runtime}
