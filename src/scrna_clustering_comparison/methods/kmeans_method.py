"""
K-means clustering (classical baseline).

Runs scikit-learn K-means on the shared PCA representation.
Key parameter: k (number of clusters).
"""

import time

from sklearn.cluster import KMeans

NAME = "kmeans"
PARAM_NAME = "k"


def run(data, k, seed=0, n_init=10, **_):
    """Cluster the cells into k groups.

    Parameters
    ----------
    data   : Dataset from utils.data_loader (uses data.X_pca)
    k      : number of clusters
    seed   : random seed
    n_init : number of K-means restarts (the best one is kept)

    Returns
    -------
    dict with "labels" (array of str, one per cell) and "runtime_sec"
    """
    X = data.X_pca
    start = time.perf_counter()
    model = KMeans(n_clusters=int(k), n_init=n_init, random_state=seed)
    labels = model.fit_predict(X)
    runtime = time.perf_counter() - start
    return {"labels": labels.astype(str), "runtime_sec": runtime}
