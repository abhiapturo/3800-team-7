"""
Evaluation metrics.

External (need reference labels):
    ARI  - Adjusted Rand Index
    NMI  - Normalized Mutual Information
Internal (no labels needed):
    Silhouette score, computed on the shared PCA representation
Stability:
    mean pairwise ARI between runs with different random seeds
"""

from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             silhouette_score)


def evaluate(pred, true, X):
    """Compute all metrics for one clustering result."""
    pred = np.asarray(pred)
    n_clusters = len(np.unique(pred))
    out = {"n_clusters": n_clusters}

    if true is not None:
        out["ari"] = adjusted_rand_score(true, pred)
        out["nmi"] = normalized_mutual_info_score(true, pred)
    else:
        out["ari"] = np.nan
        out["nmi"] = np.nan

    if 2 <= n_clusters <= len(pred) - 1:
        out["silhouette"] = silhouette_score(X, pred)
    else:
        out["silhouette"] = np.nan
    return out


def stability(label_sets):
    """Mean pairwise ARI between several clusterings of the same cells (1 = identical)."""
    if len(label_sets) < 2:
        return np.nan
    scores = [adjusted_rand_score(a, b) for a, b in combinations(label_sets, 2)]
    return float(np.mean(scores))


def summarize(metrics: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    """Average the metrics over seeds for each (method, parameter value).

    Also adds the stability across seeds for each setting.
    """
    metric_cols = ["n_clusters", "ari", "nmi", "silhouette", "runtime_sec"]
    group_cols = ["dataset", "method", "param_name", "param_value"]

    summary = (metrics.groupby(group_cols)[metric_cols]
               .agg(["mean", "std"]))
    summary.columns = [f"{m}_{s}" for m, s in summary.columns]
    summary = summary.reset_index()

    stab = []
    for _, row in summary.iterrows():
        runs = metrics[(metrics.method == row.method) &
                       (metrics.param_value == row.param_value)]["run_id"]
        stab.append(stability([assignments[r].to_numpy() for r in runs]))
    summary["stability_ari"] = stab
    summary["n_seeds"] = (metrics.groupby(group_cols).size().to_numpy())
    return summary
