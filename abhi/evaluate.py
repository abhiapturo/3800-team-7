""" This python script is a shared clustering evaluation between
K-means, Leiden, and SC3s. Measurements are in:
Adjusted Rand Index
Normalized Mutual Information
Silhouette Score
Number of clusters
Runtime
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config.yaml"
RESULTS_DIR = ROOT / "results"


def load_config():
    """Load shared experimental settings."""
    with open(CONFIG_PATH, "r") as file:
        return yaml.safe_load(file)


def load_dataset(config):
    """Load PCA features and reference cell labels."""

    data_path = ROOT / config["data_path"]
    adata = sc.read_h5ad(data_path)

    label_column = config["label_column"]

    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Missing reference column: {label_column}"
        )

    if "X_pca" not in adata.obsm:
        raise ValueError("Dataset has no X_pca representation")

    n_pcs = config["n_pcs"]

    if adata.obsm["X_pca"].shape[1] < n_pcs:
        raise ValueError("Insufficient PCA components")

    X = np.asarray(
        adata.obsm["X_pca"][:, :n_pcs]
    )

    true_labels = adata.obs[label_column]

    if true_labels.isna().any():
        raise ValueError(
            "Reference labels contain missing values"
        )

    return adata.obs_names, X, true_labels.to_numpy()


def evaluate_clustering(
    X,
    true_labels,
    pred_labels,
    method,
    runtime_seconds=np.nan,
    seed=None,
    parameter=None,
    parameter_value=None,
):
    """Evaluate one clustering run."""

    X = np.asarray(X)
    true_labels = np.asarray(true_labels)
    pred_labels = np.asarray(pred_labels)

    if X.ndim != 2:
        raise ValueError("X must be a 2D matrix")

    if (
        len(X) != len(true_labels)
        or len(X) != len(pred_labels)
    ):
        raise ValueError("Input lengths do not match")

    if not np.isfinite(X).all():
        raise ValueError("PCA features contain invalid values")

    if pd.isna(pred_labels).any():
        raise ValueError("Predicted labels contain missing values")

    # External clustering agreement
    ari = adjusted_rand_score(
        true_labels, pred_labels
    )

    nmi = normalized_mutual_info_score(
        true_labels, pred_labels
    )

    # Number of predicted clusters
    n_clusters = len(np.unique(pred_labels))

    # Silhouette requires 2 <= clusters < samples
    if 2 <= n_clusters < len(X):
        silhouette = silhouette_score(
            X,
            pred_labels,
            sample_size=min(2000, len(X)),
            random_state=0,
        )
    else:
        silhouette = np.nan

    return {
        "method": method,
        "seed": seed,
        "parameter": parameter,
        "parameter_value": parameter_value,
        "ari": ari,
        "nmi": nmi,
        "silhouette": silhouette,
        "runtime_seconds": runtime_seconds,
        "n_clusters": n_clusters,
    }


def evaluate_file(csv_path, cell_ids, X, true_labels):
    """
    Evaluate a CSV containing clustering assignments.

    Required columns:
        cell_id, cluster

    Optional columns:
        method, seed, parameter,
        parameter_value, runtime_seconds

    Multiple runs can be stored in one CSV.
    """

    df = pd.read_csv(csv_path, dtype={"cell_id": str})

    required = {"cell_id", "cluster"}

    if not required.issubset(df.columns):
        raise ValueError(
            f"{csv_path.name} must contain "
            "'cell_id' and 'cluster' columns"
        )

    # Default method from filename
    default_method = (
        csv_path.stem.replace("clusters_", "")
    )

    # Assign metadata defaults
    if "method" not in df.columns:
        df["method"] = default_method

    for column in [
        "seed",
        "parameter",
        "parameter_value",
        "runtime_seconds",
    ]:
        if column not in df.columns:
            df[column] = np.nan

    group_columns = [
        "method",
        "seed",
        "parameter",
        "parameter_value",
    ]

    results = []

    for keys, group in df.groupby(
        group_columns, dropna=False
    ):

        method, seed, parameter, parameter_value = keys

        # Ensure exactly one prediction per cell
        if group["cell_id"].duplicated().any():
            raise ValueError(
                f"Duplicate cell IDs in {csv_path.name}"
            )

        indexed = group.set_index("cell_id")

        # Check cell IDs match the dataset
        expected = pd.Index(cell_ids.astype(str))

        if (
            len(indexed) != len(expected)
            or set(indexed.index) != set(expected)
        ):
            raise ValueError(
                f"Cell IDs do not match dataset "
                f"in {csv_path.name}"
            )

        # Align predicted labels to dataset order
        pred_labels = indexed.loc[
            expected, "cluster"
        ].to_numpy()

        runtimes = group[
            "runtime_seconds"
        ].dropna().unique()

        if len(runtimes) > 1:
            raise ValueError(
                "Multiple runtimes found for one run"
            )

        runtime = (
            float(runtimes[0])
            if len(runtimes) == 1
            else np.nan
        )

        result = evaluate_clustering(
            X=X,
            true_labels=true_labels,
            pred_labels=pred_labels,
            method=method,
            runtime_seconds=runtime,
            seed=seed,
            parameter=parameter,
            parameter_value=parameter_value,
        )

        results.append(result)

    return results


def main():

    parser = argparse.ArgumentParser(
        description="Evaluate scRNA-seq clustering"
    )

    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Optional clustering CSV file",
    )

    args = parser.parse_args()

    config = load_config()

    cell_ids, X, true_labels = load_dataset(config)

    RESULTS_DIR.mkdir(
        parents=True, exist_ok=True
    )

    if args.input:
        input_files = [Path(args.input)]
    else:
        input_files = sorted(
            RESULTS_DIR.glob("clusters_*.csv")
        )

    if not input_files:
        raise FileNotFoundError(
            "No clustering CSV files found"
        )

    all_results = []

    for csv_path in input_files:

        print(f"Evaluating {csv_path.name}...")

        results = evaluate_file(
            csv_path,
            cell_ids,
            X,
            true_labels,
        )

        all_results.extend(results)

    metrics = pd.DataFrame(all_results)

    output_path = RESULTS_DIR / "metrics.csv"

    metrics.to_csv(output_path, index=False)

    print("\nEvaluation complete!\n")
    print(metrics.to_string(index=False))
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()
