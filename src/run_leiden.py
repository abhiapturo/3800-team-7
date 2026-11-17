from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import yaml


# Project root variable
ROOT = Path(__file__).resolve().parent.parent

with open(ROOT / "config.yaml", "r") as file:
    # Load shared experimental settings
    config = yaml.safe_load(file)

# Load preprocessed dataset from ROOT
adata = ad.read_h5ad(ROOT / config["data_path"])


# Verify PCA features exist
if "X_pca" not in adata.obsm:
    raise ValueError("PCA representation is missing")

# Check that enough PCA components exist
if adata.obsm["X_pca"].shape[1] < config["n_pcs"]:
    raise ValueError("Not enough PCA components")

# Verify shared neighbor graph
if "neighbors" not in adata.uns:
    raise ValueError("NOT FOUND!!! Run preprocess.py")

graph_key = adata.uns["neighbors"].get(
    "connectivities_key", "connectivities"
)

if graph_key not in adata.obsp:
    raise ValueError("Neighbor connectivity matrix missing")


# Number of reference cell types
label_column = config["label_column"]

if label_column not in adata.obs.columns:
    raise ValueError("Reference labels are missing")

reference_labels = adata.obs[label_column]

if reference_labels.isna().any():
    raise ValueError("Some reference labels are missing")

# The number of unique reference cell types in our data set
reference_k = reference_labels.nunique()

# Use the first seed for first experiments
seed = config["seeds"][0]

# Include resolution 1.0 as the baseline
resolutions = sorted(
    set(config["leiden_resolutions"] + [1.0])
)

# Temporary storage for runs
all_cluster_results = []
all_runtime_results = []

for resolution in resolutions:

    print(f"\nRunning Leiden at resolution {resolution}")

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
    labels = adata.obs["leiden"].astype(str).to_numpy()
    n_clusters = len(np.unique(labels))

    run_id = f"lr[{resolution}]seed[{seed}]"

    print(f"Clusters found: {n_clusters}")
    print(f"Runtime: {runtime:.4f} seconds")

    # Store clustering results for this run
    cluster_results = pd.DataFrame({
        "cell_id": adata.obs_names,
        "cluster": labels,
        "method": "leiden",
        "seed": seed,
        "parameter": "resolution",
        "parameter_value": resolution,
        "run_id": run_id
    })

    all_cluster_results.append(cluster_results)

    # Store runtime information for this run
    all_runtime_results.append({
        "run_id": run_id,
        "method": "leiden",
        "parameter":"resolution",
        "parameter_value": resolution,
        "seed": seed,
        "runtime_seconds": round(runtime,4),
        "n_clusters": n_clusters
    })


runtime_results = pd.DataFrame(all_runtime_results)


baseline = runtime_results.loc[
    np.isclose(runtime_results["parameter_value"], 1.0)
].iloc[0]

# Matched closest cluster count to reference_k
differences = abs(
    runtime_results["n_clusters"] - reference_k
)

matched = runtime_results.iloc[differences.argmin()]

selected_run_ids = {
    baseline["run_id"],
    matched["run_id"]
}

# Combine only selected cluster assignments
cluster_results = pd.concat(
    [
        df for df in all_cluster_results
        if df["run_id"].iloc[0] in selected_run_ids
    ],
    ignore_index=True
)

selected_runtime = runtime_results[
    runtime_results["run_id"].isin(selected_run_ids)
].copy()


output_dir = ROOT / "results"
output_dir.mkdir(parents=True, exist_ok=True)

cluster_results.to_csv(
    output_dir / "clusters_leiden.csv",
    index=False
)

print("\nSaved results/clusters_leiden.csv")

selected_runtime.to_csv(
    output_dir / "runtime_leiden.csv",
    index=False
)

print("Saved results/runtime_leiden.csv")

# Save the entire resolution search summary
runtime_results.to_csv(
    output_dir / "leiden_resolution_search.csv",
    index=False
)


# Create Umap plots

figures_dir = ROOT / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)

# Compute Umap only if it is not already available
if "X_umap" not in adata.obsm:
    sc.tl.umap(
        adata,
        random_state=config["random_state"]
    )

for run in [baseline, matched]:

    run_id = run["run_id"]

    # Retrieve this run's labels in dataset order
    run_labels = cluster_results.loc[
        cluster_results["run_id"] == run_id
    ].set_index("cell_id").loc[
        adata.obs_names, "cluster"
    ]

    adata.obs["plot_cluster"] = pd.Categorical(
        run_labels.to_numpy()
    )

    sc.pl.umap(
        adata,
        color="plot_cluster",
        title=f"Leiden resolution {run['parameter_value']}",
        show=False
    )

    plt.savefig(
        figures_dir / f"{run_id}_umap.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close("all")


print("\nReference cell types:", reference_k)

print("\nBaseline Leiden:")
print(baseline.to_string())

print("\nMatched Leiden:")
print(matched.to_string())
