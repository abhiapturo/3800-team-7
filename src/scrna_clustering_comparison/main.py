"""
Comparative evaluation of clustering methods for scRNA-seq data.
CSE 3800 / BME 4800 / CSE 5800 Final Project, Group #7.

Runs K-means, Leiden and SC3s on the same PCA representation of a dataset,
over a grid of parameters and random seeds, evaluates every run, and saves all
results so plots can be made later.

Examples
--------
# default: PBMC68k reduced (Scanpy built-in), all three methods
python main.py

# quick test run
python main.py --seeds 0 --k-values 8 10 --resolutions 0.5 1.0

# another built-in dataset
python main.py --dataset pbmc3k_processed

# your own data (raw counts in an .h5ad file, labels in obs["cell_type"])
python main.py --dataset path/to/data.h5ad --label-key cell_type --raw

# run, then make all comparison plots
python main.py --plot
"""

import warnings

# hide harmless deprecation messages from scanpy/anndata
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

from methods import METHODS
from utils.data_loader import BUILTIN_DATASETS, load_dataset
from utils.evaluation import evaluate, summarize

DEFAULT_K = [4, 6, 8, 10, 12, 14]
DEFAULT_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
DEFAULT_SEEDS = [0, 1, 2, 3, 4]

log = logging.getLogger("main")


def parse_args():
    p = argparse.ArgumentParser(description="Compare K-means, Leiden and SC3s on scRNA-seq data.")
    # data
    p.add_argument("--dataset", default="pbmc68k_reduced",
                   help=f"built-in name {list(BUILTIN_DATASETS)} or path to .h5ad/.csv/10x folder")
    p.add_argument("--label-key", default=None, help="column in adata.obs with reference labels")
    p.add_argument("--labels-file", default=None, help="CSV with cell IDs and labels (optional)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--processed", dest="processed", action="store_true", default=None,
                   help="data is already normalized/scaled (skip preprocessing)")
    g.add_argument("--raw", dest="processed", action="store_false",
                   help="data is raw counts (run normalization, HVG, scaling)")
    p.add_argument("--n-pcs", type=int, default=30, help="number of principal components")
    p.add_argument("--n-neighbors", type=int, default=15, help="neighbors for Leiden graph and UMAP")
    p.add_argument("--n-top-genes", type=int, default=2000, help="highly variable genes (raw data only)")
    # methods and parameters
    p.add_argument("--methods", nargs="+", default=list(METHODS), choices=list(METHODS))
    p.add_argument("--k-values", nargs="+", type=int, default=DEFAULT_K,
                   help="k values for K-means and SC3s")
    p.add_argument("--no-true-k", action="store_true",
                   help="do NOT add the number of reference cell types to the k values")
    p.add_argument("--resolutions", nargs="+", type=float, default=DEFAULT_RESOLUTIONS,
                   help="resolution values for Leiden")
    p.add_argument("--seeds", nargs="+", type=int, default=DEFAULT_SEEDS)
    p.add_argument("--sc3s-runs", type=int, default=5, help="SC3s realisations per consensus")
    # output
    p.add_argument("--outdir", default="results")
    p.add_argument("--save-h5ad", action="store_true", help="also save an .h5ad with all clusterings")
    p.add_argument("--plot", action="store_true", help="make all plots after running")
    return p.parse_args()


def parameter_grid(method, args, n_true_types):
    if method.PARAM_NAME == "resolution":
        return sorted(set(args.resolutions))
    ks = set(args.k_values)
    if n_true_types and not args.no_true_k:
        ks.add(n_true_types)
    return sorted(ks)


def run_all(data, args):
    """Run every method x parameter x seed. Returns (metrics_df, assignments_df)."""
    true = data.labels
    n_true = len(np.unique(true)) if true is not None else None
    rows, assignments = [], {}

    for name in args.methods:
        method = METHODS[name]
        grid = parameter_grid(method, args, n_true)
        log.info("== %s: %s = %s, seeds = %s", name, method.PARAM_NAME, grid, args.seeds)
        for value in grid:
            for seed in args.seeds:
                try:
                    res = method.run(data, value, seed=seed,
                                     n_neighbors=args.n_neighbors,
                                     n_runs=args.sc3s_runs)
                except ImportError as e:
                    log.error("Skipping %s: %s (is the package installed?)", name, e)
                    break
                scores = evaluate(res["labels"], true, data.X_pca)
                run_id = f"{name}_{method.PARAM_NAME}={value:g}_seed={seed}"
                rows.append({"run_id": run_id, "dataset": data.name, "method": name,
                             "param_name": method.PARAM_NAME, "param_value": value,
                             "seed": seed, **scores, "runtime_sec": res["runtime_sec"]})
                assignments[run_id] = res["labels"]
                log.info("%-32s clusters=%2d  ARI=%.3f  NMI=%.3f  Sil=%.3f  %.2fs",
                         run_id, scores["n_clusters"], scores["ari"], scores["nmi"],
                         scores["silhouette"], res["runtime_sec"])
            else:
                continue
            break  # only reached after an ImportError

    metrics = pd.DataFrame(rows)
    assign_df = pd.DataFrame(assignments, index=data.adata.obs_names)
    return metrics, assign_df


def save_results(data, metrics, assign_df, summary, args, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    true = data.labels

    metrics.to_csv(outdir / "metrics.csv", index=False)
    summary.to_csv(outdir / "summary.csv", index=False)

    assign_out = assign_df.copy()
    if true is not None:
        assign_out.insert(0, "true_label", true)
    assign_out.index.name = "cell"
    assign_out.to_csv(outdir / "assignments.csv")

    emb = pd.DataFrame(data.adata.obsm["X_umap"][:, :2], columns=["UMAP1", "UMAP2"],
                       index=data.adata.obs_names)
    emb["true_label"] = true if true is not None else "NA"
    emb.index.name = "cell"
    emb.to_csv(outdir / "embedding.csv")

    config = vars(args).copy()
    config.update(dataset_name=data.name, n_cells=data.n_cells,
                  n_genes=int(data.adata.n_vars), n_pcs_used=data.n_pcs,
                  label_key=data.label_key,
                  finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(outdir / "config.json", "w") as f:
        json.dump(config, f, indent=2, default=str)

    if args.save_h5ad:
        adata = data.adata.copy()
        for col in assign_df.columns:
            adata.obs[col] = pd.Categorical(assign_df[col])
        adata.write_h5ad(outdir / f"{data.name}_clustered.h5ad")

    log.info("Results saved in %s", outdir.resolve())


def print_best(summary, has_labels):
    key = "ari_mean" if has_labels else "silhouette_mean"
    best = summary.loc[summary.groupby("method")[key].idxmax()]
    cols = ["method", "param_name", "param_value", "n_clusters_mean", "ari_mean",
            "nmi_mean", "silhouette_mean", "stability_ari", "runtime_sec_mean"]
    print(f"\nBest setting per method (by {key}):")
    print(best[cols].round(3).to_string(index=False))


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")

    data = load_dataset(args.dataset, label_key=args.label_key,
                        labels_file=args.labels_file, processed=args.processed,
                        n_pcs=args.n_pcs, n_neighbors=args.n_neighbors,
                        n_top_genes=args.n_top_genes, seed=0)

    metrics, assign_df = run_all(data, args)
    if metrics.empty:
        log.error("No runs finished, nothing to save.")
        return
    summary = summarize(metrics, assign_df)

    outdir = Path(args.outdir) / data.name
    save_results(data, metrics, assign_df, summary, args, outdir)
    print_best(summary, data.labels is not None)

    if args.plot:
        from utils.plotting import make_all
        make_all(outdir)


if __name__ == "__main__":
    main()
