"""
Comparison plots made ONLY from the saved result files, so they can be
re-made at any time without re-running the clustering:

    metrics.csv, summary.csv, assignments.csv, embedding.csv

Figures (saved in <results_dir>/figures/):
    1_param_sensitivity.png    metrics vs each method's parameter (mean +/- std over seeds)
    2_metrics_vs_n_clusters.png ARI / NMI / Silhouette vs number of clusters found
    3_best_comparison.png       best setting of each method side by side
    4_runtime.png               runtime per method
    5_seed_variation.png        ARI spread across seeds at the best setting
    6_umap_comparison.png       UMAP coloured by reference labels and each method
    7_contingency_<method>.png  reference cell types vs clusters (best run)
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

METHOD_COLORS = {"kmeans": "#E76F51", "leiden": "#2A9D8F", "sc3s": "#D08C1F"}
METHOD_LABELS = {"kmeans": "K-means", "leiden": "Leiden", "sc3s": "SC3s"}
METRICS = [("ari", "ARI"), ("nmi", "NMI"), ("silhouette", "Silhouette")]


def _color(m):
    return METHOD_COLORS.get(m, "gray")


def _label(m):
    return METHOD_LABELS.get(m, m)


def load_results(results_dir):
    d = Path(results_dir)
    return dict(
        metrics=pd.read_csv(d / "metrics.csv"),
        summary=pd.read_csv(d / "summary.csv"),
        assignments=pd.read_csv(d / "assignments.csv", index_col=0, dtype=str),
        embedding=pd.read_csv(d / "embedding.csv", index_col=0),
    )


def best_settings(summary, metrics, select_by="ari"):
    """Best parameter value per method, plus the run_id of its first seed."""
    key = f"{select_by}_mean"
    if summary[key].isna().all():  # no labels -> fall back to silhouette
        key = "silhouette_mean"
    best = summary.loc[summary.groupby("method")[key].idxmax()].copy()
    run_ids = []
    for _, row in best.iterrows():
        runs = metrics[(metrics.method == row.method) &
                       np.isclose(metrics.param_value, row.param_value)].sort_values("seed")
        run_ids.append(runs.run_id.iloc[0])
    best["run_id"] = run_ids
    return best.reset_index(drop=True), key.replace("_mean", "")


def plot_param_sensitivity(summary, out):
    methods = list(summary.method.unique())
    rows = METRICS + [("n_clusters", "Number of clusters")]
    fig, axes = plt.subplots(len(rows), len(methods), figsize=(4.2 * len(methods), 3 * len(rows)),
                             squeeze=False, sharey="row")
    for j, m in enumerate(methods):
        s = summary[summary.method == m].sort_values("param_value")
        for i, (key, title) in enumerate(rows):
            ax = axes[i, j]
            ax.errorbar(s.param_value, s[f"{key}_mean"], yerr=s[f"{key}_std"].fillna(0),
                        marker="o", color=_color(m), capsize=3)
            if i == 0:
                ax.set_title(_label(m))
            if i == len(rows) - 1:
                ax.set_xlabel(s.param_name.iloc[0])
            if j == 0:
                ax.set_ylabel(title)
            ax.grid(alpha=0.3)
    fig.suptitle("Parameter sensitivity (mean ± std over seeds)", y=1.0)
    fig.tight_layout()
    fig.savefig(out / "1_param_sensitivity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_metrics_vs_n_clusters(summary, out):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, (key, title) in zip(axes, METRICS):
        for m in summary.method.unique():
            s = summary[summary.method == m].sort_values("n_clusters_mean")
            ax.plot(s.n_clusters_mean, s[f"{key}_mean"], marker="o",
                    color=_color(m), label=_label(m))
        ax.set_xlabel("Number of clusters found")
        ax.set_ylabel(title)
        ax.grid(alpha=0.3)
    axes[0].legend()
    fig.suptitle("Metrics vs number of clusters (all methods on one axis)")
    fig.tight_layout()
    fig.savefig(out / "2_metrics_vs_n_clusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_best_comparison(best, selected_by, out):
    cols = [("ari", "ARI"), ("nmi", "NMI"), ("silhouette", "Silhouette"),
            ("stability", "Stability\n(ARI across seeds)")]
    x = np.arange(len(cols))
    width = 0.8 / len(best)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for i, (_, row) in enumerate(best.iterrows()):
        vals = [row.ari_mean, row.nmi_mean, row.silhouette_mean, row.stability_ari]
        errs = [row.ari_std, row.nmi_std, row.silhouette_std, 0]
        lab = f"{_label(row.method)} ({row.param_name}={row.param_value:g})"
        ax.bar(x + i * width, vals, width, yerr=np.nan_to_num(errs), capsize=3,
               color=_color(row.method), label=lab)
    ax.set_xticks(x + width * (len(best) - 1) / 2)
    ax.set_xticklabels([c[1] for c in cols])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title(f"Best setting of each method (selected by {selected_by.upper()})")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "3_best_comparison.png", dpi=150)
    plt.close(fig)


def plot_runtime(metrics, out):
    methods = list(metrics.method.unique())
    data = [metrics[metrics.method == m].runtime_sec for m in methods]
    fig, ax = plt.subplots(figsize=(6, 4))
    bp = ax.boxplot(data, patch_artist=True)
    ax.set_xticks(range(1, len(methods) + 1))
    ax.set_xticklabels([_label(m) for m in methods])
    for patch, m in zip(bp["boxes"], methods):
        patch.set_facecolor(_color(m))
    ax.set_ylabel("Runtime per run (seconds)")
    ax.set_title("Runtime (all parameter settings and seeds)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "4_runtime.png", dpi=150)
    plt.close(fig)


def plot_seed_variation(metrics, best, selected_by, out):
    data, names, colors = [], [], []
    for _, row in best.iterrows():
        runs = metrics[(metrics.method == row.method) &
                       np.isclose(metrics.param_value, row.param_value)]
        data.append(runs[selected_by])
        names.append(f"{_label(row.method)}\n{row.param_name}={row.param_value:g}")
        colors.append(_color(row.method))
    fig, ax = plt.subplots(figsize=(6, 4))
    bp = ax.boxplot(data, patch_artist=True)
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels(names)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
    for i, d in enumerate(data, start=1):
        ax.scatter(np.full(len(d), i), d, color="black", s=12, zorder=3)
    ax.set_ylabel(selected_by.upper())
    ax.set_title("Variation across random seeds (best setting)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "5_seed_variation.png", dpi=150)
    plt.close(fig)


def _scatter_by_group(ax, emb, groups, title):
    cats = pd.Series(groups).astype(str)
    uniq = sorted(cats.unique(), key=lambda v: (len(v), v))
    cmap = plt.get_cmap("tab20" if len(uniq) > 10 else "tab10")
    for i, c in enumerate(uniq):
        mask = (cats == c).to_numpy()
        ax.scatter(emb.UMAP1[mask], emb.UMAP2[mask], s=8, color=cmap(i % cmap.N), label=c)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    return len(uniq)


def plot_umaps(embedding, assignments, best, out):
    has_labels = "true_label" in assignments.columns
    panels = ([("true_label", "Reference labels")] if has_labels else [])
    for _, row in best.iterrows():
        panels.append((row.run_id, f"{_label(row.method)} ({row.param_name}={row.param_value:g}, "
                                   f"{int(round(row.n_clusters_mean))} clusters)"))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.6), squeeze=False)
    for ax, (col, title) in zip(axes[0], panels):
        n = _scatter_by_group(ax, embedding, assignments[col].to_numpy(), title)
        if col == "true_label" or n <= 20:
            ax.legend(fontsize=6, markerscale=1.5, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(out / "6_umap_comparison.png", dpi=150)
    plt.close(fig)


def plot_contingency(assignments, best, out):
    if "true_label" not in assignments.columns:
        return
    for _, row in best.iterrows():
        tab = pd.crosstab(assignments["true_label"], assignments[row.run_id])
        tab = tab[sorted(tab.columns, key=lambda v: (len(v), v))]
        frac = tab.div(tab.sum(axis=1), axis=0)  # fraction of each cell type in each cluster
        fig, ax = plt.subplots(figsize=(0.55 * tab.shape[1] + 4, 0.45 * tab.shape[0] + 1.5))
        im = ax.imshow(frac.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(tab.shape[1]))
        ax.set_xticklabels(tab.columns)
        ax.set_yticks(range(tab.shape[0]))
        ax.set_yticklabels(tab.index, fontsize=8)
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                if tab.values[i, j] > 0:
                    ax.text(j, i, tab.values[i, j], ha="center", va="center", fontsize=7,
                            color="white" if frac.values[i, j] > 0.5 else "black")
        ax.set_xlabel(f"{_label(row.method)} cluster")
        ax.set_title(f"Reference cell types vs {_label(row.method)} clusters "
                     f"({row.param_name}={row.param_value:g}); colour = fraction of row")
        fig.colorbar(im, ax=ax, fraction=0.03)
        fig.tight_layout()
        fig.savefig(out / f"7_contingency_{row.method}.png", dpi=150)
        plt.close(fig)


def make_all(results_dir, select_by="ari"):
    results_dir = Path(results_dir)
    out = results_dir / "figures"
    out.mkdir(exist_ok=True)
    r = load_results(results_dir)
    best, selected_by = best_settings(r["summary"], r["metrics"], select_by)

    plot_param_sensitivity(r["summary"], out)
    plot_metrics_vs_n_clusters(r["summary"], out)
    plot_best_comparison(best, selected_by, out)
    plot_runtime(r["metrics"], out)
    plot_seed_variation(r["metrics"], best, selected_by, out)
    plot_umaps(r["embedding"], r["assignments"], best, out)
    plot_contingency(r["assignments"], best, out)

    best.to_csv(results_dir / "best_settings.csv", index=False)
    print(f"Figures saved in {out.resolve()}")
