"""
Make all comparison plots from saved results (no clustering is re-run).

Usage (Terminal):
    python make_plots.py                                  # default: results/pbmc68k_reduced
    python make_plots.py results/pbmc68k_reduced --select-by nmi

Usage (Spyder):
    press the green Run button                            # uses the default folder
    %run make_plots.py results/other_dataset              # another results folder
"""

import argparse
from pathlib import Path

from utils.plotting import make_all

PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULTS = PROJECT_DIR / "results" / "pbmc68k_reduced"

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot results saved by main.py")
    p.add_argument("results_dir", nargs="?", default=str(DEFAULT_RESULTS),
                   help="folder with metrics.csv, summary.csv, ... "
                        "(default: results/pbmc68k_reduced)")
    p.add_argument("--select-by", default="ari", choices=["ari", "nmi", "silhouette"],
                   help="metric used to pick the best setting of each method")
    args = p.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_absolute():
        results_dir = PROJECT_DIR / results_dir
    if not (results_dir / "metrics.csv").exists():
        raise SystemExit(f"No results found in {results_dir}. Run main.py first.")

    make_all(results_dir, select_by=args.select_by)
