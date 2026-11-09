# Clustering Methods for Single Cell RNA-Sequencing Data

A reproducible comparison of three clustering methods (K-means, Leiden, and SC3s) on the same single-cell RNA seq dataset. We measure how well each method recovers known cell types, how well separated its clusters are, how many clusters it finds, how long it takes, and how much its results change with key parameters and random seeds.

## Contents

- Background
- Methods compared
- Dataset
- Repository Structure
- Setup
- How to run
- Shared file formats
- Evaluation
- Team and responsibilities
- Timeline
- References

### Background

scRNA-seq measures gene expression in individual cells. Clustering groups cells with similar expression profiles, and these groups are then interpreted as cell types. Different clustering methods make different assumptions, so they can produce different groupings, with different stability and computational cost.

Aim

1. Compare agreement between each method's clusters and the reference cell-type labels.
2. Compare cluster separation, runtime and the number of clusters each method finds
3. Examine how results change the key parameters (number of cluster K, Leiden resolution) and random seeds.

### Methods compared

| Column  | Description                                                    |         Key Parameters |
| :------ | -------------------------------------------------------------- | ---------------------: |
| K-means | Classical centroid-based baseline, run in PCA space            | k (number of clusters) |
| Leiden  | Graph-based community detection on a k-nearest neighbour graph |             Resolution |
| SC3s    | Consensus clustering designed for scRNA-seq                    | k (number of clusters) |

All methods receive the same PCA representation.

### Dataset

PBMC68k reduced, bundled with Scanpy:

```Python
import scanpy as sc
adata = sc.datasets.pbmc68k_reduced()
```

- 700 cells × 765 genes; already normalized and scaled (do not normalize or scale again)

### Repository Structure

scrna-clustering-comparison/
├── README.md
├── environment.yml # Conda environment (all dependencies, pinned)
├── config.yaml # Shared settings: PCs, neighbours, seeds, k values, resolutions
├── .gitignore
│
├── src/
│ ├── preprocess.py # Load data, PCA, neighbour graph, UMAP → data/pbmc_pca.h5ad
│ ├── run_leiden.py # Leiden runs and resolution sweep
│ ├── run_kmeans.py # K-means runs and k sweep
│ ├── run_sc3s.py # SC3s runs and k sweep
│ ├── evaluate.py # ARI, NMI, silhouette, runtime → results/metrics.csv
│ └── make_figures.py # All final figures → figures/
│
├── data/ # Generated; not committed
│ └── pbmc_pca.h5ad
│
├── results/
│ ├── clusters_leiden.csv
│ ├── clusters_kmeans.csv
│ ├── clusters_sc3s.csv
│ └── metrics.csv
│
├── figures/ # Final figures used in the report and slides
├── notebooks/ # Exploration only; final results come from src/
└── report/ # Report drafts and final PDF

Each member does this once on their own laptop, and once in their own account on the team VM (see Running on the team VM).

1. Install the prerequisites
   Tool | Windows | Mac |Linux / VM
   | :---- | --------------------------------------------- |----------------:|
   Git | git-scm.com (includes Git Bash) |xcode-select --install | sudo apt install -y git
   Conda | Miniforge or Anaconda Miniforge Miniforge (see VM section)
   Editor | VS Code + Python and Jupyter extensions | same same

2. Connect your shell to conda (once per machine)

conda activate only works after this step. Run the line for the terminal you actually use:

| Terminal                        |               Command |
| :------------------------------ | --------------------: |
| Git Bash (prompt shows MINGW64) |       conda init bash |
| PowerShell                      | conda init powershell |
| Command Prompt                  |    conda init cmd.exe |
| Mac (zsh)                       |        conda init zsh |
| Linux / VM (bash)               |       conda init bash |

Then close and reopen the terminal. Check with conda --version.

3. Set up GitHub access with SSH (once per machine)

GitHub does not accept passwords for Git. Create an SSH key and add it to your account:

```bash
ssh-keygen -t ed25519 -C "your-email@uconn.edu"   # press Enter at each prompt
cat ~/.ssh/id_ed25519.pub                         # copy this line
```

Paste it into GitHub → Settings → SSH and GPG keys → New SSH key, then test:

```bash
ssh -T git@github.com    # should greet you by your GitHub username
4. Clone the repository
bash
git clone git@github.com:<owner>/3800-team-7.git
cd 3800-team-7
```

5. Create and activate the environment

```bash
conda env create -f environment.yml
conda activate 3800-team-7
```

Your prompt should now start with (scrna-clustering). Check that everything imports:

```bash
python -c "import scanpy, sklearn, leidenalg, sc3s; print('scanpy', scanpy.__version__)"
```

6. Use the environment in VS Code
   Scripts: Ctrl+Shift+P → Python: Select Interpreter → choose scrna-clustering.
   Notebooks: open a .ipynb, click Select Kernel (top right) → choose scrna-clustering.

7. Keep the environment up to date

When someone changes environment.yml (via a merged pull request), everyone runs:

```bash
git pull
conda env update -f environment.yml --prune
```

Rule: never pip install or conda install packages directly for project work. Add the package to environment.yml in a pull request and tell the team.

### How to run

Run the steps in order from the repository root. All settings are read from config.yaml

```bash
conda activate 3800-team-7

# 1. Preprocess: load data, PCA, neighbour graph, UMAP
python src/preprocess.py

# 2. Run each clustering method (default runs, parameter sweeps and seed repeats)
python src/run_leiden.py
python src/run_kmeans.py
python src/run_sc3s.py

# 3. Score every run
python src/evaluate.py

# 4. Make all figures
python src/make_figures.py
```

Run step 1 first. data/ is not in Git, so after cloning, everyone generates data/pbmc_pca.h5ad locally with preprocess.py. With the same config.yaml, everyone gets the same file.
Each later step can be run on its own once the previous step's output exists.
The dataset is small: each step takes seconds to a few minutes on a laptop.
Laptops are for writing and testing code. Official results (the numbers and figures in the report) are produced on the team VM.
Running on the team VM

Each member uses their own account on the VM and their own clone of the repository.

One-time setup:

```bash
# Install Miniforge in your home folder
cd ~
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
source ~/.bashrc

# Set your Git identity, add an SSH key to GitHub (Setup step 3), then clone and create the environment
git config --global user.name "Your Name"
git config --global user.email "your-email@uconn.edu"
git clone git@github.com:<owner>/3800-team-7.git
cd 3800-team-7
conda env create -f environment.yml
```

Each session:

```bash
tmux new -s work               # or: tmux attach -t work   (resume a session)
cd ~/3800-team-7
conda activate 3800-team-7
git checkout main && git pull
python src/preprocess.py       # if data/pbmc_pca.h5ad doesn't exist yet
python src/run_kmeans.py       # etc.
```

tmux keeps jobs running if your connection drops (Ctrl+B then D to detach).
Use sudo apt only for system tools; Python packages go in environment.yml.

### Shared settings

All members use these values. Change them only after the group agrees.

```yaml
n_pcs: 30
n_neighbors: 15
seeds: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
k_values: [4, 6, 8, 10, 12, 14]
leiden_resolutions: [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
label_column: bulk_labels
data_path: data/pbmc_pca.h5ad
```

### Shared file formats

Every method saves result in the same format, so one script can evaluate all of them.

** Cluster assignments **: `results/clusters_<method>.csv` one row per cell per run:

| Column    |                       Description |
| :-------- | --------------------------------: |
| cell_id   | Cell barcode from adata.obs_names |
| Cluster   |    Cluster assigned by the method |
| method    |            Kmeans, Leiden or SC3S |
| parameter |   Value of k or Leiden resolution |
| seed      |                  Random seed used |

** Metrics **: `results/metrics.csv` one row per cell per run:

| Column     |                                        Description |
| :--------- | -------------------------------------------------: |
| method     |                             kmeans, leiden or sc3s |
| parameter  |                    value of k or leiden resolution |
| seed       |                                   Random seed used |
| n_clusters |                           Number of clusters found |
| ARI        |            Adjusted Rand Index vs reference labels |
| NMI        | Normalized mutual informations vs reference labels |
| Silhouette |                      Silhouette score in PCA space |
| runtime    |                                 runtime in seconds |

### Best Practices

Make a branch for each task, not one per member. Start each one from the latest main and merge it back with a pull request when the task is done.

```bash
git checkout master
git pull                          # get everyone's latest work first
git checkout -b leiden-first-run  # new branch for this task
# ...work, then commit often:
git add .
git commit -m "Run Leiden at default resolution (#6)"
git push -u origin leiden-first-run
```

Then on GitHub:

Open a pull request into main, with Closes #6 in the description.
Review and merges it.
Delete the branch.

Good habits:

Name branches after the task: preprocess, kmeans-first-run, evaluate-script, sc3s-k-sweep.
Pull main before starting each new branch, so your work includes your teammates’ merged changes.
Keep branches short, a few days at most.
Don’t edit someone else’s files on your branch unless you’ve agreed to. Each method has its own script (run_leiden.py, run_kmeans.py, run_sc3s.py), so conflicts should be rare.
Coordinate edits to shared files like config.yaml and environment.yml in the group
Run git status before every commit and make sure no data files are included.
