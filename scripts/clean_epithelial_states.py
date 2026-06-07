import os
import scanpy as sc
import pandas as pd

os.makedirs("results/epithelial_clean", exist_ok=True)
os.makedirs("results/final_figures/epithelial_clean", exist_ok=True)

epi = sc.read_h5ad("results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad")

bad = pd.read_csv(
    "results/epithelial_subcluster/possible_contaminated_epithelial_clusters.tsv",
    sep="\t"
)

bad_clusters = bad["cluster"].astype(str).tolist()

clean = epi[~epi.obs["epi_leiden"].astype(str).isin(bad_clusters)].copy()

print("Original epithelial cells:", epi.n_obs)
print("Removed clusters:", ",".join(bad_clusters))
print("Clean epithelial cells:", clean.n_obs)

# Use existing Harmony PCA instead of recomputing from scaled expression
sc.pp.neighbors(
    clean,
    use_rep="X_pca_harmony",
    n_neighbors=15
)

sc.tl.umap(clean)

sc.tl.leiden(
    clean,
    resolution=0.3,
    key_added="clean_epi_leiden"
)

sc.settings.figdir = "results/final_figures/epithelial_clean"

sc.pl.umap(clean, color="clean_epi_leiden", save="_clusters.png", show=False)
sc.pl.umap(clean, color="subtype", save="_subtype.png", show=False)
sc.pl.umap(clean, color="condition", save="_condition.png", show=False)

clean.write("results/epithelial_clean/clean_epithelial_harmony.h5ad")

print("Clean epithelial clusters:")
print(clean.obs["clean_epi_leiden"].value_counts().sort_index())

print("Saved clean epithelial object")
