import os
import scanpy as sc
import pandas as pd

os.makedirs("results/epithelial_final", exist_ok=True)
os.makedirs("results/final_figures/epithelial_final", exist_ok=True)

clean = sc.read_h5ad("results/epithelial_clean/clean_epithelial_harmony.h5ad")

keep_clusters = ["0","4","10","14","17","19","20"]

final = clean[clean.obs["clean_epi_leiden"].astype(str).isin(keep_clusters)].copy()

print("Clean epithelial cells:", clean.n_obs)
print("Final epithelial-state cells:", final.n_obs)
print("Kept clusters:", ",".join(keep_clusters))

sc.pp.neighbors(final, use_rep="X_pca_harmony", n_neighbors=15)
sc.tl.umap(final)
sc.tl.leiden(final, resolution=0.25, key_added="final_epi_state")

sc.settings.figdir = "results/final_figures/epithelial_final"

sc.pl.umap(final, color="final_epi_state", save="_states.png", show=False)
sc.pl.umap(final, color="subtype", save="_subtype.png", show=False)
sc.pl.umap(final, color="condition", save="_condition.png", show=False)

final.write("results/epithelial_final/final_epithelial_states.h5ad")

print("\nFinal states:")
print(final.obs["final_epi_state"].value_counts().sort_index())

print("\nSubtype composition:")
print(pd.crosstab(final.obs["final_epi_state"], final.obs["subtype"], normalize="index").round(3) * 100)

print("Saved final epithelial states.")
