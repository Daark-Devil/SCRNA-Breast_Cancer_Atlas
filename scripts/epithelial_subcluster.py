import os
import scanpy as sc
import harmonypy as hm
import numpy as np

os.makedirs("results/epithelial_subcluster", exist_ok=True)
os.makedirs("results/final_figures/epithelial_subcluster", exist_ok=True)

adata = sc.read_h5ad("results/integrated_atlas/all_samples_raw_harmony_downsampled.h5ad")

epi = adata[adata.obs["cell_type"] == "Epithelial"].copy()

print("Epithelial cells:", epi.n_obs, flush=True)

sc.pp.highly_variable_genes(epi, n_top_genes=2000, flavor="seurat")
epi = epi[:, epi.var["highly_variable"]].copy()

sc.pp.scale(epi, max_value=10)
sc.tl.pca(epi, n_comps=40)

print("Running epithelial Harmony...", flush=True)

harmony_out = hm.run_harmony(
    epi.obsm["X_pca"],
    epi.obs,
    vars_use=["sample_id"]
)

Z = harmony_out.Z_corr

if Z.shape[0] == epi.n_obs:
    epi.obsm["X_pca_harmony"] = Z
elif Z.shape[1] == epi.n_obs:
    epi.obsm["X_pca_harmony"] = Z.T
else:
    raise ValueError(f"Unexpected Harmony shape: {Z.shape}")

sc.pp.neighbors(epi, use_rep="X_pca_harmony", n_neighbors=15)
sc.tl.umap(epi)
sc.tl.leiden(epi, resolution=0.6, key_added="epi_leiden")

sc.settings.figdir = "results/final_figures/epithelial_subcluster"

sc.pl.umap(epi, color="subtype", save="_subtype.png", show=False)
sc.pl.umap(epi, color="condition", save="_condition.png", show=False)
sc.pl.umap(epi, color="sample_id", save="_sample.png", show=False)
sc.pl.umap(epi, color="epi_leiden", save="_clusters.png", show=False)

epi.write("results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad")

print("Saved epithelial subcluster", flush=True)
