import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

os.makedirs("results/integrated_atlas", exist_ok=True)
os.makedirs("results/final_figures/integrated_atlas", exist_ok=True)

MAX_CELLS_PER_SAMPLE = 1500
np.random.seed(42)

meta = pd.read_csv("metadata/sample_metadata.tsv", sep="\t")
samples = meta["sample_id"].tolist()

print("Expected samples:", len(samples), flush=True)

adatas=[]

for sample in samples:
    f = f"results/{sample}/h5ad/{sample}_annotated.h5ad"

    if not os.path.exists(f):
        print("Missing:", f, flush=True)
        continue

    print("Loading:", sample, flush=True)
    a = sc.read_h5ad(f)

    a.obs["cell_type"] = a.obs["cell_type"].replace({
        "Tumor_Epithelial": "Epithelial"
    })

    if a.n_obs > MAX_CELLS_PER_SAMPLE:
        idx = np.random.choice(a.n_obs, MAX_CELLS_PER_SAMPLE, replace=False)
        a = a[idx].copy()

    a.obs_names = a.obs["sample_id"].astype(str) + "_" + a.obs_names.astype(str)
    adatas.append(a)

adata = ad.concat(adatas, join="inner", merge="same")
adata.obs_names_make_unique()

print("Combined:", adata, flush=True)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="sample_id",
    flavor="seurat"
)

adata = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

sc.external.pp.harmony_integrate(
    adata,
    key="sample_id"
)

sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=15,
    n_pcs=30
)

sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="integrated_leiden")

sc.settings.figdir = "results/final_figures/integrated_atlas"

sc.pl.umap(adata, color="cell_type", save="_downsampled_cell_type.png", show=False)
sc.pl.umap(adata, color="subtype", save="_downsampled_subtype.png", show=False)
sc.pl.umap(adata, color="condition", save="_downsampled_condition.png", show=False)
sc.pl.umap(adata, color="integrated_leiden", save="_downsampled_clusters.png", show=False)

adata.write("results/integrated_atlas/all_samples_harmony_downsampled.h5ad")

print("Saved downsampled integrated atlas", flush=True)
