import glob
import scanpy as sc
import anndata as ad
import pandas as pd
import os

os.makedirs("results/integrated_atlas", exist_ok=True)
os.makedirs("results/final_figures/integrated_atlas", exist_ok=True)

files = sorted(glob.glob("results/*/h5ad/*_annotated.h5ad"))
print("Files:", len(files))

adatas = []

for f in files:
    a = sc.read_h5ad(f)

    a.obs["cell_type"] = a.obs["cell_type"].replace({
        "Tumor_Epithelial": "Epithelial"
    })

    a.obs_names = a.obs["sample_id"].astype(str) + "_" + a.obs_names.astype(str)

    adatas.append(a)

adata = ad.concat(
    adatas,
    join="inner",
    merge="same"
)

adata.obs_names_make_unique()

print("Combined object:")
print(adata)

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

sc.pl.umap(
    adata,
    color=["cell_type"],
    save="_cell_type.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["subtype"],
    save="_subtype.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["condition"],
    save="_condition.png",
    show=False
)

sc.pl.umap(
    adata,
    color=["integrated_leiden"],
    save="_clusters.png",
    show=False
)

adata.write("results/integrated_atlas/all_samples_harmony_integrated.h5ad")

print("Saved integrated atlas")
print("results/integrated_atlas/all_samples_harmony_integrated.h5ad")
