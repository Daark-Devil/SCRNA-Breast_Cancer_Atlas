import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

os.makedirs("results/integrated_atlas", exist_ok=True)
os.makedirs("results/final_figures/integrated_atlas", exist_ok=True)

np.random.seed(42)
MAX_CELLS_PER_SAMPLE = 1500

meta = pd.read_csv("metadata/sample_metadata.tsv", sep="\t")

features = pd.read_csv(
    "GSE161529/GSE161529_features.tsv.gz",
    sep="\t",
    header=None,
    names=["gene_id", "gene_symbol", "feature_type"]
)

adatas = []

for _, row in meta.iterrows():

    sample = row["sample_id"]
    matrix_file = row["matrix_file"]
    barcode_file = row["barcode_file"]

    print("Loading raw:", sample, flush=True)

    if not os.path.exists(matrix_file) or not os.path.exists(barcode_file):
        print("Missing raw files:", sample, flush=True)
        continue

    ann_file = f"results/{sample}/h5ad/{sample}_annotated.h5ad"

    if not os.path.exists(ann_file):
        print("Missing annotated file:", sample, flush=True)
        continue

    a = sc.read_mtx(matrix_file).T

    barcodes = pd.read_csv(
        barcode_file,
        header=None,
        names=["barcode"]
    )

    a.var["gene_id"] = features["gene_id"].astype(str).values
    a.var["gene_symbol"] = features["gene_symbol"].astype(str).values
    a.var["feature_type"] = features["feature_type"].astype(str).values

    a.var_names = a.var["gene_symbol"].astype(str)
    a.var_names_make_unique()

    raw_barcodes = barcodes["barcode"].astype(str).values
    prefixed_barcodes = sample + "_" + raw_barcodes

    a.obs_names = prefixed_barcodes

    for col in [
        "gsm",
        "sample_id",
        "sample_name",
        "patient_id",
        "condition",
        "subtype",
        "tissue_fraction"
    ]:
        a.obs[col] = row[col]

    ann = sc.read_h5ad(ann_file)

    ann.obs["cell_type"] = ann.obs["cell_type"].replace({
        "Tumor_Epithelial": "Epithelial"
    })

    celltype_map = ann.obs["cell_type"].to_dict()

    a.obs["cell_type"] = a.obs_names.map(celltype_map)

    matched = a.obs["cell_type"].notna().sum()
    print("Matched cells:", matched, "of", a.n_obs, flush=True)

    a = a[~a.obs["cell_type"].isna()].copy()

    if a.n_obs == 0:
        print("No matched cells after filtering:", sample, flush=True)
        continue

    if a.n_obs > MAX_CELLS_PER_SAMPLE:
        idx = np.random.choice(a.n_obs, MAX_CELLS_PER_SAMPLE, replace=False)
        a = a[idx].copy()

    adatas.append(a)

adata = ad.concat(
    adatas,
    join="inner",
    merge="same"
)

adata.obs_names_make_unique()

print("Combined atlas object:", adata, flush=True)

sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor="seurat"
)

adata = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)

import harmonypy as hm

print("Running manual Harmony...", flush=True)

harmony_out = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    vars_use=["sample_id"]
)

Z = harmony_out.Z_corr

if Z.shape[0] == adata.n_obs:
    adata.obsm["X_pca_harmony"] = Z
elif Z.shape[1] == adata.n_obs:
    adata.obsm["X_pca_harmony"] = Z.T
else:
    raise ValueError(f"Unexpected Harmony shape: {Z.shape}, cells: {adata.n_obs}")

print("Harmony corrected shape:", adata.obsm["X_pca_harmony"].shape, flush=True)

sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=15
)

sc.tl.umap(adata)

sc.tl.leiden(
    adata,
    resolution=0.5,
    key_added="integrated_leiden"
)

sc.settings.figdir = "results/final_figures/integrated_atlas"

sc.pl.umap(
    adata,
    color="cell_type",
    save="_raw_cell_type.png",
    show=False
)

sc.pl.umap(
    adata,
    color="subtype",
    save="_raw_subtype.png",
    show=False
)

sc.pl.umap(
    adata,
    color="condition",
    save="_raw_condition.png",
    show=False
)

sc.pl.umap(
    adata,
    color="integrated_leiden",
    save="_raw_clusters.png",
    show=False
)

adata.write(
    "results/integrated_atlas/all_samples_raw_harmony_downsampled.h5ad"
)

print("Saved integrated raw atlas", flush=True)
