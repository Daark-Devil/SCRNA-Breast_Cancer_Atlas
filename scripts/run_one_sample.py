import argparse
import os
import scanpy as sc
import pandas as pd
import scrublet as scr

parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument("--raw_dir", default="GSE161529")
parser.add_argument("--features", default="GSE161529/GSE161529_features.tsv.gz")
parser.add_argument("--min_genes", type=int, default=200)
parser.add_argument("--max_genes", type=int, default=7000)
parser.add_argument("--max_mt", type=float, default=15.0)
args = parser.parse_args()

sample = args.sample
outdir = f"results/{sample}"
figdir = f"{outdir}/figures"
tabdir = f"{outdir}/tables"
h5dir = f"{outdir}/h5ad"

os.makedirs(figdir, exist_ok=True)
os.makedirs(tabdir, exist_ok=True)
os.makedirs(h5dir, exist_ok=True)
sc.settings.figdir = figdir

print(f"=== Processing {sample} ===")

matrix_file = f"{args.raw_dir}/{sample}-matrix.mtx.gz"
barcode_file = f"{args.raw_dir}/{sample}-barcodes.tsv.gz"

adata = sc.read_mtx(matrix_file).T

features = pd.read_csv(
    args.features,
    sep="\t",
    header=None,
    names=["gene_id", "gene_symbol", "feature_type"]
)

barcodes = pd.read_csv(
    barcode_file,
    header=None,
    names=["barcode"]
)

adata.var["gene_id"] = features["gene_id"].astype(str).values
adata.var["gene_symbol"] = features["gene_symbol"].astype(str).values
adata.var["feature_type"] = features["feature_type"].astype(str).values

adata.var_names = adata.var["gene_symbol"].astype(str)
adata.var_names_make_unique()
adata.var.index.name = None

adata.obs_names = sample + "_" + barcodes["barcode"].astype(str)
adata.obs.index.name = None
adata.obs["sample_id"] = sample

meta = pd.read_csv("metadata/sample_metadata.tsv", sep="\t")
row = meta[meta["sample_id"] == sample]

if len(row) == 1:
    row = row.iloc[0]
    for col in ["gsm", "sample_name", "patient_id", "condition", "subtype", "tissue_fraction"]:
        adata.obs[col] = row[col]
else:
    adata.obs["condition"] = "Unknown"
    adata.obs["subtype"] = "Unknown"
    adata.obs["tissue_fraction"] = "Unknown"

adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs["doublet_score"] = doublet_scores

if predicted_doublets is None:
    print("Scrublet automatic threshold failed. Using fallback doublet_score > 0.25")
    predicted_doublets = doublet_scores > 0.25

adata.obs["predicted_doublet"] = pd.Series(
    predicted_doublets,
    index=adata.obs_names
).astype(bool)

adata.obs.to_csv(f"{tabdir}/{sample}_qc_scrublet.tsv", sep="\t")

before_cells = adata.n_obs

adata = adata[
    (adata.obs["n_genes_by_counts"] >= args.min_genes) &
    (adata.obs["n_genes_by_counts"] <= args.max_genes) &
    (adata.obs["pct_counts_mt"] <= args.max_mt) &
    (~adata.obs["predicted_doublet"])
].copy()

after_cells = adata.n_obs

print("Before filtering:", before_cells)
print("After filtering:", after_cells)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor="seurat"
)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_0.5",
    method="wilcoxon"
)

markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(f"{tabdir}/{sample}_cluster_markers.tsv", sep="\t", index=False)

marker_sets = {
    "T_cell": ["CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "IL7R", "LTB"],
    "Activated_Tcell": ["CXCL13", "TIGIT", "PDCD1", "CTLA4", "IFNG", "GZMB", "CCL5"],
    "NK_cell": ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "FCGR3A"],
    "Macrophage": ["LYZ", "TYROBP", "FCER1G", "AIF1", "LST1", "C1QA", "C1QB", "C1QC"],
    "B_cell_APC": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA", "HLA-DPA1"],
    "Plasma_cell": ["MZB1", "JCHAIN", "XBP1", "IGHG1", "IGHG2", "IGHG3", "IGKC"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "COL6A3"],
    "Endothelial": ["PECAM1", "VWF", "ENG", "EGFL7", "KDR", "COL15A1", "GNG11"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "KRT7", "S100P", "MUC1", "CD24"],
    "Dendritic_APC": ["FCER1A", "CLEC10A", "LAMP3", "ITGAX", "CD1C", "HLA-DRA", "CD74"]
}

for celltype, genes in marker_sets.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) >= 2:
        sc.tl.score_genes(
            adata,
            gene_list=genes_present,
            score_name=f"{celltype}_score",
            use_raw=False
        )
    else:
        adata.obs[f"{celltype}_score"] = 0

score_cols = [f"{ct}_score" for ct in marker_sets]

cluster_scores = (
    adata.obs
    .groupby("leiden_0.5", observed=True)[score_cols]
    .mean()
)

cluster_scores.to_csv(
    f"{tabdir}/{sample}_cluster_celltype_scores.tsv",
    sep="\t"
)

cluster_to_celltype = {}
for cluster in cluster_scores.index:
    best = cluster_scores.loc[cluster].idxmax().replace("_score", "")
    cluster_to_celltype[str(cluster)] = best

pd.DataFrame(
    [{"cluster": k, "auto_cell_type": v} for k, v in cluster_to_celltype.items()]
).to_csv(
    f"{tabdir}/{sample}_cluster_auto_labels.tsv",
    sep="\t",
    index=False
)

adata.obs["cell_type"] = (
    adata.obs["leiden_0.5"].astype(str)
    .map(cluster_to_celltype)
    .fillna("Unknown")
)

composition = adata.obs["cell_type"].value_counts().reset_index()
composition.columns = ["cell_type", "cell_count"]
composition["percent"] = 100 * composition["cell_count"] / composition["cell_count"].sum()

composition.to_csv(
    f"{tabdir}/{sample}_celltype_composition.tsv",
    sep="\t",
    index=False
)

sc.pl.umap(
    adata,
    color=["leiden_0.5", "cell_type", "subtype", "condition", "pct_counts_mt"],
    save=f"_{sample}_summary.png",
    show=False
)

adata.write(f"{h5dir}/{sample}_annotated.h5ad")

summary = pd.DataFrame([{
    "sample": sample,
    "raw_cells": before_cells,
    "filtered_cells": after_cells,
    "retention_percent": round(100 * after_cells / before_cells, 2),
    "genes": adata.n_vars,
    "clusters": adata.obs["leiden_0.5"].nunique()
}])

summary.to_csv(
    f"{tabdir}/{sample}_sample_summary.tsv",
    sep="\t",
    index=False
)

print("Done:", sample)
print(composition)
