import scanpy as sc
import pandas as pd
import numpy as np

sample = "GSM4909281_TN-MH0126"
adata = sc.read_h5ad(f"results/{sample}/h5ad/{sample}_annotated.h5ad")

marker_sets = {
    "T_cell": ["CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "IL7R", "LTB"],
    "Activated_Tcell": ["CXCL13", "TIGIT", "PDCD1", "CTLA4", "IFNG", "GZMB", "CCL5"],
    "NK_cell": ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "FCGR3A"],
    "Macrophage": ["LYZ", "TYROBP", "FCER1G", "AIF1", "LST1", "C1QA", "C1QB", "C1QC"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74"],
    "Plasma_cell": ["MZB1", "JCHAIN", "XBP1", "IGHG1", "IGHG2", "IGHG3", "IGKC"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "COL6A3"],
    "Endothelial": ["PECAM1", "VWF", "ENG", "EGFL7", "KDR", "COL15A1", "GNG11"],
    "Tumor_Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "KRT7", "S100P", "MUC1", "CD24"],
    "Dendritic_APC": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "FCER1A", "CLEC10A", "LAMP3"]
}

# Score each cell for each marker set
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

score_cols = [f"{ct}_score" for ct in marker_sets.keys()]

cluster_scores = (
    adata.obs
    .groupby("leiden_0.5")[score_cols]
    .mean()
)

# Pick highest scoring marker set per cluster
cluster_to_celltype = {}
for cluster in cluster_scores.index:
    best_score_col = cluster_scores.loc[cluster].idxmax()
    best_celltype = best_score_col.replace("_score", "")
    cluster_to_celltype[str(cluster)] = best_celltype

print("\nCluster score matrix:")
print(cluster_scores)

print("\nAuto cluster labels:")
for k, v in cluster_to_celltype.items():
    print(k, "->", v)

adata.obs["auto_cell_type"] = adata.obs["leiden_0.5"].astype(str).map(cluster_to_celltype)

composition = adata.obs["auto_cell_type"].value_counts().reset_index()
composition.columns = ["cell_type", "cell_count"]
composition["percent"] = 100 * composition["cell_count"] / composition["cell_count"].sum()

composition.to_csv(
    f"results/{sample}/tables/{sample}_AUTO_celltype_composition.tsv",
    sep="\t",
    index=False
)

cluster_scores.to_csv(
    f"results/{sample}/tables/{sample}_cluster_celltype_scores.tsv",
    sep="\t"
)

sc.settings.figdir = f"results/{sample}/figures"

sc.pl.umap(
    adata,
    color="auto_cell_type",
    legend_loc="right margin",
    save=f"_{sample}_AUTO_celltypes.png",
    show=False
)

adata.write(f"results/{sample}/h5ad/{sample}_auto_annotated.h5ad")

print("\nAuto composition:")
print(composition)
print("\nSaved auto-annotated outputs.")
