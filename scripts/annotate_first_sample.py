import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("GSE161529_processed/TN_MH0126_markers.h5ad")

mapping = {
    "0": "Activated_Tcell",
    "1": "Tumor_Epithelial",
    "2": "Tumor_Epithelial",
    "3": "Macrophage",
    "4": "NK_cell",
    "5": "T_cell",
    "6": "Plasma_cell",
    "7": "Fibroblast",
    "8": "Dendritic_APC",
    "9": "Endothelial",
}

adata.obs["cell_type"] = adata.obs["leiden_0.5"].astype(str).map(mapping)

print("Cell type counts:")
print(adata.obs["cell_type"].value_counts())

composition = adata.obs["cell_type"].value_counts().reset_index()
composition.columns = ["cell_type", "cell_count"]
composition["percent"] = 100 * composition["cell_count"] / composition["cell_count"].sum()

composition.to_csv(
    "results/tables/TN_MH0126_celltype_composition.tsv",
    sep="\t",
    index=False
)

sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="right margin",
    save="_TN_MH0126_celltypes.png",
    show=False
)

adata.write("GSE161529_processed/TN_MH0126_annotated.h5ad")

print("Saved:")
print("GSE161529_processed/TN_MH0126_annotated.h5ad")
print("results/tables/TN_MH0126_celltype_composition.tsv")
print("figures/umap_TN_MH0126_celltypes.png")
