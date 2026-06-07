import scanpy as sc

adata = sc.read_h5ad(
    "GSE161529_processed/TN_MH0126_clustered.h5ad"
)

markers = [
    "CD3D",      # T cells
    "NKG7",      # NK cells
    "LYZ",       # Macrophages
    "MZB1",      # Plasma cells
    "COL1A1",    # Fibroblasts
    "EPCAM",     # Tumor epithelial
    "EGFL7"      # Endothelial
]

for gene in markers:
    if gene in adata.var_names:
        sc.pl.umap(
            adata,
            color=gene,
            save=f"_{gene}.png",
            show=False
        )

print("Done")
