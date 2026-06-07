import scanpy as sc

adata = sc.read_h5ad(
    "GSE161529_processed/TN_MH0126_qc_scrublet.h5ad"
)

print("Before filtering:", adata.shape)

adata = adata[
    (adata.obs["n_genes_by_counts"] >= 200) &
    (adata.obs["n_genes_by_counts"] <= 7000) &
    (adata.obs["pct_counts_mt"] <= 15) &
    (~adata.obs["predicted_doublet"])
].copy()

print("After filtering:", adata.shape)

adata.write(
    "GSE161529_processed/TN_MH0126_filtered.h5ad"
)

print("Saved filtered object")
