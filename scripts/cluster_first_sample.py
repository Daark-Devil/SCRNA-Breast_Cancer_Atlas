import scanpy as sc

adata = sc.read_h5ad("GSE161529_processed/TN_MH0126_normalized_hvg.h5ad")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")

print(adata)
print(adata.obs["leiden_0.5"].value_counts())

sc.pl.umap(
    adata,
    color=["leiden_0.5", "n_genes_by_counts", "total_counts", "pct_counts_mt"],
    save="_TN_MH0126_clusters_qc.png",
    show=False
)

adata.write("GSE161529_processed/TN_MH0126_clustered.h5ad")
print("Saved clustered object")
