import scanpy as sc

adata = sc.read_h5ad("GSE161529_processed/TN_MH0126_filtered.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor="seurat"
)

print(adata)
print(adata.var["highly_variable"].value_counts())

adata.write("GSE161529_processed/TN_MH0126_normalized_hvg.h5ad")
print("Saved normalized + HVG object")
