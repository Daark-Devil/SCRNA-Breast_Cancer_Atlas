import scanpy as sc

adata = sc.read_h5ad(
    "results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad"
)

print(adata)
print(adata.X)
print(adata.raw)
print(adata.layers.keys() if adata.layers else "No layers")
