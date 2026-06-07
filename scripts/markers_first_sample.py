import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("GSE161529_processed/TN_MH0126_clustered.h5ad")

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_0.5",
    method="wilcoxon"
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=10,
    sharey=False,
    save="_TN_MH0126_marker_genes.png",
    show=False
)

markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv("results/tables/TN_MH0126_cluster_markers.tsv", sep="\t", index=False)

print(markers.head(30))

adata.write("GSE161529_processed/TN_MH0126_markers.h5ad")
print("Saved marker results")
