import os
import scanpy as sc
import pandas as pd

os.makedirs("results/epithelial_clean/markers", exist_ok=True)
os.makedirs("results/final_figures/epithelial_clean/markers", exist_ok=True)

clean = sc.read_h5ad("results/epithelial_clean/clean_epithelial_harmony.h5ad")

print("Clean epithelial cells:", clean.n_obs)
print("Clean clusters:", clean.obs["clean_epi_leiden"].nunique())

sc.tl.rank_genes_groups(
    clean,
    groupby="clean_epi_leiden",
    method="wilcoxon"
)

markers = sc.get.rank_genes_groups_df(clean, group=None)

markers.to_csv(
    "results/epithelial_clean/markers/clean_epithelial_cluster_markers_all.tsv",
    sep="\t",
    index=False
)

top10 = (
    markers
    .dropna(subset=["names"])
    .sort_values(["group","scores"], ascending=[True,False])
    .groupby("group")
    .head(10)
)

top10.to_csv(
    "results/epithelial_clean/markers/clean_epithelial_cluster_markers_top10.tsv",
    sep="\t",
    index=False
)

print(top10[["group","names","scores","logfoldchanges","pvals_adj"]].head(120).to_string(index=False))
print("Saved clean epithelial marker tables.")
