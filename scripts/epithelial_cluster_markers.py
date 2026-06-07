import os
import scanpy as sc
import pandas as pd

os.makedirs("results/epithelial_subcluster/markers", exist_ok=True)
os.makedirs("results/final_figures/epithelial_subcluster/markers", exist_ok=True)

epi = sc.read_h5ad("results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad")

print("Epithelial cells:", epi.n_obs)
print("Clusters:", epi.obs["epi_leiden"].nunique())

sc.tl.rank_genes_groups(
    epi,
    groupby="epi_leiden",
    method="wilcoxon"
)

markers = sc.get.rank_genes_groups_df(
    epi,
    group=None
)

markers.to_csv(
    "results/epithelial_subcluster/epithelial_cluster_markers_all.tsv",
    sep="\t",
    index=False
)

top10 = (
    markers
    .dropna(subset=["names"])
    .groupby("group")
    .head(10)
)

top10.to_csv(
    "results/epithelial_subcluster/epithelial_cluster_markers_top10.tsv",
    sep="\t",
    index=False
)

print("\nTop 10 markers per epithelial cluster saved.")
print(top10[["group","names","scores","logfoldchanges","pvals_adj"]].head(100).to_string(index=False))

sc.settings.figdir = "results/final_figures/epithelial_subcluster/markers"

top_genes = top10["names"].dropna().unique().tolist()[:40]

sc.pl.dotplot(
    epi,
    var_names=top_genes,
    groupby="epi_leiden",
    save="_top_marker_dotplot.png",
    show=False
)

print("Saved marker tables and dotplot.")
