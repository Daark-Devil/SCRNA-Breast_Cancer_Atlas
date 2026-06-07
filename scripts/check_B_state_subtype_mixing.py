import os
import scanpy as sc
import pandas as pd

os.makedirs("results/epithelial_diagnostics", exist_ok=True)

adata = sc.read_h5ad(
    "results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad"
)

print("Object:")
print(adata)

print("\nSubtype counts:")
print(adata.obs["subtype"].value_counts())

print("\nEpithelial cluster counts:")
print(adata.obs["epi_leiden"].value_counts().sort_index())

count_table = pd.crosstab(
    adata.obs["epi_leiden"],
    adata.obs["subtype"]
)

percent_table = pd.crosstab(
    adata.obs["epi_leiden"],
    adata.obs["subtype"],
    normalize="index"
) * 100

count_table.to_csv(
    "results/epithelial_diagnostics/epi_cluster_subtype_counts.tsv",
    sep="\t"
)

percent_table.to_csv(
    "results/epithelial_diagnostics/epi_cluster_subtype_percent.tsv",
    sep="\t"
)

dominant = percent_table.idxmax(axis=1)
dominant_percent = percent_table.max(axis=1)

summary = pd.DataFrame({
    "epi_cluster": percent_table.index,
    "dominant_subtype": dominant.values,
    "dominant_percent": dominant_percent.values,
    "total_cells": count_table.sum(axis=1).values
})

summary = summary.sort_values(
    ["dominant_subtype", "dominant_percent"],
    ascending=[True, False]
)

summary.to_csv(
    "results/epithelial_diagnostics/epi_cluster_dominant_subtype_summary.tsv",
    sep="\t",
    index=False
)

print("\nCluster subtype percent:")
print(percent_table.round(1).to_string())

print("\nDominant subtype summary:")
print(summary.to_string(index=False))
