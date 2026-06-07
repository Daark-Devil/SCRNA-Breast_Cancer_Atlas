import scanpy as sc
import scrublet as scr
import pandas as pd

adata = sc.read_h5ad("GSE161529_processed/TN_MH0126_qc_metrics.h5ad")

print("Input object:")
print(adata)

# Scrublet expects raw count matrix
scrub = scr.Scrublet(adata.X)

doublet_scores, predicted_doublets = scrub.scrub_doublets()

adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets

print("\nDoublet prediction summary:")
print(adata.obs["predicted_doublet"].value_counts())

print("\nDoublet score summary:")
print(adata.obs["doublet_score"].describe())

adata.write("GSE161529_processed/TN_MH0126_qc_scrublet.h5ad")

adata.obs[["sample_id", "subtype", "condition", "n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score", "predicted_doublet"]].to_csv(
    "results/tables/TN_MH0126_scrublet_qc_table.tsv",
    sep="\t"
)

print("\nSaved:")
print("GSE161529_processed/TN_MH0126_qc_scrublet.h5ad")
print("results/tables/TN_MH0126_scrublet_qc_table.tsv")
