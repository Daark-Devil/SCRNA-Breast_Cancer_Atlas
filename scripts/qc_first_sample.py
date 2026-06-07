import scanpy as sc

adata = sc.read_h5ad("GSE161529_processed/TN_MH0126_raw.h5ad")

adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print(adata)
print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe())

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_TN_MH0126_qc_violin.png",
    show=False
)

sc.pl.scatter(
    adata,
    x="total_counts",
    y="pct_counts_mt",
    save="_TN_MH0126_counts_vs_mt.png",
    show=False
)

sc.pl.scatter(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    save="_TN_MH0126_counts_vs_genes.png",
    show=False
)

adata.write("GSE161529_processed/TN_MH0126_qc_metrics.h5ad")

print("Saved QC object")
