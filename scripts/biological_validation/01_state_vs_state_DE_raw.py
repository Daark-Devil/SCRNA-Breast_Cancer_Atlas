import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

OUTDIR = "results/biological_validation/state_comparisons"
os.makedirs(OUTDIR, exist_ok=True)

STATE_MAP = {
    "37": "ER_enriched_37",
    "22": "Normal_enriched_22",
    "33": "TNBC_enriched_33",
    "45": "HER2_enriched_45",
    "9":  "TNBC_enriched_9",
    "18": "HER2_enriched_18"
}

COMPARISONS = [
    ("ER_enriched_37", "Normal_enriched_22"),
    ("TNBC_enriched_33", "Normal_enriched_22"),
    ("TNBC_enriched_9", "Normal_enriched_22"),
    ("HER2_enriched_45", "Normal_enriched_22"),
    ("HER2_enriched_18", "Normal_enriched_22"),
    ("TNBC_enriched_33", "ER_enriched_37"),
    ("HER2_enriched_45", "ER_enriched_37"),
    ("HER2_enriched_18", "ER_enriched_37"),
    ("HER2_enriched_45", "TNBC_enriched_33")
]

print("Loading epithelial state object...", flush=True)

epi = sc.read_h5ad(
    "results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad"
)

epi.obs["epi_leiden"] = epi.obs["epi_leiden"].astype(str)

selected = epi[epi.obs["epi_leiden"].isin(STATE_MAP.keys())].copy()
selected.obs["state_label"] = selected.obs["epi_leiden"].map(STATE_MAP)

print("Selected cells:", selected.n_obs, flush=True)
print(selected.obs["state_label"].value_counts(), flush=True)

state_map = selected.obs["state_label"].to_dict()

meta = pd.read_csv("metadata/sample_metadata.tsv", sep="\t")

features = pd.read_csv(
    "GSE161529/GSE161529_features.tsv.gz",
    sep="\t",
    header=None,
    names=["gene_id", "gene_symbol", "feature_type"]
)

adatas = []

for _, row in meta.iterrows():

    sample = row["sample_id"]
    matrix_file = row["matrix_file"]
    barcode_file = row["barcode_file"]

    if not os.path.exists(matrix_file) or not os.path.exists(barcode_file):
        continue

    barcodes = pd.read_csv(
        barcode_file,
        header=None,
        names=["barcode"]
    )

    raw_obs_names = sample + "_" + barcodes["barcode"].astype(str)

    keep_mask = raw_obs_names.isin(state_map.keys())

    if keep_mask.sum() == 0:
        continue

    print("Processing raw sample:", sample, "cells:", keep_mask.sum(), flush=True)

    a = sc.read_mtx(matrix_file).T

    a.var["gene_id"] = features["gene_id"].astype(str).values
    a.var["gene_symbol_original"] = features["gene_symbol"].astype(str).values
    a.var["feature_type"] = features["feature_type"].astype(str).values

    a.var_names = features["gene_symbol"].astype(str).values
    a.var_names_make_unique()
    a.var.index.name = None

    a.obs_names = raw_obs_names.values

    a = a[keep_mask.values].copy()

    a.obs["sample_id"] = sample
    a.obs["subtype"] = row["subtype"]
    a.obs["condition"] = row["condition"]
    a.obs["patient_id"] = row["patient_id"]
    a.obs["state_label"] = a.obs_names.map(state_map)

    adatas.append(a)

adata = ad.concat(adatas, join="inner", merge="same")
adata.obs_names_make_unique()
adata.var.index.name = None

print("Raw state object:", adata, flush=True)

print("\nState counts:")
print(adata.obs["state_label"].value_counts(), flush=True)

print("\nSubtype counts:")
print(adata.obs["subtype"].value_counts(), flush=True)

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write(f"{OUTDIR}/state_anchor_log_normalized.h5ad")

summary_rows = []

for group, ref in COMPARISONS:

    print("\n==============================", flush=True)
    print(f"Running comparison: {group} vs {ref}", flush=True)

    sub = adata[adata.obs["state_label"].isin([group, ref])].copy()

    counts = sub.obs["state_label"].value_counts()
    print(counts, flush=True)

    if group not in counts.index or ref not in counts.index:
        print("Skipping missing state", flush=True)
        continue

    if counts[group] < 20 or counts[ref] < 20:
        print("Skipping low cell count", flush=True)
        continue

    sc.tl.rank_genes_groups(
        sub,
        groupby="state_label",
        groups=[group],
        reference=ref,
        method="wilcoxon"
    )

    deg = sc.get.rank_genes_groups_df(sub, group=group)

    safe_name = f"{group}_vs_{ref}"
    out_file = f"{OUTDIR}/{safe_name}_DE.tsv"

    deg.to_csv(out_file, sep="\t", index=False)

    sig_up = deg[
        (deg["pvals_adj"] < 0.05) &
        (deg["logfoldchanges"] > 1)
    ]

    sig_down = deg[
        (deg["pvals_adj"] < 0.05) &
        (deg["logfoldchanges"] < -1)
    ]

    summary_rows.append({
        "comparison": safe_name,
        "group_cells": int(counts[group]),
        "reference_cells": int(counts[ref]),
        "sig_up": len(sig_up),
        "sig_down": len(sig_down),
        "top_up_genes": ",".join(sig_up.head(15)["names"].astype(str)),
        "top_down_genes": ",".join(sig_down.head(15)["names"].astype(str)),
        "output_file": out_file
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    f"{OUTDIR}/state_vs_state_DE_summary.tsv",
    sep="\t",
    index=False
)

print("\nDONE")
print(summary.to_string(index=False), flush=True)
