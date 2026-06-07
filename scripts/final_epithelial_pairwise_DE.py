import os
import scanpy as sc
import pandas as pd

os.makedirs("results/epithelial_final/pairwise_DE", exist_ok=True)

adata = sc.read_h5ad(
    "results/epithelial_final/final_epithelial_states.h5ad"
)

comparisons = [
    ("TNBC", "Normal"),
    ("ER", "Normal"),
    ("HER2", "Normal"),
    ("mER", "ER"),
    ("TNBC", "ER"),
    ("HER2", "ER")
]

summary_rows = []

for group, ref in comparisons:

    print(f"\nRunning DE: {group} vs {ref}", flush=True)

    sub = adata[
        adata.obs["subtype"].isin([group, ref])
    ].copy()

    print(sub.obs["subtype"].value_counts(), flush=True)

    if sub.obs["subtype"].nunique() < 2:
        print("Skipping: missing one group", flush=True)
        continue

    sc.tl.rank_genes_groups(
        sub,
        groupby="subtype",
        groups=[group],
        reference=ref,
        method="wilcoxon"
    )

    deg = sc.get.rank_genes_groups_df(
        sub,
        group=group
    )

    out_file = f"results/epithelial_final/pairwise_DE/{group}_vs_{ref}_DE.tsv"

    deg.to_csv(
        out_file,
        sep="\t",
        index=False
    )

    sig_up = deg[
        (deg["pvals_adj"] < 0.05) &
        (deg["logfoldchanges"] > 1)
    ]

    sig_down = deg[
        (deg["pvals_adj"] < 0.05) &
        (deg["logfoldchanges"] < -1)
    ]

    summary_rows.append({
        "comparison": f"{group}_vs_{ref}",
        "cells_total": sub.n_obs,
        "group_cells": (sub.obs["subtype"] == group).sum(),
        "reference_cells": (sub.obs["subtype"] == ref).sum(),
        "significant_up": len(sig_up),
        "significant_down": len(sig_down),
        "output_file": out_file
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    "results/epithelial_final/pairwise_DE/pairwise_DE_summary.tsv",
    sep="\t",
    index=False
)

print("\nDONE")
print(summary.to_string(index=False))
