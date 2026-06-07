import os
import glob
import pandas as pd
import scanpy as sc

marker_sets = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "IL7R"],
    "Activated_Tcell": ["CXCL13", "TIGIT", "PDCD1", "CTLA4", "IFNG", "GZMB"],
    "NK_cell": ["NKG7", "GNLY", "PRF1", "KLRD1", "GZMB"],
    "Macrophage": ["LYZ", "C1QA", "C1QB", "CST3", "FCER1G", "AIF1"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "ENG", "EGFL7"],
    "B_cell_APC": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA"],
    "Plasma_cell": ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1"]
}

rows = []
h5ads = sorted(glob.glob("results/*/h5ad/*_annotated.h5ad"))

print("Found h5ad files:", len(h5ads))

for h5 in h5ads:
    sample = h5.split("/")[1]
    print("Processing:", sample)

    outdir = f"results/marker_validation/{sample}"
    os.makedirs(outdir, exist_ok=True)
    sc.settings.figdir = outdir

    adata = sc.read_h5ad(h5)

    if "Tumor_Epithelial" in adata.obs["cell_type"].astype(str).unique():
        adata.obs["cell_type"] = adata.obs["cell_type"].replace({"Tumor_Epithelial": "Epithelial"})

    present_marker_sets = {}
    for ct, genes in marker_sets.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) > 0:
            present_marker_sets[ct] = present

    sc.pl.dotplot(
        adata,
        present_marker_sets,
        groupby="cell_type",
        save=f"_{sample}_marker_dotplot.png",
        show=False
    )

    for ct, genes in present_marker_sets.items():
        if ct not in adata.obs["cell_type"].astype(str).unique():
            continue

        cells = adata[adata.obs["cell_type"].astype(str) == ct]
        n_cells = cells.n_obs

        gene_scores = {}
        for gene in genes:
            x = cells[:, gene].X
            mean_expr = float(x.mean())
            gene_scores[gene] = mean_expr

        marker_score = sum(gene_scores.values()) / len(gene_scores)
        top_gene = max(gene_scores, key=gene_scores.get)

        status = "PASS"
        if n_cells < 30:
            status = "LOW_CELL_COUNT"
        elif marker_score < 0.25:
            status = "REVIEW_LOW_MARKER_SCORE"

        rows.append({
            "sample": sample,
            "cell_type": ct,
            "n_cells": n_cells,
            "marker_score": marker_score,
            "top_marker_gene": top_gene,
            "top_marker_expression": gene_scores[top_gene],
            "status": status
        })

validation = pd.DataFrame(rows)
validation.to_csv("results/MASTER_marker_validation.tsv", sep="\t", index=False)

print("Saved: results/MASTER_marker_validation.tsv")
print(validation["status"].value_counts())
