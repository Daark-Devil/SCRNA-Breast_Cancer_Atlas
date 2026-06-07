import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse

FIGDIR = "results/biological_validation/figures"
os.makedirs(FIGDIR, exist_ok=True)

# -------------------------
# 1. UMAP normal vs cancer
# -------------------------
adata = sc.read_h5ad("results/epithelial_subcluster/epithelial_harmony_subcluster.h5ad")

adata.obs["normal_vs_cancer"] = adata.obs["subtype"].astype(str).map({
    "Normal": "Normal",
    "ER": "Cancer",
    "HER2": "Cancer",
    "TNBC": "Cancer",
    "mER": "Cancer",
    "BRCA1_preneoplastic": "Preneoplastic"
}).fillna("Other")

umap = adata.obsm["X_umap"]

for col, fname, title in [
    ("normal_vs_cancer", "umap_normal_vs_cancer_manual.png", "Normal vs Cancer epithelial cells"),
    ("subtype", "umap_subtype_manual.png", "Breast cancer subtype distribution"),
    ("epi_leiden", "umap_epithelial_states_manual.png", "Epithelial state clusters")
]:
    groups = adata.obs[col].astype(str).unique()

    plt.figure(figsize=(7,5))
    for g in groups:
        idx = adata.obs[col].astype(str) == g
        plt.scatter(
            umap[idx,0],
            umap[idx,1],
            s=2,
            alpha=0.55,
            label=g
        )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(title)
    plt.legend(markerscale=5, bbox_to_anchor=(1.05,1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{FIGDIR}/{fname}", dpi=300)
    plt.close()

# -------------------------
# 2. Stacked bar composition
# -------------------------
comp = pd.crosstab(
    adata.obs["epi_leiden"].astype(str),
    adata.obs["subtype"].astype(str),
    normalize="index"
) * 100

counts = adata.obs["epi_leiden"].astype(str).value_counts()
keep_states = counts[counts >= 50].index.tolist()
comp = comp.loc[keep_states]

comp["dominant_percent"] = comp.max(axis=1)
comp = comp.sort_values("dominant_percent", ascending=False)
comp = comp.drop(columns=["dominant_percent"])

ax = comp.plot(
    kind="bar",
    stacked=True,
    figsize=(13,6),
    width=0.85
)

plt.ylabel("Percent of epithelial state")
plt.xlabel("Epithelial state cluster")
plt.title("Subtype composition of epithelial states")
plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=8)
plt.tight_layout()
plt.savefig(f"{FIGDIR}/epithelial_state_subtype_stacked_bar.png", dpi=300)
plt.close()

# -------------------------
# 3. Manual marker bubble plot
# -------------------------
raw = sc.read_h5ad(
    "results/biological_validation/state_comparisons/state_anchor_log_normalized.h5ad"
)

marker_genes = [
    "EPCAM","KRT8","KRT18","KRT19",
    "GATA3","AGR3","TFF3","XBP1",
    "ERBB2","MIEN1",
    "VIM","KRT17","S100A6","GSTP1","MDK","SOX4",
    "MKI67","TOP2A"
]

genes = [g for g in marker_genes if g in raw.var_names]

states = [
    "Normal_enriched_22",
    "ER_enriched_37",
    "TNBC_enriched_33",
    "TNBC_enriched_9",
    "HER2_enriched_18",
    "HER2_enriched_45"
]

rows = []

for state in states:
    sub = raw[raw.obs["state_label"] == state]

    for gene in genes:
        x = sub[:, gene].X

        if sparse.issparse(x):
            vals = np.asarray(x.todense()).ravel()
        else:
            vals = np.asarray(x).ravel()

        rows.append({
            "state": state,
            "gene": gene,
            "mean_expression": vals.mean(),
            "percent_expressing": (vals > 0).mean() * 100
        })

dot = pd.DataFrame(rows)

x_map = {g:i for i,g in enumerate(genes)}
y_map = {s:i for i,s in enumerate(states)}

plt.figure(figsize=(12,5))

plt.scatter(
    dot["gene"].map(x_map),
    dot["state"].map(y_map),
    s=dot["percent_expressing"] * 4,
    c=dot["mean_expression"],
    alpha=0.85
)

plt.colorbar(label="Mean log-normalized expression")
plt.xticks(range(len(genes)), genes, rotation=45, ha="right")
plt.yticks(range(len(states)), states)
plt.xlabel("Marker gene")
plt.ylabel("Epithelial state")
plt.title("Subtype-associated epithelial marker programs")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/epithelial_state_marker_bubble_plot.png", dpi=300)
plt.close()

print("Saved publication scRNA visuals:")
print(f"{FIGDIR}/umap_normal_vs_cancer_manual.png")
print(f"{FIGDIR}/umap_subtype_manual.png")
print(f"{FIGDIR}/umap_epithelial_states_manual.png")
print(f"{FIGDIR}/epithelial_state_subtype_stacked_bar.png")
print(f"{FIGDIR}/epithelial_state_marker_bubble_plot.png")
