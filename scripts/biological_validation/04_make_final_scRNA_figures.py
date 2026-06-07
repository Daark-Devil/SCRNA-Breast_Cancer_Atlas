import os
import pandas as pd
import matplotlib.pyplot as plt

FIGDIR = "results/biological_validation/figures"
os.makedirs(FIGDIR, exist_ok=True)

summary = pd.read_csv(
    "results/biological_validation/state_comparisons/state_vs_state_DE_summary.tsv",
    sep="\t"
)

summary["total_sig"] = summary["sig_up"] + summary["sig_down"]

plt.figure(figsize=(10,5))
plt.bar(summary["comparison"], summary["total_sig"])
plt.xticks(rotation=90)
plt.ylabel("Significant genes")
plt.title("State-specific epithelial differential expression")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/state_DE_gene_counts.png", dpi=300)
plt.close()

subtype = pd.read_csv(
    "results/epithelial_diagnostics/epi_cluster_subtype_percent.tsv",
    sep="\t",
    index_col=0
)

plt.figure(figsize=(9,10))
plt.imshow(subtype, aspect="auto")
plt.colorbar(label="Percent of cluster")
plt.xticks(range(len(subtype.columns)), subtype.columns, rotation=45, ha="right")
plt.yticks(range(len(subtype.index)), subtype.index)
plt.title("Epithelial cluster subtype composition")
plt.tight_layout()
plt.savefig(f"{FIGDIR}/epithelial_cluster_subtype_heatmap.png", dpi=300)
plt.close()

print("Saved final scRNA figures")
