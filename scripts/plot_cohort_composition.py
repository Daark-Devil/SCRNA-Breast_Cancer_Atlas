import pandas as pd
import matplotlib.pyplot as plt

comp = pd.read_csv("results/MASTER_celltype_composition_69.tsv", sep="\t")
comp["cell_type"] = comp["cell_type"].replace({"Tumor_Epithelial": "Epithelial"})

mean_comp = (
    comp.groupby(["subtype", "cell_type"])["percent"]
    .mean()
    .reset_index()
)

pivot = mean_comp.pivot(index="subtype", columns="cell_type", values="percent").fillna(0)

pivot.to_csv("results/MASTER_mean_celltype_by_subtype.tsv", sep="\t")

plt.figure(figsize=(12,6))
pivot.plot(kind="bar", stacked=True, figsize=(12,6))
plt.ylabel("Mean cell type percentage")
plt.xlabel("Subtype")
plt.title("Mean cell-type composition by breast cancer subtype")
plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
plt.tight_layout()
plt.savefig("results/final_figures/celltype_composition_by_subtype.png", dpi=300)

plt.figure(figsize=(10,6))
plt.imshow(pivot, aspect="auto")
plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=90)
plt.yticks(range(len(pivot.index)), pivot.index)
plt.colorbar(label="Mean %")
plt.title("Cell-type composition heatmap")
plt.tight_layout()
plt.savefig("results/final_figures/celltype_composition_heatmap.png", dpi=300)

print(pivot)
print("Saved cohort composition plots")
