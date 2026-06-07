import pandas as pd
import matplotlib.pyplot as plt

comp = pd.read_csv("results/MASTER_celltype_composition_69.tsv", sep="\t")
comp["cell_type"] = comp["cell_type"].replace({"Tumor_Epithelial": "Epithelial"})

celltypes = ["Epithelial", "T_cell", "Activated_Tcell", "NK_cell", "Macrophage", "Fibroblast", "Endothelial", "Plasma_cell", "B_cell_APC"]

for ct in celltypes:
    df = comp[comp["cell_type"] == ct].copy()
    if df.empty:
        continue

    groups = []
    labels = []

    for subtype in sorted(df["subtype"].dropna().unique()):
        values = df[df["subtype"] == subtype]["percent"].values
        if len(values) > 0:
            groups.append(values)
            labels.append(subtype)

    plt.figure(figsize=(8,5))
    plt.boxplot(groups, labels=labels)
    plt.ylabel(f"{ct} percentage")
    plt.xlabel("Subtype")
    plt.title(f"{ct} abundance across subtypes")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"results/final_figures/boxplot_{ct}_by_subtype.png", dpi=300)
    plt.close()

print("Saved cell-type boxplots")
