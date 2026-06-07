import pandas as pd
from scipy.stats import kruskal

comp = pd.read_csv(
    "results/MASTER_celltype_composition_69.tsv",
    sep="\t"
)

rows=[]

for celltype in sorted(comp["cell_type"].unique()):

    tmp = comp[comp["cell_type"] == celltype]

    groups=[]

    for subtype in sorted(tmp["subtype"].dropna().unique()):

        vals = tmp[tmp["subtype"] == subtype]["percent"]

        if len(vals) > 1:
            groups.append(vals)

    if len(groups) >= 2:

        stat,p = kruskal(*groups)

        rows.append({
            "cell_type": celltype,
            "kruskal_stat": stat,
            "pvalue": p
        })

out = pd.DataFrame(rows)

out = out.sort_values("pvalue")

out.to_csv(
    "results/celltype_statistics.tsv",
    sep="\t",
    index=False
)

print(out)
