import pandas as pd

deg = pd.read_csv(
    "results/TNBC_vs_ER_Epithelial_DEG.tsv",
    sep="\t"
)

deg = deg[
    deg["logfoldchanges"].notna()
]

tnbc = deg[
    (deg["pvals_adj"] < 0.05)
    &
    (deg["logfoldchanges"] > 1)
]

er = deg[
    (deg["pvals_adj"] < 0.05)
    &
    (deg["logfoldchanges"] < -1)
]

tnbc["names"].to_csv(
    "results/TNBC_gene_list.txt",
    index=False,
    header=False
)

er["names"].to_csv(
    "results/ER_gene_list.txt",
    index=False,
    header=False
)

print("TNBC genes:", len(tnbc))
print("ER genes:", len(er))
