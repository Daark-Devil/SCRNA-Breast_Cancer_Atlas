import os
import pandas as pd

INDIR = "results/biological_validation/state_comparisons"
OUTDIR = "results/biological_validation/state_comparisons/final_summaries"

os.makedirs(OUTDIR, exist_ok=True)

summary = pd.read_csv(
    f"{INDIR}/state_vs_state_DE_summary.tsv",
    sep="\t"
)

all_rows = []

for _, row in summary.iterrows():

    comp = row["comparison"]
    f = row["output_file"]

    deg = pd.read_csv(f, sep="\t")

    up = deg[
        (deg["pvals_adj"] < 0.05) &
        (deg["logfoldchanges"] > 1)
    ].copy()

    down = deg[
        (deg["pvals_adj"] < 0.05) &
        (deg["logfoldchanges"] < -1)
    ].copy()

    up["comparison"] = comp
    up["direction"] = "up"

    down["comparison"] = comp
    down["direction"] = "down"

    up.head(100).to_csv(
        f"{OUTDIR}/{comp}_top100_up.tsv",
        sep="\t",
        index=False
    )

    down.head(100).to_csv(
        f"{OUTDIR}/{comp}_top100_down.tsv",
        sep="\t",
        index=False
    )

    up["names"].head(300).to_csv(
        f"{OUTDIR}/{comp}_up_gene_list.txt",
        index=False,
        header=False
    )

    down["names"].head(300).to_csv(
        f"{OUTDIR}/{comp}_down_gene_list.txt",
        index=False,
        header=False
    )

    all_rows.append(up.head(50))
    all_rows.append(down.head(50))

combined = pd.concat(all_rows, ignore_index=True)

combined.to_csv(
    f"{OUTDIR}/combined_top_state_DE_genes.tsv",
    sep="\t",
    index=False
)

print("Saved final state-DE summaries")
print(summary.to_string(index=False))
