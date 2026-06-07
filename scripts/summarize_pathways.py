import pandas as pd

files = {
    "TNBC": "results/pathway_TNBC/Enrichr.human.enrichr.reports.txt",
    "ER": "results/pathway_ER/Enrichr.human.enrichr.reports.txt"
}

rows=[]

for group, f in files.items():
    df = pd.read_csv(f, sep="\t")
    df["group"] = group
    df = df.sort_values("Adjusted P-value")
    rows.append(df)

out = pd.concat(rows, ignore_index=True)

out.to_csv(
    "results/TNBC_vs_ER_pathway_enrichment_summary.tsv",
    sep="\t",
    index=False
)

sig = out[out["Adjusted P-value"] < 0.05]

print("Significant pathways:")
print(sig[["group","Gene_set","Term","Adjusted P-value","Combined Score"]].to_string(index=False))

print("\nTop 20 each group:")
for group in ["TNBC","ER"]:
    print("\n", group)
    print(
        out[out["group"]==group]
        [["Gene_set","Term","Adjusted P-value","Combined Score"]]
        .head(20)
        .to_string(index=False)
    )
