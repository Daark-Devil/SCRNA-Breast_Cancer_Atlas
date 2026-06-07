import pandas as pd

markers = pd.read_csv(
    "results/epithelial_clean/markers/clean_epithelial_cluster_markers_top10.tsv",
    sep="\t"
)

immune = {
    "IGLC2","IGHA1","IGHM","IGKC","CD3D","CD3E","TRAC","NKG7",
    "LYZ","C1QA","C1QB","CTSS","HLA-DQA1","HLA-DQB1","HLA-DRB5",
    "LTB","CCL3L1","TRGC1"
}

stromal = {
    "COL1A1","COL1A2","DCN","LUM","VCAN","FBLN1","AEBP1",
    "POSTN","FN1","MMP3","THBS2","COL5A2"
}

endothelial = {
    "PECAM1","VWF","KDR","ENG","PLPP3","SPARCL1","RGCC"
}

rows=[]

for cl, sub in markers.groupby("group"):
    genes = set(sub["names"].astype(str))

    flags=[]
    if genes & immune:
        flags.append("immune")
    if genes & stromal:
        flags.append("stromal")
    if genes & endothelial:
        flags.append("endothelial")

    rows.append({
        "cluster": cl,
        "flags": ",".join(flags) if flags else "clean",
        "top_genes": ",".join(sub["names"].astype(str).head(10))
    })

out = pd.DataFrame(rows)

out.to_csv(
    "results/epithelial_clean/clean_epithelial_cluster_flags.tsv",
    sep="\t",
    index=False
)

print(out.to_string(index=False))
