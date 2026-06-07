import glob
import scanpy as sc
import pandas as pd
import anndata as ad

files = glob.glob("results/*/h5ad/*_annotated.h5ad")

adatas=[]

for f in files:

    a = sc.read_h5ad(f)

    keep = (
        a.obs["subtype"].isin(["TNBC","ER"])
        &
        (a.obs["cell_type"]=="Epithelial")
    )

    a = a[keep].copy()

    if a.n_obs > 0:
        adatas.append(a)

adata = ad.concat(adatas, join="outer", fill_value=0)

adata.obs["subtype"] = adata.obs["subtype"].astype(str)
adata.obs["subtype"] = adata.obs["subtype"].astype("category")

print(adata.obs["subtype"].value_counts())

sc.tl.rank_genes_groups(
    adata,
    groupby="subtype",
    groups=["TNBC"],
    reference="ER",
    method="wilcoxon"
)

deg = sc.get.rank_genes_groups_df(
    adata,
    group="TNBC"
)

deg.to_csv(
    "results/TNBC_vs_ER_Epithelial_DEG.tsv",
    sep="\t",
    index=False
)

print(deg.head(20))
print("Saved DEG table")
