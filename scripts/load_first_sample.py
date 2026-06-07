import scanpy as sc
import pandas as pd

sample = "GSM4909281_TN-MH0126"

adata = sc.read_mtx(
    f"GSE161529/{sample}-matrix.mtx.gz"
).T

features = pd.read_csv(
    "GSE161529/GSE161529_features.tsv.gz",
    sep="\t",
    header=None,
    names=["gene_id", "gene_symbol", "feature_type"]
)

barcodes = pd.read_csv(
    f"GSE161529/{sample}-barcodes.tsv.gz",
    header=None,
    names=["barcode"]
)

print("Feature file preview:")
print(features.head())
print("\nBarcode file preview:")
print(barcodes.head())
print("\nRaw AnnData:")
print(adata)

adata.var["gene_id"] = features["gene_id"].astype(str).values
adata.var["gene_symbol"] = features["gene_symbol"].astype(str).values
adata.var["feature_type"] = features["feature_type"].astype(str).values

adata.var_names = adata.var["gene_symbol"].astype(str)
adata.var_names_make_unique()
adata.var.index.name = None

adata.obs_names = sample + "_" + barcodes["barcode"].astype(str)
adata.obs.index.name = None

adata.obs["sample_id"] = sample
adata.obs["subtype"] = "TNBC"
adata.obs["condition"] = "Tumor"

print("\nFinal AnnData:")
print(adata)
print("\nFinal var preview:")
print(adata.var.head())
print("\nFinal obs preview:")
print(adata.obs.head())

adata.write("GSE161529_processed/TN_MH0126_raw.h5ad")

print("\nSaved: GSE161529_processed/TN_MH0126_raw.h5ad")
