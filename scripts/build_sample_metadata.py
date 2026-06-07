import os
import re
import pandas as pd

raw_dir = "GSE161529"
out_file = "metadata/sample_metadata.tsv"

samples = []

for fn in sorted(os.listdir(raw_dir)):
    if not fn.endswith("-matrix.mtx.gz"):
        continue

    sample_id = fn.replace("-matrix.mtx.gz", "")
    gsm = sample_id.split("_", 1)[0]
    name = sample_id.split("_", 1)[1] if "_" in sample_id else sample_id

    if name.startswith("TN-"):
        subtype = "TNBC"
        condition = "Tumor"
    elif name.startswith("HER2-"):
        subtype = "HER2"
        condition = "Tumor"
    elif name.startswith("ER-"):
        subtype = "ER"
        condition = "Tumor"
    elif name.startswith("mER-"):
        subtype = "mER"
        condition = "Metastatic_ER"
    elif name.startswith("B1-"):
        subtype = "BRCA1_preneoplastic"
        condition = "Pre_neoplastic"
    elif name.startswith("N-"):
        subtype = "Normal"
        condition = "Normal"
    else:
        subtype = "Unknown"
        condition = "Unknown"

    if "-LN" in name:
        tissue_fraction = "Lymph_node"
    elif "-Epi" in name:
        tissue_fraction = "Epithelial_enriched"
    elif "-Total" in name:
        tissue_fraction = "Total_cells"
    elif "-T" in name or condition == "Tumor":
        tissue_fraction = "Primary_tumor"
    else:
        tissue_fraction = "Unknown"

    patient = name
    patient = re.sub(r"^(TN-B1-|TN-|HER2-|ER-|mER-|B1-|N-)", "", patient)
    patient = re.sub(r"-(Total|Epi|T|LN|T2|T3)$", "", patient)

    samples.append({
        "gsm": gsm,
        "sample_id": sample_id,
        "sample_name": name,
        "patient_id": patient,
        "condition": condition,
        "subtype": subtype,
        "tissue_fraction": tissue_fraction,
        "matrix_file": os.path.join(raw_dir, sample_id + "-matrix.mtx.gz"),
        "barcode_file": os.path.join(raw_dir, sample_id + "-barcodes.tsv.gz")
    })

df = pd.DataFrame(samples)
df.to_csv(out_file, sep="\t", index=False)

print(df.head())
print("\nTotal samples:", len(df))
print("\nSubtype counts:")
print(df["subtype"].value_counts())
print("\nCondition counts:")
print(df["condition"].value_counts())
print("\nTissue fraction counts:")
print(df["tissue_fraction"].value_counts())
print("\nSaved:", out_file)
