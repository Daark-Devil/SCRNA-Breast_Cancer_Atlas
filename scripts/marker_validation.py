import scanpy as sc
import glob
import os

markers = {
    "T_cell": ["CD3D","CD3E","TRBC1"],
    "Macrophage": ["LYZ","C1QA","C1QB"],
    "Epithelial": ["EPCAM","KRT8","KRT18"],
    "Fibroblast": ["COL1A1","COL1A2","DCN"],
    "Endothelial": ["PECAM1","VWF","KDR"],
    "NK_cell": ["NKG7","GNLY","KLRD1"],
    "B_cell_APC": ["MS4A1","CD79A","CD74"],
    "Plasma_cell": ["MZB1","JCHAIN","SDC1"]
}

sample = "GSM4909293_HER2-MH0161"

adata = sc.read_h5ad(
    f"results/{sample}/h5ad/{sample}_annotated.h5ad"
)

sc.pl.dotplot(
    adata,
    markers,
    groupby="cell_type",
    save="_marker_validation.png",
    show=False
)

print("Marker validation plot saved")
