import gseapy as gp

gp.enrichr(
    gene_list="results/TNBC_gene_list_clean.txt",
    gene_sets=[
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "Reactome_2022"
    ],
    organism="human",
    outdir="results/pathway_TNBC",
    cutoff=0.05
)

gp.enrichr(
    gene_list="results/ER_gene_list_clean.txt",
    gene_sets=[
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "Reactome_2022"
    ],
    organism="human",
    outdir="results/pathway_ER",
    cutoff=0.05
)

print("Pathway enrichment complete")
