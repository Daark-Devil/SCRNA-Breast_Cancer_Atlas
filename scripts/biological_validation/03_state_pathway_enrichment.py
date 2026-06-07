import os
import glob
import gseapy as gp

INDIR = "results/biological_validation/state_comparisons/final_summaries"
OUTDIR = "results/biological_validation/pathway_enrichment"

os.makedirs(OUTDIR, exist_ok=True)

gene_lists = glob.glob(f"{INDIR}/*_gene_list.txt")

gene_sets = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
    "Reactome_2022"
]

failed = []

for glist in gene_lists:
    name = os.path.basename(glist).replace("_gene_list.txt", "")
    out = f"{OUTDIR}/{name}"

    print("Running enrichment:", name, flush=True)

    try:
        gp.enrichr(
            gene_list=glist,
            gene_sets=gene_sets,
            organism="human",
            outdir=out,
            cutoff=1.0
        )
    except Exception as e:
        print("FAILED:", name, str(e), flush=True)
        failed.append([name, str(e)])

with open(f"{OUTDIR}/failed_enrichment_jobs.tsv", "w") as f:
    f.write("gene_list\terror\n")
    for x in failed:
        f.write(x[0] + "\t" + x[1].replace("\t"," ") + "\n")

print("DONE pathway enrichment")
print("Failed:", len(failed))
