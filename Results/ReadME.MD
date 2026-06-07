# Results Directory Guide

The `results/` directory contains all outputs generated throughout the analysis workflow, from initial atlas construction to biological validation and tumor-state discovery.

```text
results/
├── integrated_atlas/
├── epithelial_clean/
├── epithelial_subcluster/
├── epithelial_final/
├── epithelial_full_pairwise_DE/
├── biological_validation/
├── marker_validation/
├── pathway_ER/
├── pathway_TNBC/
├── GSM*/
├── tables/
└── testing_and_per_sample_outputs/
```

---

## integrated_atlas/

Contains outputs generated after integration of all breast cancer samples into a unified single-cell atlas.

Key contents:

* Integrated UMAP embeddings
* Leiden clustering results
* Cell-type assignments
* Atlas-level summaries

Purpose:

Provides the global cellular landscape of the complete breast cancer cohort before tumor-specific analyses.

---

## epithelial_clean/

Contains epithelial-cell subsets after removal of low-confidence clusters and non-epithelial contamination.

Key contents:

* Clean epithelial UMAPs
* Refined cluster assignments
* Quality-control summaries

Purpose:

Creates a high-confidence epithelial population for downstream tumor-state discovery.

---

## epithelial_subcluster/

Contains the first round of epithelial re-clustering.

Key contents:

* Epithelial-only UMAPs
* Cluster marker genes
* Subcluster summaries

Purpose:

Identifies transcriptionally distinct epithelial populations within the breast cancer atlas.

---

## epithelial_final/

Contains the final epithelial-state analysis used throughout the study.

Key contents:

* Final epithelial state assignments
* Tumor-state labels
* State-specific marker genes
* Subtype composition summaries

Purpose:

Represents the primary biological result of the project.

Major states identified:

* Normal-enriched
* ER-enriched
* HER2-enriched
* TNBC-enriched

---

## epithelial_full_pairwise_DE/

Contains differential expression results for all epithelial-state comparisons.

Key contents:

* ER vs Normal
* HER2 vs Normal
* TNBC vs Normal
* TNBC vs ER
* HER2 vs ER
* mER vs ER

Files include:

* Differential expression tables
* Fold-change statistics
* Significance values

Purpose:

Identifies genes driving biological differences between epithelial tumor states.

---

## biological_validation/

Contains analyses used to validate biological relevance of identified tumor states.

Structure:

```text
biological_validation/
├── figures/
├── pathway_enrichment/
├── state_comparisons/
└── state_markers/
```

### figures/

Publication-ready visualizations.

Examples:

* Normal vs Cancer UMAP
* Tumor-state UMAP
* Marker bubble plots
* State composition plots

### pathway_enrichment/

Pathway enrichment outputs generated from state-specific marker genes.

Examples:

* Enriched signaling pathways
* Functional annotations
* Cancer-associated biological processes

### state_comparisons/

Differential expression summaries between epithelial states.

Includes:

* Top marker genes
* Upregulated gene lists
* Downregulated gene lists
* Combined state summaries

### state_markers/

Marker genes defining each epithelial state.

Purpose:

Provides biological interpretation and validation of tumor-state assignments.

---

## marker_validation/

Contains validation analyses used during annotation and cluster interpretation.

Examples:

* Marker gene scoring
* Cell-type validation
* Annotation confidence checks

Purpose:

Ensures accurate cell-type identification throughout the atlas.

---

## pathway_ER/

Pathway enrichment analyses focused on ER-associated epithelial states.

Purpose:

Characterizes biological programs enriched in ER-positive disease.

---

## pathway_TNBC/

Pathway enrichment analyses focused on TNBC-associated epithelial states.

Purpose:

Identifies aggressive and subtype-specific TNBC transcriptional programs.

---

## GSM*/

Contains per-sample outputs generated during processing of individual GEO samples.

Examples:

* Sample-specific UMAPs
* QC summaries
* Cluster markers
* Cell-type composition reports

Purpose:

Provides complete sample-level traceability and reproducibility.

---

## tables/

Atlas-level summary tables.

Examples:

* Cell-type composition
* QC statistics
* Sample summaries

Purpose:

Central location for cohort-wide quantitative summaries.

---

## testing_and_per_sample_outputs/

Supporting outputs generated during pipeline development and validation.

Contents include:

* Intermediate analyses
* Pilot runs
* Per-sample diagnostics
* Validation outputs

Purpose:

Provides transparency and reproducibility while keeping primary results separated from supporting analyses.
