# Figures Directory Guide

The repository contains both publication figures and supporting visualizations.

```text
main_figures/
supporting_figures/
```

---

# main_figures/

Contains the primary figures required to understand the major findings of the study.

These figures summarize the biological conclusions and should be viewed first.

---

## umap_normal_vs_cancer_manual.png

Normal versus cancer epithelial-cell projection.

Shows clear transcriptomic separation between healthy and malignant epithelial populations.

Key conclusion:

Tumor cells exhibit distinct transcriptional programs relative to normal breast epithelium.

---

## umap_subtype_manual.png

Breast cancer subtype visualization.

Displays:

* Normal
* ER+
* HER2+
* TNBC
* Metastatic ER

Key conclusion:

Different breast cancer subtypes occupy distinct transcriptional regions within the integrated atlas.

---

## epithelial_state_marker_bubble_plot.png

Marker-gene expression across epithelial states.

Visualization:

* Bubble size = percentage of expressing cells
* Bubble color = average expression

Key conclusion:

Each epithelial state possesses a unique molecular signature.

Examples:

ER markers:

* GATA3
* AGR3
* CRIP1

HER2 markers:

* MIEN1
* XBP1

TNBC markers:

* SOX4
* STMN1
* VCAN
* MDK

---

## epithelial_state_subtype_stacked_bar.png

Subtype composition across epithelial states.

Key conclusion:

Different tumor states are enriched for distinct breast cancer subtypes.

---

## state_DE_gene_counts.png

Differential expression summary.

Displays the number of significant genes identified between epithelial-state comparisons.

Key conclusion:

Substantial transcriptional differences exist between tumor states.

---

# supporting_figures/

Contains supplementary analyses used for validation and interpretation.

These figures support the conclusions presented in the main figures.

Categories include:

## Atlas Construction

* Integrated UMAPs
* Clustering diagnostics
* Annotation validation

## Cell-Type Composition

* Cell-type abundance plots
* Subtype comparisons
* Immune composition analyses

## Epithelial Analysis

* Epithelial-only UMAPs
* Subcluster visualizations
* Marker validation figures

## Quality Control

* QC summaries
* Filtering diagnostics
* Validation plots

## Pathway Validation

* Pathway enrichment outputs
* Functional validation analyses

Together these figures provide complete transparency and reproducibility for all analyses performed in the project.
