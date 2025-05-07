
# DEG Visualization Dashboards

This repository provides **two interactive R Shiny dashboards** to visualize and explore differential gene expression (DEG) analysis results:

- **`visu_DEG.R`** – *Pathway-Centric View*
- **`visu_DEG2.R`** – *Gene-Centric View*

These tools are designed to **complement each other**, not replace. Together, they provide two views to query DEG outputs from 1) pathway-level enrichment; 2) individual gene-level behavior

---

## What Do These two Dashboards Do?

| Feature                              | `visu_DEG.R` (Pathway-Centric)                    | `visu_DEG2.R` (Gene-Centric)                   |
|--------------------------------------|--------------------------------------------------|-----------------------------------------------|
| Focus                                | Enriched gene sets/pathways                      | Individual genes                              |
| Visualizations                       | Volcano plot, boxplot, heatmap                   | Volcano plot, boxplot                         |
| Interactive gene/pathway selection   | Enriched pathway dropdown                        | Gene selection                                |
| Output Tables                        | DEG output for significant genes in the selected pathway | Enriched pathways associated with selected gene |
| Expression Data Integration          | Heatmap + Boxplots for gene expression           | Boxplots only                                 |

---

## Inputs

Both dashboards use data from [bcbioRNASeq test data](https://github.com/bcbio/bcbioR-test-data/tree/main/rnaseq/DEG_visualization), including template outputs from DEG.rmd report:

- `full_sample_type_normal_vs_tumor_pathways.csv`: Annotated pathways and genes.
- `full_sample_type_normal_vs_tumor_deg.csv`: DEG table with gene names, logFC, adjusted p-values.
- `full_expression.csv`: Normalized expression matrix (genes x samples).
- `deseq_coldata.rds`: Sample metadata (including condition labels).

Where  `sample_type` is the column name in meta-data specifying the group label for each sample. `normal_vs_tumor` is the constrast we used for DEG analysis.

> **Note:** These are publicly hosted test datasets and will be downloaded directly from GitHub using URLs.

---

##  How to Run

### Prerequisites

Install required R packages:

```R
pkgs_needed <- c("shiny", "shinydashboard", "DT", "plotly", "ComplexHeatmap",
"ggpubr", "grid", "purrr", "dplyr", "glue", "data.table", "magrittr", "viridis","qs", "R.utils")

# Identify which packages are not installed
pkgs_to_install <- pkgs_needed[!pkgs_needed %in% rownames(installed.packages())]

# Install missing packages
if (length(pkgs_to_install) > 0) {
  install.packages(pkgs_to_install)
}
```

### Launch

In your R console or RStudio:

```R
shiny::runApp("visu_DEG.R")   # for pathway-centric view
shiny::runApp("visu_DEG2.R")  # for gene-centric view
```

---

## Features Summary

### `visu_DEG.R` – Pathway View
- Select a **gene set** (pathway) from a dropdown menu.
- View:
  - A **table of associated genes** in the selected pathway, only genes significantly contributing to the pathway enrichment.
  - **Volcano plot** highlighting selected genes.
  - **Expression boxplots** for selected gene, another dropdown menu to select each gene.
  - **Heatmap** for the selected gene by samples.

### `visu_DEG2.R` – Gene View
- Select a **gene** from dropdown.
- View:
  - Table of **associated pathways**.
  - **Volcano plot** highlighting selected gene.
  - **Expression boxplot** for selected gene.

---

## Notes

- These apps are useful in validating key DEG findings and exploring biological hypotheses interactively.
- Designed for our results from DEG.rmd template report
- For local use, modify the URLs to point to your local files if needed.

---

