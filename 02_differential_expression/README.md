
# DEG Visualization Dashboards

This repository provides **two interactive R Shiny dashboards** to visualize and explore differential gene expression (DEG) analysis results:

- **`visu_DEG.R`** – *Pathway-Centric View*
- **`visu_DEG2.R`** – *Gene-Centric View*

These tools are designed to **complement each other**, not replace. Together, they provide insights into both global pathway-level enrichment and individual gene-level behavior.

---

## What Do These two Dashboards Do?

| Feature                              | `visu_DEG.R` (Pathway-Centric)                    | `visu_DEG2.R` (Gene-Centric)                   |
|--------------------------------------|--------------------------------------------------|-----------------------------------------------|
| Focus                                | Enriched gene sets/pathways                      | Individual genes                              |
| Key Input                            | Annotated DEG pathway CSV                        | Annotated DEG pathway CSV                     |
| Visualizations                       | Volcano plot, boxplot, heatmap                   | Volcano plot, boxplot                         |
| Interactive gene/pathway selection   | Gene set dropdown, gene selection                | Gene selection only                           |
| Output Tables                        | DEG table for selected gene set                  | Pathways associated with selected gene        |
| Expression Data Integration          | Heatmap + Boxplots for gene expression           | Boxplots only                                 |

---

## Inputs

Both dashboards use data from [bcbioRNASeq test data](https://github.com/bcbio/bcbioR-test-data/tree/main/rnaseq/DEG_visualization), including:

- `full_sample_type_normal_vs_tumor_pathways.csv`: Annotated pathways and genes.
- `full_sample_type_normal_vs_tumor_deg.csv`: DEG table with gene names, logFC, adjusted p-values.
- `full_expression.csv`: Normalized expression matrix (genes x samples).
- `deseq_coldata.rds`: Sample metadata (including condition labels).

> **Note:** These are publicly hosted test datasets and will be downloaded directly from GitHub using URLs.

---

##  How to Run

### Prerequisites

Install required R packages:

```R
install.packages(c("shiny", "shinydashboard", "DT", "plotly", "ComplexHeatmap", 
                   "ggpubr", "grid", "purrr", "dplyr", "glue", "data.table", "magrittr", "viridis"))
```

You also need `qs` and `R.utils` if using `qs` or URL-based file loading:

```R
install.packages(c("qs", "R.utils"))
```

### Launch

In your R console or RStudio:

```R
shiny::runApp("visu_DEG.R")   # for pathway-centric view
shiny::runApp("visu_DEG2.R")  # for gene-centric view
```

---

## 📊 Features Summary

### `visu_DEG.R` – Pathway View
- Select a **gene set** (pathway) from dropdown.
- View:
  - A **table of genes** in the selected set.
  - **Volcano plot** highlighting selected genes.
  - **Expression boxplots** for selected gene.
  - **Heatmap** for the full gene set.

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

