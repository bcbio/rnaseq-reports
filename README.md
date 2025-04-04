
Templates with ![](https://img.shields.io/badge/status-stable-green) revision indicates that the components or processes have undergone comprehensive parameterization and testing.

Templates with ![](https://img.shields.io/badge/status-alpha-yellow) revision indicates that the components or processes are currently being tested. There is some test data available, but there are parameters that need to be set up manually within the code.

Templates with ![](https://img.shields.io/badge/status-draft-grey) revision indicates that the components or processes are not fully tested. There is no test data available, parameters need to be set up manually within the code, and specific code changes are required based on the data used.

Read [main page](https://github.com/bcbio) to know how to collaborate with us.

# Guideline for RNAseq downstream analysis

- Set the working directory to this file level. We recommend to use **Projects** in Rstudio.
- Use `install_dependencies.R` to install all packages used in these reports.

## Run data with nf-core rnaseq

This templates assume data has been processed by [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0/docs/usage).
We recommend to use the `samplesheet.csv` used with nf-core as metadata file, where other relevant columns can be there even if they are not used by the pipeline.

## Downstream analysis

- **Modify** `information.R` with the right information. You can use this file with any other Rmd to include the project/analysis information.
- **Modify** the `00_params/params.R` that goes together with the Rmd templates with the right input files location from the output of [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0/docs/output).
- `params*example.R` are parameters pointing to test data to be used as an example to test the reports.
- `run_markdown.R` is an example of code to run the Rmd with different parameters.

### Quality assessment

![](https://img.shields.io/badge/status-stable-green) `01_quality_assessment/QC.Rmd` is a template report that needs `00_params/params.R` for `nf-core/rnaseq` outputs.
 
### Differential expression

![](https://img.shields.io/badge/status-stable-green) `02_differential_expression/DEG.Rmd` is a template for comparison between two groups. It supports - multiple contrasts.
- `00_params/params.R` has the information for the input files to load. You need to point to `nf-core/rnaseq` output files.

On the `YAML` header file of the `Rmd` you can specify some parameters or just set them up in the first chunk of code of the template. This template has examples of:

- sub-setting data
- two groups comparison
- volcano plot
- MA plot
- Pathway analysis: Over-Representation Analysis and Gene-Set-Enrichment Analysis
- Tables

### Comparative analysis

- ![](https://img.shields.io/badge/status-alpha-yellow) `03_comparative/Pair-wise-comparison-analysis.Rmd` shows an example on how to compare two differential expression analysis from the `DEG.Rmd` template.
- ![](https://img.shields.io/badge/status-alpha-yellow)  `03_comparative/Intersections.Rmd` shows an example on how to compare multiple differential expression analyses from `DEG.Rmd` and find intersections.

### Functional analysis

- ![](https://img.shields.io/badge/status-draft-grey) `03_functional/GSVA.Rmd` shows an example on how to use [GSVA package](https://bioconductor.org/packages/release/bioc/html/GSVA.html) for estimating variation of gene set enrichment through the samples of a expression data set
- ![](https://img.shields.io/badge/status-draft-grey)  `03_functional/Nonmodel_Organism_Pathway_Analysis.Rmd` shows an example on how to run Gene Ontology over-representation, KEGG over-representation, and KEGG gene set enrichment analysis (GSEA) for non-model organisms using data from Uniprot. `params_nonmodel_org_pathways.R` has the information for the input files to load.
- ![](https://img.shields.io/badge/status-draft-grey)  `03_functional/Immune-deconvolution.Rmd` shows an example on how to run immune cell type deconvolution. `params_immune_deconv.R` has the information for the input files to load.

### Gene pattern analysis

- ![](https://img.shields.io/badge/status-alpha-yellow) `04_gene_patterns/WGCNA.Rmd` shows an example on how to use the [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html) package to find gene modules in the gene expression data.
- ![](https://img.shields.io/badge/status-alpha-yellow) `04_gene_patterns/DEGpatterns.Rmd` shows an example on hot to cluster set of genes across conditions and time points to identify specific profiles.



