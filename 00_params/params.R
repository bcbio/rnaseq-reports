# info params
date <- "YYYYMMDD"
basedir <- "./" # where to write down output files

# Your data
# This is the file used to run nf-core or compatible to that
coldata_fn <- "/Path/to/metadata/meta.csv"
# This file is inside star_salmon/ folder
counts_fn <- "/path/to/nf-core/output/star_salmon/salmon.merged.gene_counts.tsv"
# This folder called "multiqc_report_data" is inside the output directory star_salmon inside multiqc folder
multiqc_data_dir <- "/path/to/nf-core/output/multiqc/star_salmon/multiqc_report_data"
# This file is inside the genome folder in the output directory, use this only for non-model organism
# gtf_fn='/path/to/nf-core/output/genome/hg38.filtered.gtf'
gtf_fn <- NULL
