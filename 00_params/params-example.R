# project params
date <- "YYYYMMDD"
basedir <- "./" # where to write down output files

# Example data
coldata_fn <- "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/v3_18/coldata.csv"
counts_fn <- url("https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/v3_18/star_salmon/salmon.merged.gene_counts.tsv")
# This folder is in the output directory inside multiqc folder
multiqc_data_dir <- "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/v3_18/multiqc/star_salmon/multiqc_report_data/"
# This file is inside the genome folder in the output directory
gtf_fn <- "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/refs/heads/main/nf-core/v3_18/genome/genome.filtered.gtf.gz"
se_object <- NA
