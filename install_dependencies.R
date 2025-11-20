install.packages("BiocManager")
BiocManager::install("renv")
BiocManager::install("remotes")
BiocManager::install("omnideconv/immunedeconv")
BiocManager::install(renv::dependencies(path = ".")[["Package"]])

create_snapshot <- function(folder, extra_libs = NULL) {
  pkgs <- renv::dependencies(path = folder)[["Package"]]
  if (!is.null(extra_libs)) {
    pkgs <- c(pkgs, renv::dependencies(path = extra_libs)[["Package"]])
  }
  renv::snapshot(lockfile = file.path(folder, "renv.lock"), packages = pkgs)
}

# Snapshots for each folder
create_snapshot("01_quality_assessment", extra_libs = "00_libs")
create_snapshot("02_differential_expression", extra_libs = "00_libs")
create_snapshot("03_comparative")
create_snapshot("03_functional")
create_snapshot("04_gene_patterns", extra_libs = "00_libs")
