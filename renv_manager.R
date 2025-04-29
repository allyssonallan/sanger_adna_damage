#!/usr/bin/env Rscript
# renv_manager.R: initialize or restore R environment using renv

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Ensure renv is installed
install_if_missing("renv")

# Ensure BiocManager is installed for Bioconductor packages
install_if_missing("BiocManager")

if (file.exists("renv.lock")) {
  message("Restoring R dependencies from renv.lock...")
  renv::restore()
} else {
  message("Initializing renv and installing dependencies based on project code...")
  renv::init()
}