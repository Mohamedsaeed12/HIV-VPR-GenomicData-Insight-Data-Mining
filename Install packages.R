# Script to Install Required R Packages

# Function to check and install missing packages
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    if (package %in% c("affy", "limma", "clusterProfiler", "org.Hs.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
}

# List of packages to be installed
packages <- c("affy", "GEOquery", "tidyverse", "reshape2", "ggplot2", "ggfortify", 
              "limma", "clusterProfiler", "org.Hs.eg.db", "UpSetR", "gplots", 
              "cluster", "dbscan", "Rtsne")

# Install packages
for (pkg in packages) {
  install_if_missing(pkg)
}

# Load packages
lapply(packages, library, character.only = TRUE)
