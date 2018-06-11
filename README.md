# PracticeSpace
This file contains the data and `codes` for all the exercises and practice problems

# Necessary Dependencies in R
The following chunk consist of necessary packages required to execute the practice operation
```{r}
source("http://www.bioconductor.org/biocLite.R")
# list of packages
pkgs <- c("Biostrings", "seqinr", "DESeq2", "ArrayExpress", "Biobase", "annotate", "GEOquery", "gcookbook", "ape", "edgeR",
              "tidyverse", "hrbrthemes", "ggplot2", "reshape2",
              "dplyr", "cluster", "rafalib", "RColorBrewer", "lme4", "vegan", "gplots", "BiocInstaller")

# check wether installed or not
install_pkg <- pkgs[!pkgs %in% installed.packages()]

# install required packages one by one
for(lib in install_pkg) install.packages(lib, dependencies = T, quiet = T)

# returns list of packages installed
lapply(pkg_plot, library, character.only = T)
```
# Necessary Dependencies in Python
