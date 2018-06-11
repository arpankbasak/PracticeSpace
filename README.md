# PracticeSpace
This file contains the data and `codes` for all the exercises and practice problems

# Necessary Dependencies in R
The following chunk consist of necessary packages required to execute the practice operation
```{r}
pkg_plot <- c("Biostrings", "seqinr", "DESeq2", "ArrayExpress", "Biobase", "annotate", "GEOquery", "gcookbook", "ape", "edgeR",
              "tidyverse", "hrbrthemes", "ggplot2", "reshape2",
              "dplyr", "cluster", "rafalib", "RColorBrewer", "lme4", "vegan", "gplots", "BiocInstaller")


install_pkg <- pkg_plot[!pkg_plot %in% installed.packages()]
for(lib in install_pkg) install.packages(lib, dependencies = T, quiet = T)
lapply(pkg_plot, library, character.only = T)
```
# Necessary Dependencies in Python
