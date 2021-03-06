# PracticeSpace
This file contains the data and `codes` for all the exercises and practice problems

# Necessary Dependencies in R
The following chunk consist of necessary packages required to execute the practice operation
```{r}
source("http://www.bioconductor.org/biocLite.R")
# list of packages

# bioc package
pkgs_bioc <- c("Biostrings", "seqinr", "DESeq2", "ArrayExpress", "Biobase", "annotate", "GEOquery", "phangorn",
               "BiocInstaller")

# general packages
pkgs_gen <- c("gcookbook", "vegan", 
          "tidyverse", "hrbrthemes", "ggplot2", "reshape2",
          "dplyr", "cluster", "rafalib", "RColorBrewer", "lme4", "vegan", "gplots")

pkgs <- c(pkgs_bioc, pkgs_gen)

biocLite(pkgs_bioc, suppressUpdates = T)

# check wether installed or not
install_pkg <- pkgs_gen[!pkgs_gen %in% installed.packages()]

# install required packages one by one
for(lib in install_pkg) install.packages(lib, dependencies = T, quiet = T)

# returns list of packages installed
lapply(pkgs, library, character.only = T)

```
**Alternatively**
use the `callforBio()` package to obtain the required packages for the biocomputation.
**NOTE:** You need to have a table with only seperate columns enlisting your required packages noted during literature survey.
Make sure not to include any header and save the table as `.txt` tab-delim. Provide the `BioC packages` in the first column and `base` packages in the second. Re-arrange the code according to your ease.
Please follow the directions as follow.
```{r}
source("https://raw.githubusercontent.com/arpankbasak/PracticeSpace/master/callForBio.R")
bioc_pkg <- as.character(read.table("./req_pkgs.txt", sep = "\t")[[1]])
oper_pkgs <- as.character(read.table("./req_pkgs.txt", sep = "\t")[[2]])

# Please run/source the function callForBio.R
# bioC argument provides the T/F to call biocLite or not
callForBio(bioc_pkg, oper_pkgs, bioC = F)
```
