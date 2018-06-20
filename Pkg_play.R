setwd("~/Documents/REPOSITORY/MCB_KY_DataRepository/job_align_2018-06-16")
getwd()

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("GEOquery")
biocLite("ArrayExpress")
biocLite("biomaRt")

pkgs <- c("DESeq2", "ArrayExpress", "Biobase", "annotate", "GEOquery", "gcookbook", "ape", "edgeR", "ggsignif", "igraph",
              "tidyverse", "hrbrthemes", "ggplot2", "reshape2",
              "dplyr", "cluster", "rafalib", "RColorBrewer", "lme4", "vegan", "gplots", "BiocInstaller")

install_pkg <- pkgs[!pkgs %in% installed.packages()]
for(lib in install_pkg) install.packages(lib, dependencies = T, quiet = T)
lapply(pkgs, library, character.only = T)

# ===== Using Biomarts ======
library(biomaRt)

# Listing from Phytozome
listMarts()
m <-  useMart('phytozome_mart', 
              host="phytozome.jgi.doe.gov", 
              path ="/biomart/martservice/", 
              dataset = 'phytozome')

listAttributes(m)
obj_gen <- getBM(attributes = c("organism_name", "chr_name1", "gene_name","exon_chrom_start",
                                "exon_chrom_end", "gene_exon", "ortholog_group"), mart = m)

phytozome_db <- useMart
# ===== Using BioStrings ======
library(msa)
library(seqinr)

obj <- readAAStringSetStringSet("q1_1kp_CYP79B2.fasta")
obj_msa <- msa(obj, method = "Muscle", cluster = "upgma", type = "aa")
print(obj_msa, show = "complete")

obj_align <- msaConvert(obj_msa, type = "seqinr::alignment")
obj_dist <- dist.alignment(obj_align)

# Plotting dendogram
# obj_phy <- (obj_dist)

# ===== Using UniProts.ws =====
biocLite("UniProt.ws")
library("UniProt.ws")

availableUniprotSpecies("CYP79B2", keytype = "ENTRY NAME")
# List all possible keys of type entrez gene id
if(interactive()){
  egs = keys(up, "ENTREZ_GENE")
}
res <- select(up,
              keys = c("", ""),
              columns = c("PDB", "UNIGENE", "SEQUENCE"),
              keytype = "ENTREZ_GENE")


# ===== Using seqinr =====
library(seqinr)
