setwd("~/Documents/REPOSITORY/BioInformatics/R_BioInfo")
getwd()

# This script works on Sequence Read Archives from NCBI
biocLite("SRAdb")
library("SRAdb")
sqlfile <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')
sql_file <- getSRAdbFile()
sra_cons <- dbConnect(SQLite(), sql_file)
sra_tabs <- dbListTables(sra_con)
