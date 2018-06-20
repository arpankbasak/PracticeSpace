setwd("~/Documents/REPOSITORY/BioInformatics/R_BioInfo")
getwd()

# This program aims to search for homologous genes between two species

# ====== Loading dependencies ======
source("http://www.bioconductor.org/biocLite.R")

pkgs <- c("Biostrings", "seqinr", "DESeq2", "ArrayExpress", "Biobase", "annotate", "GEOquery", "gcookbook", "ape", "edgeR", "ggsignif", "igraph",
              "tidyverse", "hrbrthemes", "ggplot2", "reshape2",
              "dplyr", "cluster", "rafalib", "RColorBrewer", "lme4", "vegan", "gplots", "BiocInstaller")


install_pkg <- pkgs[!pkgs %in% installed.packages()]
for(lib in install_pkg) install.packages(lib, dependencies = T, quiet = T)
lapply(pkgs, library, character.only = T)

# ======= Quering from NCBI =====
choosebank("ensplants")
qr <- query("CYP83B1", "SP=Arabidopsis Thaliana")
qseq_1 <- getSequence(qr$req[[1]])
qseq_2 <- getSequence(qr$req[[2]])

write.fasta(names = "Seq-1", sequences = qseq_1, file.out = "seq_CYP83B1_AT.fasta")

# ====== Exploring the query =====
length(qseq_1)
table(qseq_1)
GC(qseq_1)


# ====== Custom made query ======

getncbiseq <- function(accession) {
  require("seqinr")
  dbs <- c("genbank", "refseq", "ensplants", "embl", "refseqViruses")
  num_db <- length(dbs)
  for( i in 1:num_db) {
    db <- dbs[i]
    choosebank(db)
    r_query <- try(query(".tmpquery", paste("AC=", accession)), silent = T)
    if(!(inherits(r_query, "try-error")))
    {
      queryname <- "query2"
      p_query <- paste("AC=", accession, sep = "")
      query(`queryname`, `p_query`)
      seq <- getSequence(query2$req[[1]])
      closebank()
      return(seq)
    }
    closebank()
  }
  print(paste("ERROR: accession", accession, "not found"))
}

test <- getncbiseq("NC_001477")

# ====== Pairwise Sequence alignment =====
choosebank("swissprot")
## Query 1

q1 <- query("brugia", "AC=A8PZ80")
qseq <- getSequence(q1$req[[1]])

## Query 2
q2 <- query("ulcerans", "AC=A0PQ23")
rseq <- getSequence(q2$req[[1]])

closebank()

dotPlot(qseq, rseq)

## Generate the sigma matrix for the score generation "Needleman-Wunsch Algorithm"
sig <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = T)
sig

## Pairwise alignment program
s1 <- "AAAATTTGGCC"
s2 <- "TTTAGGAAACC"
ga_seq <- pairwiseAlignment(s1,s2, substitutionMatrix = sig, 
                            gapOpening = -2, gapExtension = -8, scoreOnly = F)

## Peptide seq
data("BLOSUM50")
BLOSUM50
ga_seq <- pairwiseAlignment(c2s(qseq), c2s(rseq), substitutionMatrix = "BLOSUM50", 
                            gapOpening = -2, gapExtension = -8, scoreOnly = F)
ga_seq

printPairwiseAlignment(ga_seq, 60)

## Local Alignment
la_seq <- pairwiseAlignment(c2s(qseq), c2s(rseq), substitutionMatrix = "BLOSUM50", 
                            gapOpening = -2, gapExtension = -8, scoreOnly = F,
                            type = "local")
la_seq
ga_seq

## Long pairwise alignment
# ===== XXXXXX ======
# longPairAlign <- function(alignment, chunksize = 60, returnlist = F)
# {
#   require(Biostrings)
#   
#   #Get the alignment sets
#   seq1 <- pattern(alignment)
#   seq2 <- subject(alignment)
#   ln <- nchar(seq1)
#   #generate sequence based on number of characters in pattern
#   starts <- seq(1, ln, by = chunksize)
#   n <- length(starts)
#   res1 <- 0
#   res2 <- 0
#   for(i in 1:n)
#   {
#     chunk_seq1 <- substring(seq1, starts[i], starts[i] + chunksize - 1)
#     chunk_seq2 <- substring(seq2, starts[i], starts[i] + chunksize - 1)
#     
#     gaps1 <- countPattern("-", chunk_seq1)
#     gaps2 <- countPattern("-", chunk_seq2)
#     
#     res1 <- res1 + chunksize - gaps1
#     res2 <- res2 + chunksize - gaps2
#     
#     if(returnlist == 'F')
#     {
#       print(paste(chunk_seq1, res1))
#       print(paste(chunk_seq2, res2))
#       print(paste(" "))
#     }
#   }
#   if(returnlist == 'T')
#   {
#     vec1 <- s2c(substring(seq1, 1, nchar(seq1)))
#     vec2 <- s2c(substring(seq2, 1, nchar(seq2)))
#     ls <- list(vec1, vec2)
#     return(ls)
#   }
# }
# ====== Pairwise Alignment based on desired chunk =====
printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
{
  require(Biostrings)           # This function requires the Biostrings package
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:n) {
    chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
    
    # Find out how many gaps there are in chunkseq1aln:
    gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
    
    # Find out how many gaps there are in chunkseq2aln:
    gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
    
    # Calculate how many residues of the first sequence we have printed so far in the alignment:
    seq1alnresidues <- seq1alnresidues + chunksize - gaps1
    
    # Calculate how many residues of the second sequence we have printed so far in the alignment:
    seq2alnresidues <- seq2alnresidues + chunksize - gaps2
    
    if (returnlist == 'FALSE')
    {
      print(paste(chunkseq1aln,seq1alnresidues))
      print(paste(chunkseq2aln,seq2alnresidues))
      print(paste(' '))
    }
  }
  if (returnlist == 'TRUE')
  {
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
  }
}

printPairwiseAlignment(ga_seq, 20)
printPairwiseAlignment(la_seq, 20)

# ====== Compute statistical significance =====

generateSeqsWithMultinomialModel <- function(inputsequence, X) {
  # Change the input sequence into a vector of letters
  require("seqinr") # This function requires the SeqinR package.
  input_seq <- s2c(inputsequence)
  
  # Find the frequencies of the letters in the input sequence "inputsequencevector":
  ln <- length(input_seq)
  tab <- table(input_seq)
  
  # Find the names of the letters in the sequence
  letters <- rownames(tab)
  numletters <- length(letters)
  prob <- numeric() # Make a vector to store the probabilities of letters
  
  for (i in 1:numletters)
  {
    letter <- letters[i]
    count <- tab[[i]]
    prob[i] <- count/ln
  }
  
  # Make X random sequences using the multinomial model with probabilities "probabilities"
  seqs <- numeric(X)
  for (j in 1:X)
  {
    seq <- sample(letters, ln, rep = T, prob = prob) # Sample with replacement
    seq <- c2s(seq)
    seqs[j] <- seq
  }
  # Return the vector of random sequences
  return(seqs)
}
randomseq <- generateSeqsWithMultinomialModel("HEAGAWGHEE" ,1000)
randomseq[1:5]

randomscore <- double(1000)
obj <- c2s(rseq)
for (i in 1:1000) {
  score <- pairwiseAlignment(randomseq[i], obj, substitutionMatrix = "BLOSUM50", 
                  gapOpening = -2, gapExtension = -8, scoreOnly = T)
  randomscore[i] <- score
}
hist(randomscore)
sum(randomscore >= 0)


# ====== Multiple Alignment and Phylogenetic trees =====

# sequence retrieval
retrieveseqs <- function(seqnames,acnucdb)
{
  myseqs <- list()   # Make a list to store the sequences
  require("seqinr")  # This function requires the SeqinR R package
  choosebank(acnucdb)
  for (i in 1:length(seqnames))
  {
    seqname <- seqnames[i]
    print(paste("Retrieving sequence",seqname,"..."))
    queryname <- "query2"
    query <- paste("AC=",seqname,sep="")
    query2 <- query(queryname, query)
    seq <- getSequence(query2$req[[1]]) # Makes a vector "seq" containing the sequence
    myseqs[[i]] <- seq
  }
  closebank()
  return(myseqs)
}

# Call for sequences

seqnames <- c("P06747", "P0C569", "O56773", "Q5VKP1")
seqs <- retrieveseqs(seqnames, "swissprot")

# Construct a multiple sequence alignment "UPGMA Clustering"
msa_seq <- msa::msaMuscle(seqs, 
                          cluster = "upgma")

# Construct a phylogenetic tree analysis



write.fasta(seqs, seqnames, file = "prots.fasta")

aln_clus <- read.alignment(file = "prots1.phy", format = "phylip")

# ====== Filtering the sequences with conserved region ======
filterAlignment <- function(alignment, minpcnongap, minpcid)
{
  # make a copy of the alignment to store the new alignment in:
  newalignment <- alignment
  
  # find the number of sequences in the alignment
  nseqs <- alignment$nb
  
  # empty the alignment in "newalignment")
    for (j in 1:nseqs) { newalignment$seq[[j]] <- "" }
  
  # find the length of the alignment
  alignment_ln <- nchar(alignment$seq[[1]])
  
  # look at each column of the alignment in turn:
  for (i in 1:alignment_ln)
  {
    # see what percent of the letters in this column are non-gaps:
    nongap <- 0
    for (j in 1:nseqs)
    {
      seq_j <- alignment$seq[[j]]
      letter_ij <- substr(seq_j,i,i)
      if (letterij != "-") { nongap <- nongap + 1} # Add a score
    }
    
    pcnongap <- (nongap*100) / nseqs
    
    # Only consider this column if at least minpcnongap % of the letters are not gaps:
    if (pcnongap >= minpcnongap)
    {
      # see what percent of the pairs of letters in this column are identical:
      npairs <- 0; nid <- 0
      
      # find the letters in all of the sequences in this column:
      for (j in 1:(nseqs-1))
      {
        seq_j <- alignment$seq[[j]]
        letter_ij <- substr(seq_j,i,i)
        for (k in (j+1):nseqs)
        {
          seq_k <- alignment$seq[[k]]
          letter_kj <- substr(seq_k,i,i)
          if (letter_ij != "-" && letter_kj != "-")
          {
            npairs <- npairs + 1
            if (letter_ij == letter_kj) { nid <- nid + 1}
          }
        }
      }
      pcid <- (nid*100)/(npairs)
      
      # Only consider this column if at least %minpcid of the pairs of letters are identical:
      if (pcid >= minpcid)
      {
        for (j in 1:nseqs)
        {
          seq_j <- alignment$seq[[j]]
          letter_ij <- substr(seq_j,i,i)
          newalignment_j <- newalignment$seq[[j]]
          newalignment_j <- paste(newalignment_j,letteri_j,sep = "")
          newalignment$seq[[j]] <- newalignment_j
        }
      }
    }
  }
  return(newalignment)
}

# ====== Unrooted phylogenetic tree for protein seq ======
unrootedNJtree <- function(alignment,type)
{
  # this function requires the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    return(mytree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  mytree <- makemytree(mymat)
  # bootstrap the tree
  myboot <- boot.phylo(mytree, mymat, makemytree)
  # plot the tree:
  plot.phylo(mytree,type="u")   # plot the unrooted phylogenetic tree
  nodelabels(myboot,cex=0.7)    # plot the bootstrap values
  mytree$node.label <- myboot   # make the bootstrap values be the node labels
  return(mytree)
}
# ====== Rooted phylogenetic tree for protein seq ======
rootedNJtree <- function(alignment, theoutgroup, type)
{
  # load the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat, outgroup=`theoutgroup`)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    myrootedtree <- root(mytree, outgroup, r=TRUE)
    return(myrootedtree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  myrootedtree <- makemytree(mymat, outgroup=theoutgroup)
  # bootstrap the tree
  myboot <- boot.phylo(myrootedtree, mymat, makemytree)
  # plot the tree:
  plot.phylo(myrootedtree, type="p")  # plot the rooted phylogenetic tree
  nodelabels(myboot,cex=0.7)          # plot the bootstrap values
  mytree$node.label <- myboot   # make the bootstrap values be the node labels
  return(mytree)
}
# ====== Computational Genefinding =====
library(seqinr)
s1 <- "aaaatgcagtaacccatgccc"
matchPattern("atg", s1)

## Obtaining strat and stop codon
# Custom made potential start and stop seqs
findPotentialStartsAndStops <- function(sequence)
{
  # Define a vector with the sequences of potential start and stop codons
  codons <- c("atg", "taa", "tag", "tga") # List of potential start and stop codons
  
  # Find the number of occurrences of each type of potential start or stop codon
  for (i in 1:4)
  {
    codon <- codons[i]
    
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences <- matchPattern(codon, sequence)
    
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions <- attr(occurrences,"start")
   
     # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences <- length(codonpositions)
    
    if (i == 1)
    {
      # Make a copy of vector "codonpositions" called "positions"
      positions <- codonpositions
     
       # Make a vector "types" containing "numoccurrences" copies of "codon"
      types <- rep(codon, numoccurrences)
    }
    else
    {
      # Add the vector "codonpositions" to the end of vector "positions":
      positions <- append(positions, codonpositions, after = length(positions))
      
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types <- append(types, rep(codon, numoccurrences), after = length(types))
    }
  }
  
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  
  # Return a list variable including vectors "positions" and "types":
  mylist <- list(positions,types)
  
  return(mylist)
}
findPotentialStartsAndStops(s1)





# ====== Entrez: database query ======
library(rentrez)

# PMIDs matching citations
cite_obj <- c("proc natl acad sci u s a|1991|88|3248|mann bj|test1|",
              "science|1987|235|182|palmenberg ac|test2|")

entrez_citmatch(cite_obj, db = "pubmed")

# Available databases
entrez_dbs()

# Record from a given database

##Searching for the database and fetching the ids
taxid <- entrez_search(db="taxonomy", term="Osmeriformes")$ids
tax_links <- entrez_db_links("taxonomy")
tax_links

##Searching for the database and fetching the ids and storing in a dataframe
entrez_link(dbfrom="taxonomy", db="pmc", id=taxid)
sra_links <- entrez_db_links("sra")
as.data.frame(sra_links)

# List available search fields for a given database
pmc_fields <- entrez_db_searchable("pmc")
## Check for the field
pmc_fields[["AFFL"]]

## Search for the key in the database
entrez_search(db="pmc", term="Otago[AFFL]", retmax=0)
entrez_search(db="pmc", term="Auckland[AFFL]", retmax=0)

## Store the fields obtained
sra_fields <- entrez_db_searchable("sra")
as.data.frame(sra_fields)

# Retrieve summary information about an NCBI database
entrez_db_summary("swissprot")

# Download data from NCBI databases
katipo <- "Latrodectus katipo[Organism]"
katipo_search <- entrez_search(db="nuccore", term = katipo)

fly_id <- entrez_search(db="taxonomy", term="Drosophila")
## Limititing search
(tax_fields <- entrez_db_searchable("taxonomy"))

## Look for promising field to limit the search, ex - RANK
tax_fields$RANK
entrez_search(db="taxonomy", term="Drosophila & Genus[RANK]")

## Set return type search ENTREZ
kaitpo_seqs <- entrez_fetch(db="nuccore", id = katipo_search$ids, rettype="fasta")
class(kaitpo_seqs)

# Records that match a given term across all NCBI Entrez databases
re_query <- entrez_global_query(term="Heliconius")

# Links to datasets related to records from an NCBI database
pubmed_search <- entrez_search(db = "pubmed", term ="10.1016/j.ympev.2010.07.013[doi]")
linked_dbs <- entrez_db_links("pubmed")
linked_dbs


## Obtain the nucleotide data from pubmed
nucleotide_data <- entrez_link(dbfrom = "pubmed", id = pubmed_search$ids, db ="nuccore") #Sources for the full text of the paper
res <- entrez_link(dbfrom="pubmed", db="", cmd="llinks", id=pubmed_search$ids)
## Retrieve the links associated with the data
linkout_urls(res)

# ====== Computational gene finding =====
tablecode()
genesearch <- function(sequence, codons)
{
  code <- as.list(codons)
  for(i in length(code))
  {
    codon <- code[i]
    occurrence <- matchPattern(codon, sequence)
    codon_pos <- attr(occurrence, "start")
    num_occurrences <- length(codon_pos)
    if(i == 1)
    {
      pos <- codon_pos
      types <- rep(codon, num_occurrences)
    }
    else
    {
      pos <- append(pos, codon_pos, after = length(pos))
      types <- append(types, rep(codon, num_occurrences), after = length(types))
    }
    index <- order(pos)
    pos <- pos[index]
    types <- types[index]
    re_list <- list(pos, types)
    return(re_list)
  }
}