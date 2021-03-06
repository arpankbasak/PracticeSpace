---
title: "Phylogenetic Analysis - Archaebacteria"
author: "Arpan Kumar Basak"
date: '2018-07-18'
---

# Objective:
1. Phylogenetic analysis between the protein samples
2. Phylogenetic analysis across *Archaebacteria* database
3. Phylogenetic statistical inference


```{r Dependencies}
setwd("~/Repository/Artifact/")
```

```{r Data}
query <- readAAStringSet("~/Repository/Artifact/arch_query.fasta", format = "fasta", seek.first.rec = T, use.names = T)
print(query, show = "complete")
seq_name <- data.frame(names(query), row.names = as.character(seq(1:25)))
head(seq_name)
write.table(seq_name, "metadata_2018-07-22.txt", sep = "\t")
query
```
### Computing Multiple sequence alignment
Here, we are aligninging the query sequence sets with `Default` parameters. However, the parameters can be altered for optimum alignment. The multiple sequence alignment is computed by `ClustalOmega` and the output is displayed in `tex` format. The output of the multiple sequence alignment shows the consensus sequence using the logo graph, the shading shows the structural and functional property.

```{r MultipleSequenceAlignment}
align_obj <- msa(inputSeqs = query, method = "ClustalOmega", 
                 type = "protein", order = "input")
# Default printing options
# print(align_obj, show = "complete")

# Print multiple-sequence alignment with properties
msaPrettyPrint(align_obj, output = "tex",
               showNames = "left", shadingMode = "functional", showLogo = "top",
               shadingModeArg = "structure", showConsensus = "bottom", consensusColors = "RedBlue", 
               askForOverwrite = F)

```
```NOTE: The code down below represents the latex source for visulaising the multiple sequence alignment.```
```
$$
\relax
\else
\pdfpagewidth=\paperwidth
\pdfpageheight=\paperheight
\fi
\oddsidemargin=-0.9in
\topmargin=-0.7in
\textwidth=10.8in
\textheight=7.9in

\pagestyle{empty}

\begin{document}
\begin{texshade}{/var/folders/vt/jynnhgxs26n1tjcppvm__zzh0000gn/T//RtmptPauiv/seq1876154b24f.fasta}
\seqtype{P}
\shadingmode[structure]{functional}
\threshold{50}
\showconsensus[RedBlue]{bottom}
\shadingcolors{blues}
\showsequencelogo[chemical]{top}
\hidelogoscale
\shownames{left}
\nameseq{1}{ Archaeoglobus fulgidus ferritin}
\nameseq{2}{ Thermotoga maritima}
\nameseq{3}{ Pyrococcus furiosus}
\nameseq{4}{ tr O26261 O26261 METTH Ferritin like protein (RsgA)}
\nameseq{5}{ E.coli bacterioferritin(BFR ECOLI)}
\nameseq{6}{ Helicobacter pylori ferritin UniProtKB - Q9ZLI1 (FTN HELPJ)}
\nameseq{7}{ Human H-ferritin; UniProtKB - P02794 (FRIH HUMAN)}
\nameseq{8}{ Human L  chain; UniProtKB - P02792 (FRIL HUMAN)}
\nameseq{9}{ Human Mitochondrial ferritin UniProtKB - Q8N4E7 (FTMT HUMAN) }
\nameseq{10}{ Phycomyces blakesleeanus UniProtKB - A0A162U993 (filamentous fungus)}
\nameseq{11}{ sp Q948P5 FRI4 SOYBN Ferritin-4, chloroplastic}
\nameseq{12}{ sp Q39101 FRI1 ARATH Ferritin-1, chloroplastic}
\nameseq{13}{ tr B6DMH6 B6DMH6 9STRA Ferritin (Fragment)(Diatom, Ocean)}
\nameseq{14}{ tr Q52SA8 Q52SA8 TRINI Ferritin (Fragment)(Moth insect)}
\nameseq{15}{ sp P07797 FRI3 LITCT Ferritin, lower subunit}
\nameseq{16}{ sp P22759 BFR AZOVI Bacterioferritin}
\nameseq{17}{ sp Q0P891 DPS CAMJE DNA protection during starvation protein}
\nameseq{18}{ sp P9WNE5 BFRB MYCTU Ferritin BfrB}
\nameseq{19}{ sp P07798 FRI2 LITCT Ferritin, middle subunit}
\nameseq{20}{ tr D0LZ73 D0LZ73 HALO1 Uncharacterized protein}
\nameseq{21}{ tr A0A075ML49 A0A075ML49 CHAVR Ferritin}
\nameseq{22}{ tr Q9HY79 Q9HY79 PSEAE Bacterioferritin}
\nameseq{23}{ tr Q8KBP5 Q8KBP5 CHLTE Ferritin}
\nameseq{24}{ tr H7CGH2 H7CGH2 ULVPE Ferritin}
\nameseq{25}{ tr Q9KVR1 Q9KVR1 VIBCH Ferritin OS Vibrio cholerae serotype O1}
\shownumbering{right}
\showlegend
\end{texshade}
\end{document}
$$
```
Here we have computed multiple sequence alignment by `DECIPHER` algorithm. The algorithm has two levels, one without staggered and staggered alignment respectively. This computation involves most stringent mode of alignment with a threshold of `0.01`. However, these parameters can be altered for optimum alignment structure.

```{r MSA: DECIPHER}
align_obj_dec <- AlignSeqs(query, processors = 3)
align_obj_dec_stg <- StaggerAlignment(align_obj_dec, processors = 3, threshold = 0.01)
# Showing by the consensus sequence
BrowseSeqs(align_obj_dec, highlight = 0)
BrowseSeqs(align_obj_dec_stg, highlight = 0)
```

### Computing the phylogeny
Here, we are constructing the phylogeny from the alignment data. We use the `maximum likelihood method` and created a `neigbour-joining tree`.

```{r Phylogeny}
phy_obj <- msaConvert(align_obj, type = "seqinr::alignment")
phy_obj <- as.phyDat(phy_obj, type = "AA")
dist_obj <- dist.ml(phy_obj, model = "JTT")
phy_tree <- NJ(dist_obj)

# Test for the best possible phylogenetic model
mod <- modelTest(phy_obj, model = "all", multicore = T, tree = phy_tree)
env <- attr(mod, "env")

# Evaluate the best model
fit_start <- eval(get(mod$Model[which.min(mod$BIC)], env), env)
fit_start$tree

# Fit a model suited best for computation (20% invariability)
fit_mod <- pml(phy_tree, phy_obj, model = "Blosum62", k = 4, inv = .2)

# Optimisation of the maximum likelihood (Optimising the distribution and parameters)
fit_opt <- optim.pml(fit_mod, rearrangement  = "stochastic", optInv = T, optGamma = T)
fit_opt$tree$tip.label <- as.character(rownames(seq_name))
fit_opt$tree$tip.label <- names(query)

```

Now that we have an optimised phylogenetic model, we can proceed computing bootstrap values for constructing the phylogenetic tree, `phylogenetic fan`.

```{r Bootstrap Analysis}
# Compute boostrap values
bs_obj = bootstrap.pml(fit_opt, bs = 1000, optNni = T, multicore = T)
bs_obj[[1]]$tip.label

# Plotting the tree
plotBS(midpoint(fit_opt$tree, node.labels = "deleted"), bs_obj, cex = 0.5, p = 70, 
       type = "fan", frame = "circle") # Fan

plotBS(midpoint(fit_opt$tree, node.labels = "deleted"), bs_obj, cex = 0.5, p = 70, 
       type = "phylogram", frame = "circle") # Phylogram

plotBS(midpoint(fit_opt$tree, node.labels = "deleted"), bs_obj, cex = 0.5, p = 70, 
       type = "cladogram", frame = "circle") # Cladogram

plotBS(midpoint(fit_opt$tree, node.labels = "deleted"), bs_obj, cex = 0.5, p = 70, 
       type = "unrooted", frame = "circle") # Unrooted


```

Here, we can plot publication quality phylogenetic plots, by using `ggtree` package.
```{r GGTREE}
# Beatufication of the Phylogenetic tree
ggtree(fit_opt$tree, layout = 'circular') + 
  geom_tiplab(size = 1, angle = 45) 
```

### DECIPHER Pipeline - Phylogenetic analysis
Here, we have computed the phylogeny based on `DECIPHER` pipeline. Considering the stringent alignment, we have computed the distances, considering Terminal Gaps, penalizing `Gap letter match` and `Gap gap match`. However, the parameters can be changed for optimum distance matrix calculation. Following the distance matrix, we have clustered the aligned sequence based on their similarities, using `UPGMA` method.

```{r Phylogeny: Decipher}
# By DECIPHER
con_obj <- ConsensusSequence(align_obj_dec)
names(align_obj_dec_stg) <- as.character(seq(1:25))

# Distance matrix
obj_stg_dist <- DistanceMatrix(myXStringSet = align_obj_dec_stg, type = "dist", 
                               includeTerminalGaps = T, 
                               penalizeGapLetterMatches = T,
                               penalizeGapGapMatches = T,
                               correction = "none")

# Clustering
clust_obj <- IdClusters(obj_stg_dist, method = "UPGMA", cutoff = 0.01, showPlot = T, type = "both", model = "ML")
print(clust_obj)
seq_name$cluster <- clust_obj[[1]]$cluster

```

### Statistical Inference
Here, we are computing the statistical inference of the phylogenetic tree. There are many options to compute at the level of nodes, clades, and tree balance. Also, we can compute the stairs obtained in the phylogenetic tree.

```{r StatisticalInference}

# Computing statistics of the tree

# Average score
avg_mod <- pitchforks(phy_tree, normalise = T)

# Sacking Index
sack_ind <- sackin.phylo(fit_opt$tree, normalise = T)

# Node fraction
node_frac <- nodeImbFrac(fit_opt$tree, threshold = 2)

# Tree Imbalance
tree_imb <- treeImb(fit_opt$tree)

# Staircase fractions
str_obj <- stairs(fit_opt$tree)

slops <- splitTop(fit_opt$tree, dist = 6)

```
---
#####Session Info
```{r SessionInfo}
sessionInfo()
```
