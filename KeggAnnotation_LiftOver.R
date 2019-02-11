# Script for annotation operation
# Arpan Kumar Basak
# KEGG operations especially for Metabolome annotations

# Function calls

# Returns possible lits of kegg compounds on a query in dataframe format with important columns and returns additional dataset including max and min mass of the query search
annoKEGG <- function(query = NULL, type = "compound", output = NULL, bin.width = 1) {
    require("KEGGREST")
    require("tidyverse")
  
  # KEGG query on target    
      comp <- data.frame(compound = keggFind(database = type, query = q))
      
  # Setting parameters
      description <- NULL
      formula <- NULL
      entry <- NULL
      exact.mass <- NULL
      mol.wt <- NULL
      reaction <- NULL
      pathway <- NULL
      enzyme <- NULL
      brite <- NULL
      dblinks <- NULL
      atom <- NULL
      bond <- NULL
      i = 1
      
      # Fetch target metadata for spectral alignment and annotation and returning NA for not annotated
      for(i in 1:nrow(comp)) {
        comp_annot <- keggGet(dbentries = rownames(comp)[i])
        
        formula[i] <- c(ifelse(length(comp_annot[[1]]$FORMULA) != 0, comp_annot[[1]]$FORMULA, "NA"))
        entry[i] <- c(ifelse(length(comp_annot[[1]]$ENTRY) != 0, comp_annot[[1]]$ENTRY, "NA"))
        exact.mass[i] <- c(ifelse(length(comp_annot[[1]]$EXACT_MASS) != 0, comp_annot[[1]]$EXACT_MASS, "NA"))
        mol.wt[i] <- c(ifelse(length(comp_annot[[1]]$MOL_WEIGHT) != 0, comp_annot[[1]]$MOL_WEIGHT, "NA"))
        reaction[i] <- c(ifelse(length(comp_annot[[1]]$REACTION) != 0, comp_annot[[1]]$REACTION, "NA"))
        pathway[i] <- c(ifelse(length(comp_annot[[1]]$PATHWAY) != 0, comp_annot[[1]]$PATHWAY, "NA"))
        enzyme[i] <- c(ifelse(length(comp_annot[[1]]$ENZYME) != 0, comp_annot[[1]]$ENZYME, "NA"))
        brite[i] <- c(ifelse(length(comp_annot[[1]]$BRITE) != 0, comp_annot[[1]]$BRITE, "NA"))
        dblinks[i] <- c(ifelse(length(comp_annot[[1]]$DBLINKS) != 0, comp_annot[[1]]$DBLINKS, "NA"))
        atom[i] <- c(ifelse(length(comp_annot[[1]]$ATOM) != 0, comp_annot[[1]]$ATOM, "NA"))
        bond[i] <- c(ifelse(length(comp_annot[[1]]$BOND) != 0, comp_annot[[1]]$BOND, "NA"))
        comp_annot <- NULL
      }
      
      # Table for the query target
      temp <- cbind.data.frame(comp, formula, entry, exact.mass = as.numeric(exact.mass), 
                               mol.wt = as.numeric(mol.wt), reaction, pathway, enzyme, brite, dblinks, atom, bond)
      
      # For long range of Exact mass within the hits
      min_r <- min(as.numeric(temp$exact.mass), na.rm = T)
      max_r <- max(as.numeric(temp$exact.mass), na.rm = T)
      
      # set bin width for scanning exact mass within the range
      # parameters bin.width
      bin_space <- as.integer(seq(from = min_r[1], to = max_r[1], by = bin.width))
      
      # Obtain a long range mass table within the minimum and maximum range of targets
      compound_exact_mass = keggFind("compound", bin_space, option = "exact_mass")
      long_r <- as.data.frame(compound_exact_mass)
      row.names(long_r) <- str_replace(row.names(long_r), pattern = "cpd:", replacement = "")
      long_r <- cbind(compound_id = row.names(long_r), long_r)
      
      # Returns an object with target ersults, and range results
      res <- list(target = temp, 
                     range_min = min_r, 
                     range_max = max_r, 
                      range_data = long_r,
                     range_hits = table(long_r))
      
      # Documentation for the downloaded informations
      if(!is.null(output)){
        
        write_delim(temp, delim = "\t", path = paste0(output,"/",query,".txt"))
        write_delim(long_r, delim = "\t", path = paste0(output,"/",query,"_range.txt")) 
      }
      return(res)
}

# Test run
# f <- annoKEGG(query = q, type = "compound", output = "./GSL_neg_ion/")

# Fetch kegg compounds within the range of mass
fetchKEGG <- function(msdata, feature_data, mass_col){
    require("KEGGREST")
    require("tidyverse")
  obj <- NULL
  form_obj <- NULL
  info_obj <- NULL
  loc <- NULL
  i = 1
  j = 1
  for (i in 1:nrow(msdata)) {
    for (j in 1:nrow(feature_data)) {
      if (round(msdata[i,mass_col], digits = 0) == round(as.numeric(feature_data$exact.mass[j]), digits = 0))
      {
        obj[i] <- c(feature_data$entry[j])
        form_obj[i] <- c(feature_data$formula[j])
        info_obj[i] <- c(feature_data$compound[j])
        break();
      }
    }
  }
  hits_obj <- data.frame(kegg_id = obj, kegg_compounds = info_obj, row.names = rownames(msdata))
  return(hits_obj)
}

# keg_hits_tar <- fetchKEGG(msdata = mdat, feature_data = fdat, mass_col = "Detected.Mass.Ave.")
