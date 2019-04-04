# Loading the list of packages required for bio-computation in R
# Prerequisite a table of bioc packages and base packages
# Follow the example

callForBio <- function(bioc_pkg = NULL, oper_pkgs = NULL, bioC = T)
  {
  if(is.null(bioc_pkg) != T || is.null(oper_pkgs) != T)
  {
    # Sourcing Bioconductor
    source("http://www.bioconductor.org/biocLite.R")
    bioc_pkg <- as.character(bioc_pkg)
    
    # Checking for biocLite call
    if(bioC != F) biocLite(bioc_pkg, suppressUpdates = T)
    
    # Sourcing the operational package
    oper_pkgs <- as.character(oper_pkgs)
    
    # Combining the packages for calling the library function
    pkgs <- c(oper_pkgs, bioc_pkg)
    
    # Installing the operational packages
    install_pkg <- oper_pkgs[!oper_pkgs %in% installed.packages()]
    
    # Calling the library function for each packages
    for(lib in install_pkg) install.packages(lib, dependencies = T, quiet = T)
    lapply(pkgs, library, character.only = T)
    
    print("Successfully loaded the above packages.")
  }
  else
  {print("operational package or bioC package is an empty vector")}
}

# Saturation of intensities especially for plotting heatmaps
saturate <- function(x){
    max <- quantile(x, .99)
    min <- quantile(x, .01)
    
    idx <- (x > max)
    x[idx] <- max
    
    idx <- x < min
    x[idx] <- min
    return(x)
}
