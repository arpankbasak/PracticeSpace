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

# Image Feature detection in tidy format
FeatureImage <- function(mask, ref_image, channel = 2, display = FALSE, path){

    lapply(list("EBImage", "tideyverse", "parallel"), require, character.only = T)
    
    z = channel
    
    # Obain binary labels for mask
    counts.temp <- bwlabel(mask)
    so <- stackObjects(counts.temp, img[,,z], combine = T)
    
    if(display != FALSE) display(colorLabels(counts.temp), method = "raster")

    # Compute the segment spatial features
    a <- as.data.frame(computeFeatures.shape(counts.temp, img[,,z]))
    b <- as.data.frame(computeFeatures.moment(counts.temp, img[,,z]))
    c <- as.data.frame(computeFeatures.haralick(counts.temp, img[,,z]))
    d <- as.data.frame(computeFeatures.basic(counts.temp, img[,,z]))

    dat <- cbind(a, b, c, d)
    dat$total_counts <- max(counts.temp)

    # Annotate SampleID
    sample_id <- str_replace(x, pattern = ".png", replacement = "")
    sample_id <- str_replace(sample_id, pattern = paste0(as.character(path), "/image_stock/"), replacement = "")
    sample_id <- str_replace(sample_id, pattern = "/", replacement = "")

    if(nrow(dat) != 0) {
        dat <- dat %>% add_column(SampleID = sample_id, .before = 1)
        dat <- dat %>% add_column(FeatureID = paste(sample_id, row.names(dat), sep = "_"), .before = 1)
        
        png(filename = paste(path, "/segmented_images/binary_segmentation_", sample_id, ".png", sep = ""), bg = "transparent")
        display(colorLabels(counts.temp), method = "raster")
        dev.off()
        return(dat)
        
    }else {
        message(paste("Sample ", x, " has non features (ERB) detected. Change the blur parameter for stringent masking."))
    }

}

# Setup environment
all_env <- function() {
  ipkg = as.data.frame(installed.packages())
  nix = grep(x = ipkg$NeedsCompilation, pattern = "no", ignore.case = T)
  lapply(unique(ipkg$Package[nix]), library, character.only = T)
}

# Function to extract background from Images normalises the image data for extracting features without considering the background
extractBackground = function(x) {
  require(genefilter)
  x    = log(x)
  loc  = half.range.mode( x )
  left = (x - loc)[ x < loc ]
  wid  = sqrt( mean(left^2) )
  c(loc = loc, 
    wid = wid, 
    thr = loc + 6 * wid)
}
