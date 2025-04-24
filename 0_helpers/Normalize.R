#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: edgeR normalization helper function
# Author: Amanda Zacharias
# Date: 2023-02-27
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
library(edgeR) # 3.28.1


# Pathways -----------------------------------------



# Load data -----------------------------------------
Normalize <- function(data, toLog = FALSE) {
  #' Perform normalization with edgeR
  #'
  #' @param data Input count matrix
  #' @param toLog Whether to log transform the data
  #' @return Returns a transformed count matrix
  dgeList <- DGEList(data)
  calc <- calcNormFactors(dgeList, method = "TMM")
  if (toLog == FALSE) {
    cpm <- cpm(calc, normalized.lib.sizes = TRUE)
  } else {
    cat("\nNormalizing and log2 transforming counts")
    cpm <- cpm(calc, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  }
  rownames(cpm) <- rownames(data)
  return(cpm)
}
