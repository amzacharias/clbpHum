#-------------------------------------------------
# Title: Get a named list of paths given a directory and pattern
# Author: Amanda Zacharias
# Date: 2022-12-31
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Packages -----------------------------------------

# Read dfs function  ---------------------------------
GetNamedPaths <- function(dir, ptrn) {
  #' Get a named list of paths given a directory and pattern
  #'
  #' @param dir path to directory to search (string)
  #' @param ptrn regex pattern for filenames
  #' @return Returns a named list of paths (strings)
  #' @example
  #' pathsList <- GetNamedPaths("./data", ".csv")
  #'
  pathsList <- list.files(dir, pattern = ptrn, full.names = TRUE)
  names(pathsList) <- gsub(ptrn, "", basename(pathsList))
  return(pathsList)
}
