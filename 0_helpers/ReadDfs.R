#-------------------------------------------------
# Title: Reads in dataframes from a named list of paths
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
ReadDfs <- function(paths) {
  #' Read in dataframes from a named list of paths
  #'
  #' @param paths list of paths to files (strings)
  #' @return Returns a named list of dataframes
  #' @example
  #' dfsList <- ReadDfs(list("foo" = "./foo.csv", "bar" = "./bar.csv"))
  #'
  dfsList <- list()
  for (set in names(paths)) {
    dfsList[[set]] <- read.csv(
      paths[[set]],
      row.names = 1,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  return(dfsList)
}
