#-------------------------------------------------
# Title: Write the lists of gtfs
# Author: Amanda Zacharias
# Date: 2022-12-27
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
# Input ============
baseDir <- file.path(getwd(), "3_stringtie")
pass1GtfsDir <- file.path(baseDir, "pass1Gtfs")
pass2GtfsDir <- file.path(baseDir, "pass2Gtfs")

# Output ============
gtfListsDir <- file.path(baseDir, "gtfLists")
system(paste("mkdir", gtfListsDir))

# Load data -----------------------------------------
# Coldata
coldata <- read.csv("coldata.csv",
  row.names = 1, stringsAsFactors = FALSE
)

# Function to write paths -----------------------------------------
WritePaths <- function(coldat, gtfDir, newFilename) {
  lines <- c()
  for (idx in 1:nrow(coldat)) {
    filename <- coldat$filename[idx]
    path <- file.path(gtfDir, paste(filename, "gtf", sep = "."))
    if (grepl("2", newFilename)) {
      lines <- c(lines, paste(filename, path))
    } else {
      lines <- c(lines, path)
    }
  } # end loop through samples
  cat("Writing lines with", length(lines), "lines\n")
  write(lines, file.path(gtfListsDir, newFilename))
}

# Execute -----------------------------------------
WritePaths(coldata, pass1GtfsDir, "pass1List.txt")
WritePaths(coldata, pass2GtfsDir, "pass2List.txt")
