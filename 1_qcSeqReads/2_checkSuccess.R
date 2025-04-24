#!/usr/bin/env Rscript
#-------------------------------
# Title: Check success of scripts
# Author: Amanda Zacharias
# Date: 2022-12-17
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# If job not yet completed successfully, write "sbatch <scriptn>" to a text file
# Run these jobs manually, the output file will just
# save you some typing and make it easier to monitor progress !
#
# Options -----------------------------------------

# Packages -----------------------------------------
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
# Input ============
baseDir <- file.path(getwd(), "1_qcSeqReads")
fastqcOutDir <- file.path(baseDir, "fastqcOut")

# Output ============
outFilePath <- file.path(baseDir, "jobsToRun.sh")

# Load data -----------------------------------------
# Coldata
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

# Expected files -----------------------------------------
# Make list of expected out files
MakeExpectedList <- function(directory, coldat) {
  paths <- c()
  for (filename in coldata$filename) { # loop through samples
    for (idx in 1:2) { # loop through forward and rev reads
      readname <- paste(filename, idx, "fastqc", sep = "_")
      path <- file.path(
        directory,
        paste(readname, "html", sep = ".")
      )
      paths <- c(paths, path)
    } # finish looping through 1 and 2
  } # finish looping through samples
  return(paths)
}
expectedFiles <- MakeExpectedList(fastqcOutDir, coldata)

# File check ---------------------------------------------------------
# Check whether file exists and is sufficient in size
# If file isn't good, save to a list of lines to be written
CheckFiles <- function(filesList, fileEnd, tool) {
  badLines <- c("To run:")
  for (path in filesList) { # loop through samples
    if (file.exists(path) == FALSE | file.size(path) < 1000) {
      filename <- gsub(fileEnd, "", basename(path))
      newLine <- paste("sbatch ", tool, filename, ".sh", sep = "")
      badLines <- c(badLines, newLine)
    } # end if statement
  } # end loop through samples
  return(badLines)
}
lines <- CheckFiles(expectedFiles, "_fastqc.html", "fq")

# Write lines ---------------------------------------------------------
fileConn <- file(outFilePath)
writeLines(lines, fileConn)
close(fileConn)
