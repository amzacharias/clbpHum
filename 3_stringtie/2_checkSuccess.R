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
baseDir <- file.path(getwd(), "3_stringtie")
pass1GtfsDir <- file.path(baseDir, "pass1Gtfs")
pass2GtfsDir <- file.path(baseDir, "pass2Gtfs")

# Output ============
outFilePath <- file.path(baseDir, "jobsToRun.sh")

# Load data -----------------------------------------
# Coldata
coldata <- read.csv("coldata.csv",
  row.names = 1, stringsAsFactors = FALSE
)

# Make list of expected out files -----------------------------------------
MakeExpectedList <- function(directory, coldat) {
  paths <- c()
  for (filename in coldat$filename) { # loop through samples
    path <- file.path(
      directory, paste(filename, "gtf", sep = ".")
    )
    paths <- c(paths, path)
  }
  return(paths)
}
pass1ExpectedFiles <- MakeExpectedList(pass1GtfsDir, coldata)
pass2ExpectedFiles <- MakeExpectedList(pass2GtfsDir, coldata)

# Check whether file exists and is sufficient in size -----------------------------------------
# If file isn't good, save to a list of lines to be written
CheckFiles <- function(filesList, passNum, fileEnd, tool) {
  badLines <- c(paste("Pass number is", passNum))
  for (path in filesList) { # loop through samples
    if (file.exists(path) == FALSE | file.size(path) < 1000) {
      sampleId <- gsub(fileEnd, "", basename(path))
      newLine <- paste(
        "sbatch",
        paste(tool,
          paste(sampleId, "sh", sep = "."),
          sep = ""
        )
      )
      badLines <- c(badLines, newLine)
    } # end if statement
  } # end loop through samples
  return(badLines)
}
pass1Lines <- CheckFiles(pass1ExpectedFiles, 1, ".gtf", "st")
pass2Lines <- CheckFiles(pass2ExpectedFiles, 2, ".gtf", "st")

# Write lines -----------------------------------------
fileConn <- file(outFilePath)
writeLines(c(pass1Lines, pass2Lines), fileConn)
close(fileConn)
