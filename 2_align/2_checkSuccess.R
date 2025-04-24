#-------------------------------------------------
# Title: Check success of hisat scripts
# Author: Amanda Zacharias
# Date: 2023-01-01
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
library(dplyr) # 1.0.7

# Pathways -----------------------------------------
baseDir <- file.path(getwd(), "2_align")
alignedDir <- file.path(baseDir, "aligned")
outFilePath <- file.path(baseDir, "jobsToRun.sh")

# Load data -----------------------------------------
# Get sample names / filenames
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

# Make list of expected out files -----------------------------------------
makeExpectedList <- function(directory, coldata) {
  paths <- c()
  for (sample in coldata$filename) { # loop through samples
    path <- file.path(
      directory,
      paste(sample, "sort.mrkdup.bam.bai", sep = ".")
    )
    paths <- c(paths, path)
  } # finish looping through samples
  return(paths)
}
expectedFiles <- makeExpectedList(alignedDir, coldata)

# Check files -----------------------------------------
# Check whether file exists and is sufficient in size
# If file isn't good, save to a list of lines to be written
checkFiles <- function(paths, fileEnd, tool) {
  badLines <- c("Jobs to run")
  for (path in paths) { # loop through samples
    if (file.exists(path) == FALSE | file.size(path) < 100000) {
      filename <- gsub(fileEnd, "", basename(path))
      new_line <- paste(
        "sbatch",
        paste(tool,
          paste(filename, "sh", sep = "."),
          sep = ""
        )
      )
      badLines <- c(badLines, new_line)
    } # end if statement
  } # end loop through samples
  return(badLines)
}
runLines <- checkFiles(expectedFiles, ".sort.mrkdup.bam.bai", "hi")

# Write lines -----------------------------------------
fileConn <- file(outFilePath)
writeLines(runLines, fileConn)
close(fileConn)
