#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual fastqc scripts
# Author: Amanda Zacharias
# Date: 2022-12-17
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# Running fastqc can be finicky, refer to the following link when troubleshooting:
# https://cac.queensu.ca/wiki/index.php/Software:Frontenac
#
#
# Options -----------------------------------------

# Packages -----------------------------------------
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
# Input ============
dataDir <- file.path(getwd(), "0_data")
rawFilepaths <- list.files(file.path(dataDir, c("BGI", "BGI_2"), "CleanData"),
  full.names = TRUE, recursive = TRUE, pattern = "fq.gz"
)

baseDir <- file.path(getwd(), "1_qcSeqReads")
baseScriptPath <- file.path(baseDir, "baseFastqcScript.sh")

# Output ============
scriptDir <- file.path(baseDir, "fastqcScripts")
fastqcOutDir <- file.path(baseDir, "fastqcOut")

system(paste("mkdir", scriptDir, fastqcOutDir))

# Load data -----------------------------------------
# Coldata
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)

# Base script
baseScript <- read.table(
  baseScriptPath,
  sep = "\n", blank.lines.skip = FALSE,
  comment.char = "", quote = "\'",
  stringsAsFactors = FALSE
)$V1 %>%
  as.character()

# Modify lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, inPath, outPath, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param inPath path to input fastq file (string)
  #' @param outPath path to output fastqc folder (string)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[27] <- paste(script[27], inPath, sep = "")
  script[28] <- paste(script[28], outPath, sep = "")
  # Return
  return(script)
}
ProcessInfo <- function(coldat) {
  #' Prepare information for modifying files,
  #'  so don't have to copy and paste code
  #' @param coldat Dataframe with sample metadata
  #' @return NA, writes bash scripts and messages to console
  #' @example
  #'
  # Loop through samples
  for (idx in 1:nrow(coldat)) {
    # Path info
    filename <- coldat[idx, ]$filename
    readNames <- paste(filename, paste0(c("1", "2")), sep = "_")
    readPaths <- rawFilepaths[
      match(paste(readNames, "fq.gz", sep = "."), basename(rawFilepaths))
    ]
    for (jdx in 1:length(readPaths)) {
      readName <- gsub(".fq.gz", "", basename(readPaths[jdx]))
      cat("\n\t", readName)
      # Modify base script
      sampleLines <- ModifyScript(
        readName, readPaths[jdx], fastqcOutDir, baseScript
      )
      # Save
      outFilename <- paste("fq", readName, ".sh", sep = "")
      fileConn <- file(file.path(scriptDir, outFilename))
      writeLines(sampleLines, fileConn)
      close(fileConn)
    } # end loop through raw files
  } # end loop through samples in a tissue
}

# Execute ===================================
ProcessInfo(coldata)
