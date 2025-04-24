#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual multiqc scripts
# Author: Amanda Zacharias
# Date: 2022-12-17
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
baseDir <- file.path(getwd(), "1_qcSeqReads")
fastqcOutDir <- file.path(baseDir, "fastqcOut")
baseScriptPath <- file.path(baseDir, "baseMultiqcScript.sh")

# Output ============
outScriptPath <- file.path(baseDir, "multiqc.sh")
multiqcOutDir <- file.path(baseDir, "multiqcOut")

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
ModifyScript <- function(name, outDir, inPath, script) {
  #' Modify the base script
  #'
  #' @param name output base filename (string)
  #' @param outDir string; output directory (string)
  #' @param inPath paths to input files (string, sep paths by space)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[28] <- paste(script[28], outDir, sep = "")
  script[29] <- paste(script[29], name, sep = "")
  script[30] <- paste(script[30], inPath, sep = "")
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
  inPath <- fastqcOutDir
  # Modify base script
  sampleLines <- ModifyScript(
    "all", multiqcOutDir, inPath, baseScript
  )
  # Save
  filename <- paste("mq", "all", ".sh", sep = "")
  fileConn <- file(outScriptPath)
  writeLines(sampleLines, fileConn)
  close(fileConn)
}

# Execute ===================================
ProcessInfo(coldata)
