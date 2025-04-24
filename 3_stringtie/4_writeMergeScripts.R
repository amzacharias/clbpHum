#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make individual merge scripts
# Author: Amanda Zacharias
# Date: 2023-07-12
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
baseDir <- file.path(getwd(), "3_stringtie")
# Input ========
refGtfPath <- file.path(
  getwd(), "0_resources", "gencode.v41.annotation.gtf"
)
gtfListsDir <- file.path(baseDir, "gtfLists")
gtfListPath <- list.files(gtfListsDir, pattern = "pass1", full.names = TRUE)

baseScriptPath <- file.path(baseDir, "4_baseMergeScript.sh")

# Output ========
mergeDir <- file.path(baseDir, "4_merge")
system(paste("mkdir", mergeDir))

# Load data -----------------------------------------
# Coldata
coldata <- read.csv("coldata.csv", stringsAsFactors = FALSE, row.names = 1)

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
modifyScript <- function(name, refPath, outPath, inPath, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param refPath path to reference genome gtf file (string)
  #' @param outPath path to output gtf file (string)
  #' @param inPath path in txt file with paths to gtf files (string)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[25] <- paste(script[25], inPath, sep = "")
  script[26] <- paste(script[26], outPath, sep = "")
  script[27] <- paste(script[27], refPath, sep = "")
  # Return
  return(script)
}

processInfo <- function(coldat) {
  #' Prepare information for modifying files,
  #'  so don't have to copy and paste code
  #' @param coldat Dataframe with sample metadata
  #' @return NA, writes bash scripts and messages to console
  #' @example
  #'
  # Modify base script
  sampleLines <- modifyScript(
    name = "all",
    refPath = refGtfPath,
    outPath = file.path(mergeDir, "merged.gtf"),
    inPath = gtfListPath,
    script = baseScript
  )
  # Save
  outFilename <- paste("all", ".merge", ".sh", sep = "")
  fileConn <- file(file.path(mergeDir, outFilename))
  writeLines(sampleLines, fileConn)
  close(fileConn)
}

# Execute ===================================
processInfo(coldataMRNA)
