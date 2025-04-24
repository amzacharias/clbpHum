#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual scripts for StringTie pass 1
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
dataDir <- file.path(getwd(), "0_data")
refGtfPath <- file.path(
  getwd(), "0_resources", "gencode.v41.annotation.gtf"
)
alignedDir <- file.path(getwd(), "2_align", "aligned")
rawFilepaths <- list.files(alignedDir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = ".sort.mrkdup.bam$"
)

baseScriptPath <- file.path(baseDir, "1_baseStringtie.sh")

# Output ========
scriptsDir <- file.path(baseDir, "pass1IndivScripts")
gtfsDir <- file.path(baseDir, "pass1Gtfs")
system(paste("mkdir", scriptsDir, gtfsDir))

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
ModifyScript <- function(name, refPath, outPath, inPath, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param refPath path to reference genome gtf file (string)
  #' @param outPath path to output gtf file (string)
  #' @param inPath path to input .bam file (string)
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
    bamPath <- rawFilepaths[
      match(paste(filename, ".sort.mrkdup.bam", sep = ""), basename(rawFilepaths))
    ]
    # Modify base script
    sampleLines <- ModifyScript(
      name = filename,
      refPath = refGtfPath,
      outPath = file.path(gtfsDir, paste(filename, "gtf", sep = ".")),
      inPath = bamPath,
      script = baseScript
    )
    # Save
    outFilename <- paste("st", filename, ".sh", sep = "")
    fileConn <- file(file.path(scriptsDir, outFilename))
    writeLines(sampleLines, fileConn)
    close(fileConn)
  } # end loop through samples
}

# Execute ===================================
ProcessInfo(coldata)
