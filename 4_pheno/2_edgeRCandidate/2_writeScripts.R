#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual bash scripts
# Author: Amanda Zacharias
# Date: 2023-04-19
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
projDir <- "4_pheno"
rscriptNames <- list(
  "sexPheno_FACTOR" = "0_sexPheno_FACTOR.R"
)

# Packages -----------------------------------------
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
baseDir <- file.path(getwd(), projDir, "2_edgeRCandidate")
candidatesDir <- file.path(baseDir, "candidates", "cleanCandidates")
candPaths <- list.files(candidatesDir, full.names = TRUE)
names(candPaths) <- gsub(".txt", "", basename(candPaths))

bashDir <- file.path(baseDir, "bash")
system(paste("mkdir", bashDir))

# Load data -----------------------------------------
baseScript <- as.character(
  read.table(file.path(baseDir, "2_baseScript.sh"),
    sep = "\n", blank.lines.skip = FALSE,
    comment.char = "", quote = "\'",
    stringsAsFactors = FALSE
  )$V1
)

# Modify Lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, rscriptPath, candPath, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param rscriptPath Path to Rscript that will be run (string)
  #' @param candPath Path to candidate genes txt file (string)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[24] <- paste(script[24], rscriptPath, sep = "")
  script[25] <- paste(script[25], candPath, sep = "")
  script[26] <- paste(script[26], name, sep = "")
  script[34] <- paste(script[34], getwd(), sep = "")
  # Return
  return(script)
}
ProcessInfo <- function(candPaths, rscriptNames) {
  #' Prepare information for modifying files,
  #'  so don't have to copy and paste code
  #' @param coldat Dataframe with sample metadata
  #' @return NA, writes bash scripts and messages to console
  #' @example
  #'
  # Loop through samples
  for (candName in names(candPaths)) {
    cat("\n", candName)
    for (rName in names(rscriptNames)) {
      cat("\n\t", rName)
      sampleLines <- ModifyScript(
        rscriptPath = rscriptNames[[rName]],
        candPath = candPaths[[candName]],
        name = paste(candName, rName, sep = "."),
        script = baseScript
      )
      # Save
      outFilename <- paste(candName, rName, "sh", sep = ".")
      fileConn <- file(file.path(bashDir, outFilename))
      writeLines(sampleLines, fileConn)
      close(fileConn)
    } # end loop through rNames
  } # end loop through candNames
} # end function

# Execute -----------------------------------------
ProcessInfo(candPaths, rscriptNames)

# To Run
toRun <- tidyr::crossing(names(candPaths), names(rscriptNames))
toRunStr <- paste(toRun$`names(candPaths)`, toRun$`names(rscriptNames)`, "sh", sep = ".")
toRunLines <- paste("sbatch", toRunStr)

fileConn <- file(file.path(baseDir, "jobsToRun.sh"))
writeLines(toRunLines, fileConn)
close(fileConn)
