#-------------------------------------------------
# Title: Write Hisat2 Scripts
# Author: Amanda Zacharias
# Date: 2023-01-01
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------


# Pathways -----------------------------------------
dataDir <- file.path(getwd(), "0_data")
refIdxPath <- file.path(getwd(), "0_resources", "index", "idx")

baseDir <- file.path(getwd(), "2_align")
scriptsDir <- file.path(baseDir, "indivScripts")
alignedDir <- file.path(baseDir, "aligned")
summariesDir <- file.path(baseDir, "summaries")

system(paste("mkdir", scriptsDir, alignedDir, summariesDir))

# Load data -----------------------------------------
# Get sample names / filenames
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = F
)
# Raw paths
rawPaths <- list.files(file.path(dataDir, c("BGI", "BGI_2"), "CleanData"),
  full.names = TRUE, recursive = TRUE, pattern = "fq.gz"
)

# Base script that we will modify
baseScript <- read.table(
  file.path(baseDir, "1_baseScript.sh"),
  sep = "\n", blank.lines.skip = FALSE,
  comment.char = "", quote = "\'",
  stringsAsFactors = FALSE
)$V1 %>%
  as.character()

# Modify lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, fwdPath, revPath, idxPath,
                         alignPath, sumPath, batchNum, laneNum, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param fwdPath path to forward input fastq file (string)
  #' @param revPath path to reverse input fastq file (string)
  #' @param idxPath path to reference genome index (string)
  #' @param alignPath path to output aligned files (string)
  #' @param sumPath path to output summary file (string)
  #' @param batchNum sequencing batch (numeric)
  #' @param laneNUm sequencing lane (numeric)
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], name, sep = "")
  script[9] <- paste(script[9], name, ".out", sep = "")
  script[10] <- paste(script[10], name, ".err", sep = "")
  # Content
  script[26] <- paste(script[26], fwdPath, sep = "")
  script[27] <- paste(script[27], revPath, sep = "")
  script[28] <- paste(script[28], idxPath, sep = "")
  script[29] <- paste(script[29], alignPath, sep = "")
  script[30] <- paste(script[30], sumPath, sep = "")
  script[31] <- paste(script[31], batchNum, sep = "")
  script[32] <- paste(script[32], laneNum, sep = "")
  script[33] <- paste(script[33], name, sep = "")
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
    readPaths <- rawPaths[
      match(paste(readNames, "fq.gz", sep = "."), basename(rawPaths))
    ]
    fwdPath <- readPaths[grep("_1", basename(readPaths))]
    revPath <- readPaths[grep("_2", basename(readPaths))]
    # Modify base script
    sampleLines <- ModifyScript(
      name = filename,
      fwdPath = fwdPath,
      revPath = revPath,
      idxPath = refIdxPath,
      alignPath = file.path(alignedDir, filename),
      sumPath = file.path(summariesDir, paste(filename, "txt", sep = ".")),
      batchNum = coldat$batch[idx],
      laneNum = coldat$lane[idx],
      script = baseScript
    )
    # Save
    outFilename <- paste("hi", filename, ".sh", sep = "")
    fileConn <- file(file.path(scriptsDir, outFilename))
    writeLines(sampleLines, fileConn)
    close(fileConn)
  } # end loop through samples
}

# Execute ===================================
ProcessInfo(coldata)
