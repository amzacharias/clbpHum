#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Prepare lists of candidate genes
# Author: Amanda Zacharias
# Date: 2023-04-07
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
# Options -----------------------------------------
projDir <- "4_pheno"

# Packages -----------------------------------------
library(dplyr) # 1.1.0
library(org.Hs.eg.db) # 3.10.0

# Pathways -----------------------------------------
candidatesDir <- file.path(projDir, "2_edgeRCandidate", "candidates")
goTerms <- data.frame(matrix(ncol = 2,nrow = 0, dimnames = list(NULL, c("goId","goName"))))

## Add neutrophil pro-inflam pathways ===========
neutPaths <- read.table(file.path(candidatesDir, "GONeutInflamTerms.txt"),
  stringsAsFactors = FALSE,
  header = FALSE, sep = "\t", row.names = NULL
)
colnames(neutPaths) <- c("goId", "goName")
goTerms <- rbind(goTerms, neutPaths)

# Get GO candidates using org.Hs.eg.db -----------------------------------------
goGenes <- list()
for (idx in seq_len(nrow(goTerms))) {
  goId <- goTerms[idx, ]$goId
  cat("\n", goId)
  annoRes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keytype = "GOALL", keys = goId,
    columns = c("SYMBOL", "GENENAME", "ENSEMBL")
  )
  if (length(unique(annoRes$SYMBOL)) >= 2) { # some go terms don't have human genes
    goGenes[[goId]] <- annoRes
    fGoId <- gsub(":", "", goId) # colons in filenames can make computers throw errors
    write.csv(annoRes, file.path(candidatesDir, "rawCandidates", paste(fGoId, "csv", sep = ".")))
  } else {
    cat("\n\t", goTerms[idx, ]$goName, "!= >= 2 transcripts")
  }
}
cat(length(goGenes))

keptGoIds <- names(goGenes)
keptGoIds <- goTerms[goTerms$goId %in% keptGoIds, ]
write.csv(keptGoIds, file.path(candidatesDir, "keptGoIds.csv"))

# Pathways -----------------------------------------
candidateFiles <- list.files(file.path(candidatesDir, "rawCandidates"),
  full.names = TRUE, pattern = "\\."
)
candidateNames <- gsub(".txt", "", basename(candidateFiles))
candidateNames <- gsub(".csv", "", candidateNames)

# Load data -----------------------------------------
dfsList <- list()
for (name in candidateNames) {
  if (grepl("GO", name) || grepl("hsa", name)) {
    dfsList[[name]] <- read.csv(
      candidateFiles[grepl(name, candidateFiles)],
      stringsAsFactors = FALSE, row.names = 1, header = TRUE
    )
  } else {
    dfsList[[name]] <- read.table(
      candidateFiles[grepl(name, candidateFiles)],
      stringsAsFactors = FALSE, sep = "\t", row.names = NULL,
      header = TRUE
    )
  }
}

# Unify column name for genes -----------------------------------------
for (dfName in names(dfsList)) {
  cat("\n", dfName)
  if (grepl("GO", dfName) == TRUE) {
    cat("\tGO")
    dfsList[[dfName]] <- dfsList[[dfName]] %>%
      dplyr::rename("gene_name" = "SYMBOL")
  }
}

# Get genes -----------------------------------------
genesList <- list()
for (dfName in names(dfsList)) {
  genes <- dfsList[[dfName]]$gene_name
  # Format and remove duplicates
  genesF <- unique(toupper(genes))
  genesList[[dfName]] <- genesF
  cat("\nNumber genes in", dfName, length(unique(genes)))
}

# Save -----------------------------------------
for (dfName in names(genesList)) {
  write.table(genesList[[dfName]],
    file = file.path(candidatesDir, "cleanCandidates", paste(dfName, "txt", sep = ".")),
    quote = FALSE,
    sep = "\n",
    row.names = FALSE,
    col.names = FALSE
  )
}
