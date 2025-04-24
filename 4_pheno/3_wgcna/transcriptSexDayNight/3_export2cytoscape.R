#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Export modules to cytoscape
# Author: Amanda Zacharias
# Date: 2023-05-30
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
setwd("absolutePath/clbpHum") # use if running remotely
prefix <- "transcriptSexDayNight"
projDir <- "4_pheno"

# Packages -----------------------------------------
library(WGCNA) # 1.70-3
library(dplyr) # 1.1.0

# Pathways -----------------------------------------
# Inputs =====
vstPath <- file.path(projDir, "1_dataPrep", "combatSeq", "vst.corrCond.csv")

wgcnaDir <- file.path(projDir, "3_wgcna", prefix)
rdataDir <- file.path(wgcnaDir, "rdata")
dissTOMpath <- file.path(rdataDir, paste(prefix, "dissTOMList.rds", sep = "."))
mergedListPath <- file.path(rdataDir, paste(prefix, "mergedList.rds", sep = "."))
geneInfoDir <- file.path(wgcnaDir, "geneInfo")
geneInfoPaths <- list.files(geneInfoDir, full.names = TRUE)
names(geneInfoPaths) <- gsub(".csv", "", basename(geneInfoPaths))
modTraitResPaths <- list.files(file.path(wgcnaDir, "modTrait"),
  pattern = "isRhythmic", full.names = TRUE
)
names(modTraitResPaths) <- gsub(".isRhythmic.csv", "", basename(modTraitResPaths))

# Outputs =====
cytoscapeDir <- file.path(wgcnaDir, "cytoscape")
system(paste("mkdir", cytoscapeDir))

# Load data -----------------------------------------
vstCounts <- read.csv(vstPath, row.names = 1, stringsAsFactors = FALSE) %>% t()
dissTOMList <- readRDS(dissTOMpath) # Can't be loaded into RStudio because too large !!!!!
mergedList <- readRDS(mergedListPath)
geneInfoDfs <- lapply(
  1:length(geneInfoPaths),
  function(n) read.csv(geneInfoPaths[[n]], stringsAsFactors = FALSE, row.names = 1)
)
names(geneInfoDfs) <- names(geneInfoPaths)

modTraitDfs <- lapply(
  1:length(modTraitResPaths),
  function(n) read.csv(modTraitResPaths[[n]], stringsAsFactors = FALSE, row.names = 1)
)
names(modTraitDfs) <- names(modTraitResPaths)

# Export -----------------------------------------
ExportCyto <- function(dissTOM, modColors, mod2plot, counts, geneInfo, set, plotThresh = 0.005) {
  # TOM Matrix
  TOM <- 1 - dissTOM
  # Select_modules
  targModLogic <- modColors %in% mod2plot
  cat("\n\tNumber of genes in modules", nrow(geneInfo %>% subset(module_color == mod2plot)))
  cat("\n\tNumber of genes in matched", sum(targModLogic))
  # Select corresponding features
  modTranscripts <- colnames(counts)[targModLogic]
  matchId2Name <- match(modTranscripts, geneInfo$isoform_id) # where is x in y
  modGeneNames <- geneInfo$gene_name[matchId2Name]
  isDet <- geneInfo$is_det[matchId2Name]
  isHub <- geneInfo$is_hub[matchId2Name]
  # Subset TOM
  modTOM <- TOM[targModLogic, targModLogic]
  dimnames(modTOM) <- list(modTranscripts, modTranscripts)
  # Export
  cyt <- exportNetworkToCytoscape(
    modTOM,
    edgeFile = file.path(cytoscapeDir, paste(set, mod2plot, "CytoscapeInput-edges.txt", sep = ".")),
    nodeFile = file.path(cytoscapeDir, paste(set, mod2plot, "CytoscapeInput-nodes.txt", sep = ".")),
    weighted = TRUE, threshold = plotThresh,
    nodeNames = modTranscripts, altNodeNames = modGeneNames,
    nodeAttr = data.frame(
      "modColor" = rep(mod2plot, nrow(modTOM)),
      "isDet" = isDet,
      "isHub" = isHub
    )
  )
  cat("\nDone!")
}

# Individual modules
for (set in names(geneInfoDfs)) {
  for (mod in modTraitDfs[[set]] %>%
    subset(pval < 0.05) %>%
    pull(modName)) {
    pThresh <- ifelse((set == "D" & mod == "MEbrown") |
      (set == "N" & mod == "MElightcyan1"),
    0.05, 0.005
    )
    cat("\n", set, mod, pThresh)
    ExportCyto(
      dissTOM = dissTOMList[[set]],
      modColors = mergedList[[set]]$mergedColors,
      mod2plot = sub("ME", "", mod),
      counts = vstCounts,
      geneInfo = geneInfoDfs[[set]],
      set = set,
      plotThresh = pThresh
    )
  } # end loop through individual modules
} # end loop through networks


cat("\nDone!")
