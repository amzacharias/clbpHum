#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Pathway analysis of DETs
# Author: Amanda Zacharias
# Date: 2023-04-10
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
trait <- "sexPheno_FACTOR" # or sexArythOpioid_FACTOR
projDir <- "4_pheno"
setwd("absolutePath/clbpHum")

# Packages -----------------------------------------
library(dplyr) # 1.1.0
library(gprofiler2) # 0.2.1
library(ggplot2) # 3.4.1
library(UpSetR) # 1.4.0
library(ComplexHeatmap) # 2.2.0

# Source -----------------------------------------------
source("0_helpers/GetNamedPaths.R")
source("0_helpers/ReadDfs.R")
source("0_helpers/enrichPlotHelpers.R")
source("0_helpers/MakeUpset.R")

# Pathways -----------------------------------------
# Input =============
edgeRDir <- file.path(projDir, "2_edgeR")
resDfsPaths <- GetNamedPaths(file.path(edgeRDir, "dataframes", trait), ".sigRes.csv")

# Output =============
gprofilerDir <- file.path(edgeRDir, "gprofiler")
traitGprofDir <- file.path(gprofilerDir, trait)
system(paste("mkdir", gprofilerDir, traitGprofDir))

# Load data -----------------------------------------
resDfs <- ReadDfs(resDfsPaths)

# G:profiler -----------------------------------------
RunGost <- function(queryList, newPath, newUniqPath = "", ToGetLink = FALSE) {
  gostRes <- gost(
    query = queryList,
    organism = "hsapiens", ordered_query = FALSE, # not ordered!
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = TRUE,
    user_threshold = 0.05, correction_method = "g_SCS",
    domain_scope = "annotated", custom_bg = NULL,
    numeric_ns = "", sources = NULL, as_short_link = ToGetLink
  )
  # Save output as file
  if (ToGetLink == FALSE) {
    # Let's see the failed names
    cat("\n\tfailed names:", gostRes$meta$genes_metadata$failed)
    cat("\n\twriting df")
    resDf <- as.data.frame(gostRes$result)
    resDfChr <- apply(resDf, 2, as.character)
    write.csv(resDfChr, newPath)

    # get subset of only interpretable size terms
    resSubDf <- resDf %>% subset(term_size >= 10 & term_size <= 500)
    gostRes$resultSubset <- resSubDf
    uniqResDfChr <- apply(resSubDf, 2, as.character)
    write.csv(uniqResDfChr, newUniqPath)

    return(gostRes) # return to the gostRes object
  } else {
    cat("\n\twriting link")
    write(gostRes, newPath)
  } # done writing/returning results
} # end function

gostList <- list()
for (idx in seq_along(resDfs)) {
  resName <- names(resDfs)[idx]
  query <- resDfs[[resName]] %>%
    pull(gene_name)
  if (length(query) < 6000) { # gprofiler throws error with extremely large queries
    cat("\nQuery length for", resName, "is", length(query))
    gostList[[resName]] <- RunGost(
      query = query,
      ToGetLink = FALSE,
      newPath = file.path(traitGprofDir, paste(resName, "csv", sep = ".")), 
      newUniqPath = file.path(traitGprofDir, paste(resName, "subSize", "csv", sep = "."))
    )
    RunGost(
      query = query,
      ToGetLink = TRUE,
      newPath = file.path(traitGprofDir, paste(resName, "txt", sep = "."))
    )
  } # end if statement
} # end loop through dfs

# RData ------------------------------------------------------------
saveRDS(gostList, file = file.path(traitGprofDir, "gostList.rds"))
gostList <- readRDS(file.path(traitGprofDir, "gostList.rds"))

# Plots ------------------------------------------------------------
# All =========
allDfs <- list()
for (resName in names(gostList)) {
  allDfs[[resName]] <- gostList[[resName]]$result %>%
    arrange(p_value) %>%
    slice_head(n = 10)
}
allLimits <- GetPlotLimits(allDfs)
for (resName in names(allDfs)) {
  MakeBubbleplot(
    df = allDfs[[resName]] %>%
      arrange(p_value) %>%
      slice_head(n = 10),
    newFilename = paste(resName, "all", sep = "."),
    newPath = traitGprofDir,
    title = paste("Top 10 pathways in", resName),
    colorLimits = allLimits$p,
    sizeLimits = allLimits$i,
    wrapNum = 20,
    newHeight = 120,
    newWidth = 100
  )
}
# Subset =======
subDfs <- list()
for (resName in names(gostList)) {
  subDfs[[resName]] <- gostList[[resName]]$resultSubset %>%
    arrange(p_value) %>%
    slice_head(n = 10)
}
subLimits <- GetPlotLimits(subDfs)
for (resName in names(subDfs)) {
  cat("\n", resName, ": ")
  toPlotDf <- subDfs[[resName]]
  if (nrow(toPlotDf) >= 1) {
    cat(nrow(toPlotDf))
    MakeBubbleplot(
      df = toPlotDf,
      newFilename = paste(resName, "subSize", sep = "."),
      newPath = traitGprofDir,
      title = paste("Pathways in", resName),
      colorLimits = subLimits$p,
      sizeLimits = subLimits$i,
      wrapNum = 20,
      newHeight = 120,
      newWidth = 100
    )
  }
}
