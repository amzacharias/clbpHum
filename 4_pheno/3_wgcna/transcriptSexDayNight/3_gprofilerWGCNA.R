#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: G: profiler analysis of modules
# Author: Amanda Zacharias
# Date: 2023-05-29
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
prefix <- "transcriptSexDayNight"
projDir <- "4_pheno"
setwd("absolutePath/clbpHum")
source("0_helpers/enrichPlotHelpers.R")
source("0_helpers/GetNamedPaths.R")
source("0_helpers/ReadDfs.R")

# Packages -----------------------------------------
library(gprofiler2) # 0.2.1
library(ggplot2) # 3.4.1

# Pathways -----------------------------------------
wgcnaDir <- file.path(projDir, "3_wgcna", prefix)
geneInfoPaths <- GetNamedPaths(file.path(wgcnaDir, "geneInfo"), ".csv")
modTraitPaths <- GetNamedPaths(file.path(wgcnaDir, "modTrait"), ".isRhythmic.csv")

# Output paths =====
plotsDir <- file.path(wgcnaDir, "plots")
gprofDir <- file.path(wgcnaDir, "gprofiler")
gprofPlotsDir <- file.path(plotsDir, "gprofiler")

# Make output directories
system(paste("mkdir", plotsDir, gprofDir, gprofPlotsDir))

# Load data -----------------------------------------
geneInfoList <- ReadDfs(geneInfoPaths)
modTraitRes <- ReadDfs(modTraitPaths)

sigModTrait <- lapply(
  1:length(modTraitRes),
  function(n) modTraitRes[[n]] %>% subset(pval < 0.05)
)
names(sigModTrait) <- names(modTraitRes)

# G: profiler -----------------------------------------
RunGost <- function(query_list, to_get_link, out_path) {
  gost_res <- gost(
    query = query_list,
    organism = "hsapiens", ordered_query = FALSE, # not ordered!
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = TRUE,
    user_threshold = 0.05, correction_method = "g_SCS",
    domain_scope = "annotated", custom_bg = NULL,
    numeric_ns = "", sources = NULL, as_short_link = to_get_link
  )
  # Save output as file
  if (to_get_link == FALSE) {
    # Let's see the failed names
    cat("\n\tfailed names:", gost_res$meta$genes_metadata$failed)
    cat("\n\twriting df")
    res_df <- as.data.frame(gost_res$result)
    res_df_chr <- apply(res_df, 2, as.character)
    write.csv(res_df_chr, out_path)

    # get subset of only interpretable size terms
    resSub_df <- res_df %>% subset(term_size >= 10 & term_size <= 500)
    gost_res$resultSubset <- resSub_df

    return(gost_res) # return to the gost_res object
  } else {
    cat("\n\twriting link")
    write(gost_res, out_path)
  } # done writing/returning results
} # end function

gostList <- list()
for (set in names(sigModTrait)) {
  setGostList <- list()
  for (idx in 1:length(unique(sigModTrait[[set]]$modName))) {
    modName <- gsub("ME", "", sigModTrait[[set]]$modName[idx])
    query <- geneInfoList[[set]] %>%
      subset(module_color == modName) %>%
      pull(gene_name)
    cat("\nQuery length for", set, modName, "is", length(query))
    gostRes <- RunGost(
      query = query,
      to_get_link = FALSE,
      out_path = file.path(gprofDir, paste(set, modName, "csv", sep = "."))
    )
    RunGost(
      query = query,
      to_get_link = TRUE,
      out_path = file.path(gprofDir, paste(set, modName, "txt", sep = "."))
    )
    # Output result, to be added to setGostList
    setGostList[[modName]] <- gostRes
  } # end loop through modules
  gostList[[set]] <- setGostList
} # end loop through sets

# RData -----------------------------------------
saveRDS(gostList, file = file.path(gprofDir, "gostList.rds"))
gostList <- readRDS(file.path(gprofDir, "gostList.rds"))

# Make plots -----------------------------------------
# Function ==========
MakePlots <- function(gostList, dfType, prefix) {
  #' Description
  #'
  #' @param gostList Object with results from G: Profiler
  #' @param dfType Whether to extract "results" or "resultsSubset" dataframes
  #' @param prefix Prefix for output filenames\
  #' @return Writes pdf files with bubble plots
  #' @example
  #' MakePlots(gostList, "results", "all")
  # Make plot dataframes ======
  dfsList <- list()
  for (set in names(gostList)) {
    setDfs <- list()
    for (modName in names(gostList[[set]])) {
      setDfs[[modName]] <- gostList[[set]][[modName]][[dfType]] %>%
        arrange(p_value) %>%
        slice_head(n = 10)
    }
    dfsList[[set]] <- setDfs
  }
  limitsList <- lapply(
    1:length(dfsList),
    function(n) GetPlotLimits(dfsList[[n]])
  )
  names(limitsList) <- names(dfsList)
  # Make plots =====
  system(paste("mkdir", file.path(gprofPlotsDir, prefix)))
  for (set in names(dfsList)) {
    for (modName in names(dfsList[[set]])) {
      MakeBubbleplot(
        df = dfsList[[set]][[modName]],
        newFilename = paste(set, modName, prefix, sep = "."),
        newPath = file.path(gprofPlotsDir, prefix),
        title = paste("Pathways in ", gsub("_", " ", set), ", ", modName, sep = ""),
        colorLimits = limitsList[[set]]$p,
        sizeLimits = limitsList[[set]]$i,
        wrapNum = 25,
        newHeight = 120,
        newWidth = 100
      )
    }
  }
}

# Execute ============
MakePlots(gostList, "result", "all")
MakePlots(gostList, "resultSubset", "subSize")
