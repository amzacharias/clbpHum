#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Get the number of transcripts with more non-rhythmic FC
# Author: Amanda Zacharias
# Date: 2024-06-25
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0

# Options -----------------------------------------------
trait <- "sexPheno_FACTOR"
projDir <- "4_pheno"
setwd("absolutePath/clbpHum")

# Packages -----------------------------------------------
library(dplyr) # 1.1.4
# Source -----------------------------------------------
source("0_helpers/GetNamedPaths.R")
source("0_helpers/ReadDfs.R")

# Pathways -----------------------------------------------
# Input ===========
edgeRDir <- file.path(projDir, "2_edgeR")
resDfsPaths <- GetNamedPaths(file.path(edgeRDir, "dataframes", trait), ".sigRes.csv")

# Output ===========

# Load data -----------------------------------------------
resDfs <- ReadDfs(resDfsPaths)

# Subset -----------------------------------------------
posFC.dfs <- list()
for (set in names(resDfs)) {
  cat("\n", set)
  posFC.dfs[[set]] <- resDfs[[set]] %>%
    filter_at(vars(starts_with("logFC")), any_vars( . > 0))
  cat(sprintf(" nrow filt: %s nrow sig: %s", nrow(posFC.dfs[[set]]), nrow(resDfs[[set]])))
}
# Sex Pheno FACTOR:
# RhythmicVsAll nrow filt: 114 nrow sig: 169
#  RhythmicVsCH nrow filt: 9 nrow sig: 39
#  RhythmicVsCL nrow filt: 26 nrow sig: 55
#  RhythmicVsConstant nrow filt: 78 nrow sig: 126
#  RhythmicVsMixed nrow filt: 12 nrow sig: 37
