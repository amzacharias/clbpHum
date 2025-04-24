#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Perform batch effect correction
# Author: Amanda Zacharias
# Date: 2023-02-13
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
# Options -----------------------------------------
setwd("absolutePath/clbpHum")
projDir <- "4_pheno"
source("0_helpers/Normalize.R")
source("0_helpers/vstAndPCA.R")

# Packages -----------------------------------------
library(tidyverse) # 1.3.2
library(data.table) # 1.14.2
library(dplyr) # 1.1.0
library(reshape2) # 1.4.4
library(DESeq2) # 1.26.0
library(ggplot2) # 3.3.6

# system('wget "https://www.bioconductor.org/packages/release/bioc/src/contrib/sva_3.46.0.tar.gz"')
# install.packages("0_resources/pkgs/sva_3.46.0.tar.gz", repos = NULL, type = "source")
library(sva) # 3.36.0

# Pathways -----------------------------------------
baseDir <- file.path(projDir, "1_dataPrep")
cleanDataDir <- file.path(baseDir, "cleanData")
filtDataDir <- file.path(baseDir, "filtData")

combatSeqDir <- file.path(baseDir, "combatSeq")
pcaDir <- file.path(combatSeqDir, "pca")
system(paste("mkdir", combatSeqDir, pcaDir))

# Load data -----------------------------------------
# Read in dfs ======
rawCounts <- fread(
  list.files(filtDataDir, full.names = TRUE, pattern = "raw.counts"),
  stringsAsFactors = FALSE, check.names = FALSE
) %>%
  tibble::column_to_rownames(var = "V1")
coldata <- read.csv(file.path(cleanDataDir, "coldata.csv"), row.names = 1) %>%
  mutate(
    phenotype = as.factor(phenotype),
    sex = as.factor(sex)
  )
# Make sure order of samples is the same ======
rawCounts <- rawCounts %>%
  dplyr::select(all_of(coldata$filename))

# Run combatseq -----------------------------------------
# Including condition ======
correctedCond <- ComBat_seq(
  counts = data.matrix(rawCounts),
  batch = coldata$batch,
  group = coldata$phenotype
)
write.csv(correctedCond, file.path(combatSeqDir, "raw.correctedCond.csv"))
# Found 2 batches
# Using full model in ComBat-seq.
# Adjusting for 3 covariate(s) or covariate level(s)
# Estimating dispersions
# Fitting the GLM model
# Shrinkage off - using GLM estimates for parameters
# Adjusting the data
max(correctedCond) # 8449426

# Read -----------------------------------------
correctedCond <- read.csv(file.path(combatSeqDir, "raw.correctedCond.csv"),
  row.names = 1, stringsAsFactors = FALSE
)

# PCA Plot -----------------------------------------
intGrps <- c(
  "phenotype"
)
vstCorrCond <- GetVst(correctedCond, coldata, "phenotype")
write.csv(assay(vstCorrCond), file.path(combatSeqDir, "vst.corrCond.csv"))
MakePCA(
  vstCorrCond,
  intGrps,
  rep(FALSE, ncol(vstCorrCond)),
  "corrCond.pdf",
  pcaDir
)

# Normalize --------------
write.csv(
  Normalize(correctedCond),
  file.path(combatSeqDir, "norm.correctedCond.csv")
)
