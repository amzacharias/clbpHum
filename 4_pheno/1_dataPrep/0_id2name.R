#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make id2name file from t2g
# Author: Amanda Zacharias
# Date: 2023-07-06
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# This is basically unnecessary because we have t2g.
# Swapping column names to be more consistent with other code in this project.
#
# Options -----------------------------------------
projDir <- "4_pheno"

# Packages -----------------------------------------
library(dplyr) # 1.1.0


# Pathways -----------------------------------------
# Input ===========
t2gPath <- file.path("3_stringtie", "isoformAnalyzeR", "t2g.csv")

# Output ===========
id2namePath <- file.path(projDir, "1_dataPrep", "id2name.csv")

# Load data -----------------------------------------
t2g <- read.csv(t2gPath, row.names = 1, stringsAsFactors = FALSE)

# Rename columns  -----------------------------------------
id2name <- t2g %>%
  dplyr::rename("isoform_id" = "isoform_id", "gene_id" = "gene_id", "gene_name" = "gene_name")
write.csv(id2name, id2namePath)
