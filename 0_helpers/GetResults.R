#!/usr/bin/env Rscript
#-------------------------------------------------
# Title:
# Author: Amanda Zacharias
# Date: 2023-07-17
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------



# Function-----------------------------------------
GetResults <- function(fitObj, nGenesInput, id2nameObj, idColumn="isoform_id") {
  # Extract results
  topgenes <- topTags(fitObj,
    n = nGenesInput,
    adjust.method = "bonferroni",
    sort.by = "p.value"
  )
  # Extract table of results and format
  # 1. Add molecule's id from rownames
  # 2. Add molecule's name; easier for humans to interpret
  if (idColumn == "isoform_id"){
    topgenesIdName <- topgenes$table %>%
      tibble::rownames_to_column("isoform_id") %>%
      left_join(id2nameObj, by = "isoform_id") # keep all obs in x only
  } else {
    topgenesIdName <- topgenes$table %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(id2nameObj, by = "gene_id") # keep all obs in x only
  }
  # Subset significant results
  sigTable <- topgenesIdName %>%
    subset(FWER < 0.05)
  # Return
  return(list(
    "topgenes" = topgenes,
    "topgenesIdName" = topgenesIdName,
    "sigTable" = sigTable
  ))
}
