#-------------------------------------------------
# Title: Helper: Determine whether 'marked' outliers actually should be removed
# Author: Amanda Zacharias
# Date: 2022-12-31
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Packages -----------------------------------------

# Function -----------------------------------------
ExtractOutliers <- function(combinedOuts) {
  #' Determine if samples are "true" outliers
  #'
  #' @param combinedOuts list of lists
  #' @return Returns a normalized dataframe
  #' @example
  #' outs <- extractOutliers(list('raw' = list, 'norm' = list))
  #' @description Samples are outliers if...
  #   - According to the same outlier metric, they are an
  #       outlier before and after normalization, and/or,
  #   - After normalization, they are marked as an outlier by multiple metrics
  setRawOuts <- combinedOuts[["raw"]]
  setNormOuts <- combinedOuts[["norm"]]
  # 1. Find outliers in raw and normal in same test
  outRawNormNames <- c() # outliers in both raw and normalized counts
  allNormOuts <- c() # collecting for part 2
  for (test in names(setRawOuts)) {
    rawOuts <- names(setRawOuts[[test]]@which) # get outliers in raw
    normOuts <- names(setNormOuts[[test]]@which) # get outliers in norm
    outRawNormNames <- c(outRawNormNames, intersect(rawOuts, normOuts))
    allNormOuts <- c(allNormOuts, normOuts)
  } # loop through tests for part 1
  # 2. Find outliers in 2+ tests of normalized counts
  outNorm2Names <- c()
  for (normOut in allNormOuts) {
    # if outlier occurs 2+ times, add to outliers list
    if (sum(allNormOuts %in% normOut) >= 2) {
      outNorm2Names <- c(outNorm2Names, normOut)
    } # end if statement
  } # end looping through all norm outliers

  # merge outliers vectors
  extractedOutliers <- c(outRawNormNames, outNorm2Names)
  # Return
  return(extractedOutliers)
} # end function
