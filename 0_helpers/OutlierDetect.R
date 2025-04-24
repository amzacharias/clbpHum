#-------------------------------------------------
# Title: Helper: Detect outliers with arrayQualityMetrics
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
OutlierDetect <- function(counts, coldat, outDirname, intGroup) {
  #' Perform outlier detection using arrayQualityMetrics
  #'
  #' @param counts count matrix
  #' @param coldat metadata for samples
  #' @param outDirname folder for output files
  #' @param intGroup the variable of interest (string);
  #'  matches associated coldat column name
  #' @return Returns information about labelled "outliers"
  #' @example
  #' outliers <- outlierDetect(counts, coldata, "./qcOut")
  # Do outlier detection
  phenoData <- AnnotatedDataFrame(data = coldat)
  eset <- ExpressionSet(assayData = counts, phenoData = phenoData)
  arrayQualityMetrics(
    expressionset = eset,
    intgroup = intGroup,
    outdir = outDirname,
    do.logtransform = FALSE,
    force = TRUE
  )
  # Advanced outlier detection
  preparedData <- prepdata(
    expressionset = eset,
    intgroup = intGroup,
    do.logtransform = FALSE
  )
  box <- aqm.boxplot(preparedData)
  heat <- aqm.heatmap(preparedData)
  ma <- aqm.maplot(preparedData)
  qmetrics <- list("Boxplot" = box, "Heatmap" = heat, "MA" = ma)
  outliers <- list(
    "Boxplot" = box@outliers,
    "Heatmap" = heat@outliers,
    "MA" = ma@outliers
  )
  capture.output(outliers, file = file.path(outDirname, "advancedOutlierData.txt"))
  return(outliers)
}
