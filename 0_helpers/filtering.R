#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Filtering helper funcions
# Author: Amanda Zacharias
# Date: 2023-02-27
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------



# Pathways -----------------------------------------



# Load data -----------------------------------------
GetMads <- function(df) {
  #' Calculate median absolute deviation (MAD) for each row
  #' @return Returns a dataframe with median absolute deviations
  mads <- data.frame(apply(X = df, MARGIN = 1, FUN = mad))
  colnames(mads) <- "mad"
  return(mads)
}
GetMadCutoff <- function(counts, newPath) {
  #' Plot MADs, so can visualize a good filtering threshold
  #' Adapted from Shreyansh Anand email: `anand@queensu.ca`
  #' @return Writes a pdf file
  # Visualize quantiles to determine the optimal cutoff
  madsAll <- data.frame(id = rownames(counts), mad = GetMads(counts))
  madQuantiles <- quantile(madsAll$mad, probs = seq(0, .9, by = .01))

  madNumGenes <- c()
  quantileLabels <- c()
  for (val in madQuantiles) {
    madNumGenes <- c(madNumGenes, length(subset(madsAll, mad >= val)$mad))
    quantileLabels <- c(quantileLabels, names(madQuantiles)[match(val, madQuantiles)])
  }
  madPlotQuantiles <- data.frame(quantileLabels, madNumGenes)

  # Plot relationship between threshold and number of genes kept
  gplot <- madPlotQuantiles %>%
    ggplot(aes(x = as.numeric(gsub("%", "", quantileLabels)), y = madNumGenes)) +
    geom_point() +
    scale_x_continuous(name = "Percentile", breaks = seq(0, 100, 10), limits = c(0, 100)) +
    ylab("Number of remaning genes") +
    ggtitle("Visualizing MAD cutoffs") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(
    plot = gplot, path = newPath,
    filename = "madQuantiles.pdf",
    width = 170, height = 85, units = "mm"
  )
}
