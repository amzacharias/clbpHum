#-------------------------------------------------
# Title: Enrichment plot helpers
# Author: Amanda Zacharias
# Date: 2023-01-03
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Packages -----------------------------------------
library(dplyr) # 1.0.9
library(ggplot2) # 3.3.6
library(stringr) # 1.4.0
library(cividis) # 0.2.0

# Get plotting limits ------------------------------------------------------------
GetPlotLimits <- function(dfsList) {
  #' Get limits for color and size scales of plots, so legends are unified
  #'
  #' @param dfsList A named list of g:profiler result dataframes
  #' @return Return a named list of vectors
  #' @example
  #'
  # p_value ======
  pMin <- sapply(dfsList, function(n) min(n$p_value, na.rm = TRUE)) %>% min(na.rm = TRUE)
  pMin <- -log10(pMin)
  pMax <- sapply(dfsList, function(n) max(n$p_value, na.rm = TRUE)) %>% max(na.rm = TRUE)
  pMax <- -log10(pMax)
  # intersection size ========
  iMin <- sapply(dfsList, function(n) min(n$intersection_size, na.rm = TRUE)) %>% min(na.rm = TRUE)
  iMax <- sapply(dfsList, function(n) max(n$intersection_size, na.rm = TRUE)) %>% max(na.rm = TRUE)
  # return
  return(list("p" = c(pMax, pMin), "i" = c(iMin, iMax)))
}

# Bubbleplot -----------------------------------------
MakeBubbleplot <- function(df, newFilename, newPath, title,
                           yName = "Term name", newWidth = 120, newHeight = 185,
                           toAddAsterisk = FALSE, wrapNum = 30,
                           colorLimits = c(0, 5), sizeLimits = c(0, 100)) {
  #' Make a bubbleplot to display enrichment analysis results
  #'
  #' @param df Enrichment analysis results dataframe
  #' @param newFilename Prefix for output file (string)
  #' @param newPath path to output file excluding filename (string)
  #' @param title Title for plot (string)
  #' @param yName Y-axis label for plot (string)
  #' @param newWidth Width of output file (numeric)
  #' @param newHeight Height of output file (numeric)
  #' @param toAddAsterisk Whether to label significant results with an asterisk
  #' @param wrapNum Number of characters to wrap the y-axis names with (numeric)
  #' @param colorLimits Range of values to scale bubbles' colouring by (numerics)
  #' @param sizeLimits Range of values to scale bubbles' sizing by (numerics)
  #' @return Returns a ggplot object
  #' @example
  #' MakeBubbleplot(df, "name", "./plots", "title", "Term name", 5, 6, FALSE, 25)
  if (toAddAsterisk == TRUE) {
    df$sig <- ifelse(df$p_adj < 0.05, "*", "")
  }
  gplot <- df %>%
    # mutate("recall" = intersection_size / term_size) %>%
    ggplot() +
    geom_point(aes(
      x = recall, # recall = intersection / term size
      y = reorder(term_name, -p_value),
      size = intersection_size,
      color = -log10(p_value)
    )) +
    scale_color_cividis(name = "-Log10 adj.\np-value", limits = colorLimits) +
    scale_size(range = c(3, 11), name = "Intersection\nsize", limits = sizeLimits) +
    xlab("Intersection size / term size") +
    scale_y_discrete(
      labels = function(x) str_wrap(x, width = wrapNum),
      name = yName
    ) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    )
  if (toAddAsterisk == TRUE) {
    gplot <- gplot +
      geom_text(aes(
        x = recall, y = reorder(term_name, -p_value),
        label = sig, vjust = -3.2
      ), size = 8)
  }
  ggsave(
    plot = gplot,
    filename = paste("bubble", newFilename, "pdf", sep = "."),
    path = newPath,
    width = newWidth, height = newHeight, units = "mm"
  )
  return(gplot)
}
