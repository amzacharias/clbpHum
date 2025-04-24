#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: UpSet Plot Helper function
# Author: Amanda Zacharias
# Date: 2023-04-14
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
# install.packages("rjson_0.2.20.tar.gz", repos = NULL, type = "source")
# BiocManager::install("ComplexHeatmap")
library(UpSetR) # 1.4.0
library(ComplexHeatmap) # 2.2.0

# Pathways -----------------------------------------


# Function  -----------------------------------------
MakeUpset <- function(genesList, nSsLabels = 3, nCsLabels = 3,
                      basename = "upSet", newpath = getwd()) {
  #' Make an upSet plot
  #'
  #' @param genesList A named list of genes ex. list("a" = c(1,2), "b" = c(2,3))
  #' @param nSsLabels  How many labels for the set size axis? Will be nSsLabels + 1 (numeric)
  #' @param nCsLabels same as nSsLabels, but for combination size (numeric)
  #' @param newFilename Base filename for output plot file (string)
  #' @param newPath Output directory for file (string)
  #' @return Returns a dataframe with enrichment results
  #' @example
  #' inList <- list("a" = c("a", "a", "b", "c"),
  #'                "b" = c("b", "b", "c", "d"),
  #'                "c" = c("c", "c", "d", "e", "f", "g"))
  #' MakeUpset(genesList = inList, nSsLabels = 3, nCsLabels = 3,
  #'           basename = "upSet", newpath = getwd())
  # Data preparation
  m <- make_comb_mat(genesList, mode = "distinct")
  ss <- set_size(m)
  cs <- comb_size(m)
  ssLabels <- pretty(x = c(0, ss), n = nSsLabels)
  csLabels <- pretty(x = cs, n = nCsLabels)

  # Basic plot
  up <- UpSet(m,
    set_order = order(ss),
    comb_order = order(-cs),
    # Right annotation ======
    top_annotation = HeatmapAnnotation(
      "Intersecting transcripts" = anno_barplot(
        cs,
        border = FALSE,
        gp = gpar(fill = "black"),
        height = unit(60, "mm"),
        axis_param = list(
          side = "left",
          at = csLabels,
          labels = csLabels,
          labels_rot = 0
        )
      ), # end anno_plot()
      annotation_name_side = "left",
      annotation_name_offset = unit(15, "mm"),
      annotation_name_rot = 90
    ), # end HeatmapAnnotation()
    # Left annotation ======
    left_annotation = rowAnnotation(
      "# of DETs" = anno_barplot(-ss,
        baseline = 0,
        border = FALSE,
        gp = gpar(fill = "black"),
        width = unit(40, "mm"),
        axis_param = list(
          at = -ssLabels,
          labels = ssLabels,
          labels_rot = 0
        )
      ), # end anno_barplot()
      set_name = anno_text(set_name(m),
        location = 0.5,
        just = "center"
      )
    ), # end rowAnnotation()
    right_annotation = NULL,
    show_row_names = FALSE
  ) # end UpSet()
  # Save
  pdf(file.path(newpath, paste(basename, "pdf", sep = ".")),
    width = 7, height = 5
  )
  up <- draw(up)
  od <- column_order(up)
  decorate_annotation("Intersecting transcripts", {
    grid.text(cs[od],
      x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
      default.units = "native", just = c("left", "bottom"),
      gp = gpar(fontsize = 10, col = "#404040"), rot = 45
    )
  })
  dev.off()
  # Return
  return(up)
} # end function
