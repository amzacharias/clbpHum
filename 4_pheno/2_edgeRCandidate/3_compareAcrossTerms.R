#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Compare and contrast all candidate analysis results
# Author: Amanda Zacharias
# Date: 2023-04-29
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# Since this UpSet will be very large, custom making instead of using the helper function
#
#
# Options -----------------------------------------
projDir <- "4_pheno"
analysisName <- "sexPheno_FACTOR"
compareName <- "RhythmicVsAll"

# Packages -----------------------------------------
library(dplyr) # 1.1.0
library(org.Hs.eg.db) # 3.10.0

library(UpSetR) # 1.4.0
library(ComplexHeatmap) # 2.2.0

# Pathways -----------------------------------------
# Input =======
edgeRCandDir <- file.path(projDir, "2_edgeRCandidate")
resDfPaths <- list.files(
  file.path(edgeRCandDir, "dataframes"),
  full.names = TRUE, recursive = TRUE, pattern = paste(compareName, "sigRes", sep = "_")
)
resDfPaths <- resDfPaths[! grepl("archive", resDfPaths)]
resDfPaths <- resDfPaths[
  gsub(".*\\.", "", basename(dirname(resDfPaths))) == analysisName
]
names(resDfPaths) <- gsub("\\..*", "", basename(dirname(resDfPaths)))

id2namePath <- file.path(projDir, "1_dataPrep", "id2name.csv")

keepResDfs <- resDfPaths

# Output =======
compareDir <- file.path(edgeRCandDir, "compare")
analysisCompareDir <- file.path(compareDir, analysisName)
system(paste("mkdir", compareDir, analysisCompareDir))

# Load data -----------------------------------------
keptGoIds <- read.csv(file.path(edgeRCandDir, "candidates", "keptGoIds.csv"),
  row.names = 1, stringsAsFactors = FALSE
)

resDfs <- lapply(
  seq_along(keepResDfs),
  function(n) read.csv(keepResDfs[[n]], row.names = 1, stringsAsFactors = FALSE)
)
names(resDfs) <- names(keepResDfs)

genesList <- lapply(
  seq_along(resDfs),
  function(n) resDfs[[n]]$isoform_id
)
names(genesList) <- names(resDfs)

id2name <- read.csv(id2namePath, row.names = 1, stringsAsFactors = FALSE)

# Format GO ids -----------------------------------------
unformattedGO <- names(genesList)[grepl("GO", names(genesList))]
formattedGO <- gsub("^(.{2})(.*)$", "\\1:\\2", unformattedGO)
names(genesList)[names(genesList) %in% unformattedGO] <- formattedGO

# Get the 11 dataframes with most results -----------------------------------------
# Upset package has a maximum of 15 sets; there are only 11 GO terms, so not a real problem
nrowEach <- vapply(genesList, length, numeric(1))
orderedGenesList <- genesList[rev(order(nrowEach))]
top15genesList <- orderedGenesList[1:11]

# Replace GO ids with GO names  -----------------------------------------
names(top15genesList) <- paste(
  "GO: ", keptGoIds$goName[match(names(top15genesList), keptGoIds$goId)],
  sep = ""
)

# Make upset function  -----------------------------------------
MakeUpset <- function(genesList, nSsLabels = 3, nCsLabels = 3,
                      basename = "upSet", newpath = getwd()) {
  # genesList = A named list of genes ex. list("a" = c(1,2), "b" = c(2,3))
  # nSsLabels = How many labels for the set size axis? Will be nSsLabels + 1
  # nCsLabels = same as nSsLabels, but for combination size
  # basename = base filename for output plot
  # newpath = output directory of file

  # Data preparation
  names(genesList) <- gsub("_", " ", names(genesList)) # format the names of the named list
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
        location = 0,
        just = "left"
      )
    ), # end rowAnnotation()
    right_annotation = NULL,
    show_row_names = FALSE
  ) # end UpSet()
  # Save
  pdf(file.path(newpath, paste(basename, "pdf", sep = ".")),
    width = 9, height = 5
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

MakeUpset(top15genesList,
  basename = paste("upset", compareName, sep = "_"),
  newpath = analysisCompareDir
)

# Get all possible intersection combinations  -----------------------------------------
# Helpful link: https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements
getCombos <- function(inList, nCombo) {
  comboNames <- combn(
    x = names(inList), m = nCombo, FUN = paste0, collapse = "", simplify = FALSE
  )
  combos <- combn(inList, nCombo, simplify = FALSE)
  intersects <- lapply(
    seq_along(combos),
    function(n) Reduce(intersect, combos[[n]]) # use reduce so code is more flexible
  )
  names(intersects) <- comboNames
  return(intersects)
}
## Execute ===========
intersectsList <- lapply(
  seq_along(resDfs), function(n) getCombos(top15genesList, n)
)

for (idx in 2:length(resDfs)) {
  cat("\n", idx)
  for (comboName in names(intersectsList[[idx]])) {
    cat("\n   ", comboName, "\n\t", intersectsList[[idx]][[comboName]])
  }
}
## Save ===========
saveRDS(intersectsList,
  file = file.path(analysisCompareDir, paste(compareName, "intersects.rds", sep = "."))
)
capture.output(intersectsList,
  file = file.path(analysisCompareDir, paste(compareName, "intersects.txt", sep = "."))
)
