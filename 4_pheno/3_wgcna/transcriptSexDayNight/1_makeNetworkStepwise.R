#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make a signed network with WGCNA
# Author: Amanda Zacharias
# Date: 2023-07-18
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
setwd("absolutePath/clbpHum")
prefix <- "transcriptSexDayNight"
projDir <- "4_pheno"

# Packages -----------------------------------------
library(dplyr) # 1.1.10
library(ggplot2) # 3.3.6
library(cividis) # 0.2.0
library(cowplot) # 1.1.1
library(WGCNA) # 1.70-3

# Pathways -----------------------------------------
# Input ===========
dataPrepDir <- file.path(projDir, "1_dataPrep")
coldataPath <- file.path(dataPrepDir, "cleanData", "coldata.csv")
vstPath <- file.path(dataPrepDir, "combatSeq", "vst.corrCond.csv")
id2namePath <- file.path(dataPrepDir, "id2name.csv")
deaResPath <- file.path(
  projDir, "2_edgeR", "dataframes",
  "sexPheno", "RhythmicVsAll_sigRes.csv"
)

# Output ===========
wgcnaDir <- file.path(projDir, "3_wgcna", prefix)
rdataDir <- file.path(wgcnaDir, "rdata")
geneInfoDir <- file.path(wgcnaDir, "geneInfo")
plotsDir <- file.path(wgcnaDir, "plots")
sftPlotsDir <- file.path(plotsDir, "sft")
netPlotsDir <- file.path(plotsDir, "net")
system(paste(
  "mkdir", rdataDir, plotsDir,
  sftPlotsDir, netPlotsDir, geneInfoDir
))

# Load data -----------------------------------------
# Coldata ======
coldata <- read.csv(coldataPath, row.names = 1, stringsAsFactors = FALSE)

# Counts results ========
counts <- read.csv(vstPath, row.names = 1, check.names = FALSE) %>%
  dplyr::select(all_of(coldata$filename))

# Subset coldata and counts into batch groups -----------------------------------------
coldataList <- list()
for (set in unique(coldata$timepoint)){
  coldataList[[set]] <- coldata %>% subset(timepoint == set)
}

# Screen genes after stratification -----------------------------------------
badGenes <- c()
for (set in names(coldataList)) {
  df <- counts %>%
    dplyr::select(all_of(coldataList[[set]]$filename)) %>%
    t() %>%
    as.data.frame()
  setBadGenes <- colnames(df)[!goodSamplesGenes(df)$goodGenes]
  badGenes <- c(badGenes, setBadGenes)
  cat("\n", set, setBadGenes, "\n")
}
# N
# D

countsList <- list()
for (set in names(coldataList)) {
  df <- counts %>%
    dplyr::select(all_of(coldataList[[set]]$filename)) %>%
    t() %>%
    as.data.frame()
  countsList[[set]] <- df %>% dplyr::select(-all_of(badGenes))
}

# This count format is more compatible with WGCNA functions
multiExpr <- list(
  "N" = list("data" = countsList[["N"]]),
  "D" = list("data" = countsList[["D"]])
)

# saveRDS -----------------------------------------
saveRDS(countsList, file = file.path(rdataDir, paste(prefix, "countsList.rds", sep = ".")))
saveRDS(multiExpr, file = file.path(rdataDir, paste(prefix, "multiExpr.rds", sep = ".")))
rm(coldata, counts, coldataList, multiExpr)
countsList <- readRDS(file.path(rdataDir, paste(prefix, "countsList.rds", sep = ".")))

# Soft-threshold power -----------------------------------------
# Functions =========
EvalThresholds <- function(countDf) {
  #' Evalulate soft-thresholding powers for a signed network
  #'
  #' @param countDf A gene count matrix (rows = samples)
  #' @return Returns an object with results
  #' @example
  #' sftRes <- EvalThresholds(countMat)
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  # Analyze impact of soft thresholds
  sft <- pickSoftThreshold(countDf,
    powerVector = powers, RsquaredCut = 0.85,
    verbose = 5, networkType = "signed",
    corFnc = "cor",
    corOptions = list(use = "p")
  )
  cat("\nPower estimate is", sft$powerEstimate, "\n")
  return(sft)
}
PlotThresholds <- function(sft, prefix) {
  #' Evalulate soft-thresholding powers for a signed network
  #'
  #' @param sft Soft-thresholding object
  #' @param prefix Prefix for output file (string)
  #' @return Writes a pdf file
  #' @example
  #' PlotThresholds(sftObj, "dataset")
  # 1. Scale independence
  siPlot <- sft$fitIndices %>%
    ggplot(aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
    geom_hline(yintercept = 0.85, color = cividis(20)[10], linetype = "dashed") +
    geom_point() +
    geom_text(label = sft$fitIndices$Power, hjust = 0.5, vjust = -0.7) +
    scale_y_continuous(breaks = seq(-1, 1, by = 0.1)) +
    xlab("Soft threshold (power)") +
    ylab("Scale free topology model fit, signed R^2") +
    ggtitle("Scale independence") +
    theme_bw()
  # 2. Mean connectivity
  mcPlot <- sft$fitIndices %>%
    ggplot(aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_text(label = sft$fitIndices$Power, hjust = 0.5, vjust = -0.7) +
    xlab("Soft threshold (power)") +
    ylab("Mean connectivity") +
    ggtitle("Mean connectivity") +
    theme_bw()
  # 3. Combine into a grid
  comboPlot <- plot_grid(siPlot, mcPlot, labels = c("i", "ii"), label_size = 12)
  # 4. Save the combined plot
  filename <- paste(prefix, "pdf", sep = ".")
  ggsave(
    plot = comboPlot, filename = filename, path = sftPlotsDir,
    width = 185, height = 138.75, units = "mm"
  )
}

sftResList <- list()
for (set in names(countsList)){
  cat("\n", rep("-", 10), "Working on", set, rep("-", 10), "\n")
  sftResList[[set]] <- EvalThresholds(countsList[[set]])
  PlotThresholds(sftResList[[set]], set)
}

# After inspection, define the thresholds for each tissue ===========
# Manual says, "if fit index fails to reach values above 0.8 for reasonable powers
# (< 15 for unsigned or signed-hybrid networks, & < 30 for signed networks) ....
powerList <- list(
  "N" = 12,
  "D" = 12
)

# Save data ============
saveRDS(sftResList, file = file.path(rdataDir, paste(prefix, "sftResList.rds", sep = ".")))
rm(sftResList)
saveRDS(powerList, file = file.path(rdataDir, paste(prefix, "powerList.rds", sep = ".")))
powerList <- readRDS(file.path(rdataDir, paste(prefix, "powerList.rds", sep = ".")))

# Make network ----------------------------------------------------------
# Similarity and Adjacency =======================
# Similarity = (correlation + 1) / 2
# Adjacency = similarity ^ power
cat("\nCalculating adjacency matrix")
countsList <- readRDS(file.path(rdataDir, paste(prefix, "countsList.rds", sep = ".")))
adjList <- list()
for (set in names(countsList)) {
  cat("\n", rep("-", 10), "Working on", set, rep("-", 10), "\n")
  adjList[[set]] <- adjacency(
    datExpr = countsList[[set]],
    power = powerList[[set]],
    type = "signed"
  )
}
rm(countsList)
saveRDS(adjList, file = file.path(rdataDir, paste(prefix, "adjList.rds", sep = ".")))
adjList <- readRDS(file.path(rdataDir, paste(prefix, "adjList.rds", sep = ".")))

# Topological Overlap Matrix ===============
cat("\nCalculating TOM Matrix")
dissTOMList <- list()
for (set in names(adjList)) {
  cat("\n", rep("-", 10), "Working on", set, rep("-", 10), "\n")
  dissTOMList[[set]] <- 1 - TOMsimilarity(
    adjMat = adjList[[set]],
    TOMType = "signed"
  )
}
rm(adjList)
saveRDS(dissTOMList, file = file.path(rdataDir, paste(prefix, "dissTOMList.rds", sep = ".")))
dissTOMList <- readRDS(file.path(rdataDir, paste(prefix, "dissTOMList.rds", sep = ".")))

# Clustering ========================
cat("\nClustering")
DoClustering <- function(dissTOM, prefix) {
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  # Plot the resulting clustering tree (dendrogram)
  filename <- paste(prefix, "geneTree", "pdf", sep = ".")
  pdf(file.path(netPlotsDir, filename), width = 9, height = 6)
  plot(geneTree,
    xlab = "", sub = "",
    main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04
  )
  dev.off()
  return(geneTree)
}
geneTreeResList <- list()
for (set in names(dissTOMList)) {
  geneTreeResList[[set]] <- DoClustering(dissTOMList[[set]], set)
}
saveRDS(geneTreeResList, file = file.path(rdataDir, paste(prefix, "geneTreeResList.rds", sep = ".")))
geneTreeResList <- readRDS(file.path(rdataDir, paste(prefix, "geneTreeResList.rds", sep = ".")))

# Module identification using dynamic tree cut ========
# Dynamic tree cut to identify modules
cat("\nDynamic tree cut starting")
dissTOMList <- readRDS(file.path(rdataDir, paste(prefix, "dissTOMList.rds", sep = ".")))
dynModsList <- list()
for (set in names(geneTreeResList)) {
  dynModsList[[set]] <-
    cutreeDynamic(
      dendro = geneTreeResList[[set]],
      distM = dissTOMList[[set]],
      deepSplit = 2, # 4 is too sensitive, produces > 300 modules
      pamRespectsDendro = FALSE,
      minClusterSize = 30
    )
}
rm(dissTOMList)
saveRDS(dynModsList, file = file.path(rdataDir, paste(prefix, "dynModsList.rds", sep = ".")))
dynModsList <- readRDS(file.path(rdataDir, paste(prefix, "dynModsList.rds", sep = ".")))

# Convert numeric lables into colors ==================
cat("\nNum 2 label starting")
Num2Label <- function(dynMods, geneTree, prefix) {
  dynamicColors <- labels2colors(dynMods)
  # Plot the dendrogram and colors underneath
  filename <- paste(prefix, "colorDendro", "pdf", sep = ".")
  pdf(file.path(netPlotsDir, filename), width = 8, height = 6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
  dev.off()
  return(dynamicColors)
}

dynColorsList <- list()
for (set in names(dynModsList)) {
  dynColorsList[[set]] <- Num2Label(
    dynModsList[[set]], geneTreeResList[[set]], set
  )
}
rm(geneTreeResList, dynModsList)
saveRDS(dynColorsList, file = file.path(rdataDir, paste(prefix, "dynColorsList.rds", sep = ".")))
dynColorsList <- readRDS(file.path(rdataDir, paste(prefix, "dynColorsList.rds", sep = ".")))

length(unique(dynColorsList$N)) # 134
length(unique(dynColorsList$D)) # 117

# Calculate eigengenes ============================
cat("\nCalculating eigengenes")
CalcEigen <- function(expr, colours, MEDissThres, prefix) {
  MEList <- moduleEigengenes(expr, colors = colours)
  MEs <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1 - cor(MEs)
  # Cluster module eigengenes
  METree <- hclust(as.dist(MEDiss), method = "average")
  # Plot the result
  filename <- paste(prefix, "METree", "pdf", sep = ".")
  pdf(file.path(netPlotsDir, filename), width = 7, height = 6)
  plot(METree,
    main = "Clustering of module eigengenes",
    xlab = "", sub = ""
  )
  abline(h = MEDissThres, col = "red")
  dev.off()
  return(MEs)
}
countsList <- readRDS(file.path(rdataDir, paste(prefix, "countsList.rds", sep = ".")))
MEsList <- list()
for (set in names(dynColorsList)) {
  MEsList[[set]] <- CalcEigen(
    expr = countsList[[set]],
    colours = dynColorsList[[set]],
    MEDissThres = 0.3,
    prefix = set
  )
}
saveRDS(MEsList, file = file.path(rdataDir, paste(prefix, "MEsList.rds", sep = ".")))
rm(MEsList)

# Merge similar modules with automatic function  =======================
cat("\nMerging modules")
MergeMods <- function(expr, colours, MEDissThres) {
  merge <- mergeCloseModules(expr,
    colours,
    cutHeight = MEDissThres, verbose = 3
  )
  # The merged module colors
  mergedColors <- merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs <- merge$newMEs
  return(list("merge" = merge, "mergedColors" = mergedColors, "mergedMEs" = mergedMEs))
}
countsList <- readRDS(file.path(rdataDir, paste(prefix, "countsList.rds", sep = ".")))
mergedList <- list()
for (set in names(dynColorsList)) {
  cat("\n", rep("-", 10), "Working on", set, rep("-", 10), "\n")
  mergedList[[set]] <- MergeMods(
    expr = countsList[[set]],
    colours = dynColorsList[[set]],
    MEDissThres = 0.3
  )
}
rm(countsList)
saveRDS(mergedList, file = file.path(rdataDir, paste(prefix, "mergedList.rds", sep = ".")))
mergedList <- readRDS(file.path(rdataDir, paste(prefix, "mergedList.rds", sep = ".")))

ncol(mergedList$N$mergedMEs) # 82 modules after merging
ncol(mergedList$D$mergedMEs) # 77 modules after merging

formattedFreqList <- list()
for (set in names(mergedList)) {
  cat("\n", set, "\n")
  nonGreyColors <- mergedList[[set]]$mergedColors[
    !is.finite(match("grey", mergedList[[set]]$mergedColor))
  ]
  formattedFreq <- as.data.frame(table(nonGreyColors))$Freq
  formattedFreqList[[set]] <- formattedFreq
  print(summary(formattedFreq))
}

# D
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 33     130     282    1312    1185   20012
#
# N
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 35.0   101.8   224.0  1232.1   725.0 13193.0

# Plot final networks -------------------------
cat("\nFinal plot")
FinalNetPlot <- function(modColors, tree, prefix) {
  # Plot the dendrogram and colors underneath
  filename <- paste(prefix, "finalDendro", "pdf", sep = ".")
  pdf(file.path(netPlotsDir, filename), width = 8, height = 6)
  plotDendroAndColors(
    tree, modColors, c("Merged dynamic"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
  dev.off()
}
geneTreeResList <- readRDS(file.path(rdataDir, paste(prefix, "geneTreeResList.rds", sep = ".")))
for (set in names(geneTreeResList)) {
  FinalNetPlot(
    modColors = mergedList[[set]]$mergedColors,
    tree = geneTreeResList[[set]],
    prefix = set
  )
}
rm(geneTreeResList)

# Module membership -------------------------
cat("\nModule membership")
GetMM <- function(mergedObj, countDf) {
  # Labels to color
  ME_color_names <- names(mergedObj$mergedMEs)
  ME_colors <- substring(ME_color_names, 3)

  geneMM <- as.data.frame(cor(countDf, mergedObj$mergedMEs, use = "p"))
  pvalMM <- as.data.frame(corPvalueStudent(as.matrix(geneMM),
    nSamples = nrow(countDf)
  ))
  colnames(geneMM) <- paste("MM", ME_colors, sep = "")
  colnames(pvalMM) <- paste("p.MM", ME_colors, sep = "")
  return(list("geneMM" = geneMM, "pvalMM" = pvalMM))
}
countsList <- readRDS(file.path(rdataDir, paste(prefix, "countsList.rds", sep = ".")))
MMsList <- list()
for (set in names(mergedList)) {
  MMsList[[set]] <- GetMM(
    mergedObj = mergedList[[set]],
    countDf = countsList[[set]]
  )
}
saveRDS(MMsList, file = file.path(rdataDir, paste(prefix, "MMsList.rds", sep = ".")))
MMsList <- readRDS(file.path(rdataDir, paste(prefix, "MMsList.rds", sep = ".")))

# Hub genes -------------------------------------
cat("\nHub genes")
hubGenesList <- list()
for (set in names(mergedList)) {
  hubGenesList[[set]] <- chooseTopHubInEachModule(
    datExpr = countsList[[set]],
    colorh = mergedList[[set]]$mergedColors,
    power = powerList[[set]],
    type = "signed"
  )
}
saveRDS(hubGenesList, file = file.path(rdataDir, paste(prefix, "hubGenesList.rds", sep = ".")))
hubGenesList <- readRDS(file.path(rdataDir, paste(prefix, "hubGenesList.rds", sep = ".")))

# Make geneInfo ------------------------------
cat("\nMaking gene info")
# Function ========
MakeGeneInfo <- function(mergedObj, mm, countDf, hubGenes, deaRes, id2name) {
  mod_names <- mergedObj$mergedColors
  if (length(mod_names) == ncol(countDf)) {
    cat("\nNumber of mod names matches number of features")
  }
  # id to gene name
  mmInId2name <- match(rownames(mm$geneMM), id2name$isoform_id)
  geneNames <- id2name$gene_name[mmInId2name]
  # is a det
  mmInDea <- match(rownames(mm$geneMM), deaRes$isoform_id)
  isDet <- is.finite(mmInDea)
  # Hub genes
  isHub <- is.finite(match(rownames(mm$geneMM), hubGenes)) # hub gene info
  # Make gene info
  geneInfo0 <- data.frame(
    "isoform_id" = rownames(mm$geneMM), # isoform_id
    "gene_name" = geneNames, # gene name
    "is_hub" = isHub, # is a hub gene
    "is_det" = isDet, # is differentially expressed
    "module_color" = mod_names
  ) # module color
  # Add module membership
  for (idx in 1:ncol(mm$geneMM)) { # For every module
    old_colnames <- colnames(geneInfo0)
    # Append MM and MM p-value
    geneInfo0 <- data.frame(geneInfo0, mm$geneMM[, idx], mm$pvalMM[, idx])
    # The new columns' names are "MM." and "p.MM." + the module color
    colnames(geneInfo0) <- c(
      old_colnames,
      paste("MM.", colnames(mm$geneMM)[idx], sep = ""),
      paste("p.MM.", colnames(mm$geneMM)[idx], sep = "")
    )
  } # end loop through modules
  # Order by color
  geneInfo <- geneInfo0[order(geneInfo0$module_color), ]
  # Remove grey module
  geneInfo <- geneInfo %>% subset(module_color != "grey")
  # Return
  return(geneInfo)
} # end function

# Load metadata ========
id2name <- read.csv(id2namePath, row.names = 1, stringsAsFactors = FALSE) %>%
  mutate(gene_name = coalesce(gene_name, gene_id))
deaRes <- read.csv(deaResPath, row.names = 1, stringsAsFactors = FALSE)

# Execute ========
geneInfoList <- list()
for (set in names(mergedList)) {
  cat("\n", set)
  geneInfoList[[set]] <- MakeGeneInfo(
    mergedObj = mergedList[[set]],
    mm = MMsList[[set]],
    countDf = countsList[[set]],
    hubGenes = hubGenesList[[set]],
    deaRes = deaRes,
    id2name = id2name
  )
  write.csv(geneInfoList[[set]], file.path(geneInfoDir, paste(set, "csv", sep = ".")))
}

# End --------------------
cat("\nDone!")
