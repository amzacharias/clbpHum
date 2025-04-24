#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Module trait relationship
# Author: Amanda Zacharias
# Date: 2023-04-10
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
prefix <- "transcriptSexDayNight"
projDir <- "4_pheno"
setwd("absolutePath/clbpHum")
source("0_helpers/GetNamedPaths.R")
source("0_helpers/ReadDfs.R")

# Packages -----------------------------------------
library(dplyr) # 1.1.0
library(WGCNA) # 1.70.3
library(ggplot2) # 3.4.1
library(broom) # 1.0.3

# Pathways -----------------------------------------
# Input paths ======
dataPrepDir <- file.path(projDir, "1_dataPrep")
coldataPath <- file.path(dataPrepDir, "cleanData", "coldata.csv")

wgcnaDir <- file.path(projDir, "3_wgcna", prefix)
geneInfoDir <- file.path(wgcnaDir, "geneInfo")
geneInfoPaths <- GetNamedPaths(file.path(wgcnaDir, "geneInfo"), ".csv")
rdataDir <- file.path(wgcnaDir, "rdata")
MEsPath <- file.path(rdataDir, paste(prefix, "mergedList.rds", sep = "."))

# Output paths =====
plotsDir <- file.path(wgcnaDir, "plots")
modTraitDir <- file.path(wgcnaDir, "modTrait")
modTraitPlotsDir <- file.path(plotsDir, "modTrait")
indivRelationPlotDir <- file.path(modTraitPlotsDir, "indivRelation")

# Make output directories
system(paste("mkdir", modTraitDir, plotsDir, modTraitPlotsDir, indivRelationPlotDir))

# Load data -----------------------------------------
# Coldata ======
coldata <- read.csv(coldataPath, row.names = 1, stringsAsFactors = TRUE)
coldata <- coldata %>%
  mutate(
    "phenotype" = factor(coldata$phenotype, levels = unique(coldata$phenotype))
  )
coldataD <- coldata %>% subset(timepoint == "D")
coldataN <- coldata %>% subset(timepoint == "N")
rm(coldata)

# WGCNA data ========
geneInfoList <- ReadDfs(geneInfoPaths)
mergedMEs <- readRDS(MEsPath)
MEsD <- mergedMEs$D$mergedMEs
MEsN <- mergedMEs$N$mergedMEs

# Add yes no Rhythmic column -----------------------------------------
# Rhythmic = 1; other = 0
coldataD <- coldataD %>%
  mutate("isRhythmic" = if_else(phenotype == "Rhythmic", 1, 0)) %>%
  mutate_at("isRhythmic", as.factor)
coldataN <- coldataN %>%
  mutate("isRhythmic" = if_else(phenotype == "Rhythmic", 1, 0)) %>%
  mutate_at("isRhythmic", as.factor)

# Match order of samples in coldata and MEs -----------------------------------------
MEsD <- MEsD[rownames(MEsD) == coldataD$filename, ]
MEsN <- MEsN[rownames(MEsN) == coldataN$filename, ]

# Plot relationships -----------------------------------------
PlotRelation <- function(modelData, response = "phenotype", time = "D") {
  #' Plot module eigengene counts across the variable of interest as violin plot
  #'
  #' @param modelData results from an association test as dataframe
  #' @param response The variable of interest to plot
  #' @param time The variable by which data is stratified, impacts output filename
  #' @return Returns nothing, writes a pdf file.
  #' @example PlotRelation(modelData, "phenotype", "D")
  # Make folder for output files if folder doesn't already exist
  outDir <- file.path(indivRelationPlotDir, response)
  if (file.exists(outDir) == FALSE) {
    system(paste("mkdir", outDir))
  }
  # Create plot
  modName <- colnames(modelData)[2]
  gplot <- modelData %>%
    ggplot(aes(x = .data[[response]], y = .data[[modName]], color = .data[[response]])) +
    geom_violin() +
    geom_jitter(alpha = 0.5) +
    coord_flip() +
    ggtitle(modName) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  # Save plot
  ggsave(
    plot = gplot, path = outDir,
    filename = paste(time, response, modName, "pdf", sep = "."),
    width = 138.75, height = 92.5, units = "mm"
  )
}

# Binomial / Logistic regression -----------------------------------------
RunGlm <- function(MEs, coldat, geneInfo, formulaStart,
                   formulaEnd = "", response = "base", time = "D") {
  #'Use a generalized linear model to find association between module eigengene and response variable
  #'
  #'@param MEs dataframe with module eigengenes' counts across samples
  #'@param coldat dataframe with metadata about samples
  #'@param geneInfo dataframe with information about genes in modules
  #'@param formulaStart the beginning of the glm forumula, a string
  #'@param formulaEnd the end of the glm formula, a string
  #'@param response The name of the response variable of interest
  #'@param time The variable by which data is stratified, impacts output file names
  #'@return Returns a named list with the glm results object and
  #'        a dataframe that summarizes results
  #'@example res <- RunGlm(MEs, coldat, geneInfo, "pheno ~", "", "pheno", "D")
  fitsList <- list()
  initDf <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(initDf) <- c(
    "term", "modSize", "estimate", "std.error", "statistic",
    "p.value", "conf.low", "conf.high",
    "OR", "OR.low", "OR.high", "warning"
  )
  for (idx in seq_len(MEs)) { # for every module
    # Variables
    modName <- colnames(MEs)[idx]
    modSize <- nrow(geneInfo %>% subset(module_color == gsub("ME", "", modName)))
    # Extract needed data
    modelData <- MEs %>%
      tibble::rownames_to_column("filename") %>%
      dplyr::select(all_of(c("filename", modName))) %>%
      merge(coldat)
    formula <- paste(formulaStart, modName, formulaEnd)
    w <- length(warnings()) # get warnings length
    # Run analysis
    withCallingHandlers(
      {
        # Fit model
        fit <- glm(as.formula(formula), data = modelData, family = binomial())
        # Convert to tabular
        fitTabular <- tidy(fit, conf.int = TRUE) # gets a two-sided p-value
        fitClean <- fitTabular %>%
          # Remove intercept row
          filter(term != "(Intercept)") %>%
          # Add odds ratios
          mutate(
            OR = exp(estimate),
            OR.conf.low = exp(conf.low),
            OR.conf.high = exp(conf.high)
          ) %>% # Add module size
          mutate("modSize" = modSize, .before = "estimate")
      },
      warning = function(w) {
        cat("\nWarning for", modName, w$message)
        initDf[idx, ]$warning <<- w$message
        invokeRestart("muffleWarning")
      }
    ) # end warning handling
    if (fitClean$p.value < 0.05) {
      # If nominally significant, plot
      PlotRelation(modelData, response = response, time = time)
    } # end if pval if else statement
    # Add to output objects
    fitsList[[modName]] <- fit
    initDf[idx, 1:11] <- fitClean
  } # end loop through MEs
  # Finalize output dataframe
  allRes <- initDf %>%
    mutate("response" = response, .before = "term")
  cat("\nNumber nominal", nrow(allRes %>% subset(p.value < 0.05)))
  # Write results
  write.csv(
    allRes,
    file.path(modTraitDir, paste(time, response, "csv", sep = "."))
  )
  # Return
  return(list("fits" = fitsList, "allRes" = allRes))
}
isRhythmic ~ MEs
levels(coldataD$isRhythmic) # "0", "1"
resIsCircD <- RunGlm(
  MEs = MEsD,
  coldat = coldataD,
  geneInfo = geneInfoList$D,
  formulaStart = "isRhythmic ~",
  formulaEnd = "",
  response = "isRhythmic"
)
# Number sig 3
levels(coldataN$isRhythmic) # "0", "1"
resIsCircN <- RunGlm(
  MEs = MEsN,
  coldat = coldataN,
  geneInfo = geneInfoList$N,
  formulaStart = "isRhythmic ~",
  formulaEnd = "",
  response = "isRhythmic"
)
# Number sig 5

# Save and load results, so don't have to re-run analyses  -----------------------------------------
# Save ======
saveRDS(resIsCircD, file = file.path(modTraitDir, "resIsCircD.rds"))
saveRDS(resIsCircN, file = file.path(modTraitDir, "resIsCircN.rds"))
# Load =======
resIsCircD <- readRDS(file.path(modTraitDir, "resIsCircD.rds"))
resIsCircN <- readRDS(file.path(modTraitDir, "resIsCircN.rds"))

# Forest plots  -----------------------------------------
# Prepare dataframes for plotting =======
plotResD <- resIsCircD$allRes %>%
  mutate(
    text = ifelse(p.value < 0.05, p.value, NA),
    term = sub("ME", "", term)
  )
plotResN <- resIsCircN$allRes %>%
  mutate(
    text = ifelse(p.value < 0.05, p.value, NA),
    term = sub("ME", "", term)
  )

# Forest plot function =======
MakeForest <- function(plotDf, xVar, minVar, maxVar, xLabel, time, newTitle = "",
                       newFilename = "glm.pdf", newWidth = 185, newHeight = 250) {
  # Intercept =======
  if (xVar == "estimate") {
    xIntercept <- 0
  } else {
    xIntercept <- 1
  }
  # Make plot ======
  allPlot <- plotDf %>%
    ggplot(aes(y = term, x = .data[[xVar]])) +
    geom_vline(xintercept = xIntercept, linetype = "dashed") +
    geom_point() +
    geom_errorbarh(
      aes(
        xmin = .data[[minVar]], xmax = .data[[maxVar]],
        height = 0
      ),
      linewidth = 0.6
    ) +
    geom_text(aes(label = signif(text, 2), vjust = -0.5), show.legend = FALSE) +
    xlab(label = xLabel) +
    ylab(label = "Module name") +
    ggtitle(newTitle) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, -10, 0, 10),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(
    plot = allPlot, filename = paste(time, newFilename, sep = "."),
    path = modTraitPlotsDir,
    width = newWidth, height = newHeight, units = "mm"
  )
}
# Execute for all results ======
MakeForest(
  plotResD, "estimate", "conf.low", "conf.high", expression(beta), "D",
  newFilename = "glmRhythm.pdf", newTitle = "Day",
  newWidth = 92
)
MakeForest(
  plotResN, "estimate", "conf.low", "conf.high", expression(beta), "N",
  newFilename = "glmRhythm.pdf", newTitle = "Night",
  newWidth = 92
)
# Execute for only sig results ======
MakeForest(
  plotResD %>% subset(p.value < 0.05),
  "estimate", "conf.low", "conf.high", expression(beta), "D",
  newFilename = "glmSigRhythm.pdf", newTitle = "Day",
  newWidth = 92, newHeight = 75
)
MakeForest(
  plotResN %>% subset(p.value < 0.05),
  "estimate", "conf.low", "conf.high", expression(beta), "N",
  newFilename = "glmSigRhythm.pdf", newTitle = "Night",
  newWidth = 92, newHeight = 75
)
