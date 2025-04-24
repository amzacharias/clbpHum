#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Module trait relationship - outcome is categorical with multiple levels
# Author: Amanda Zacharias
# Date: 2023-05-08
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
library(tidyverse) # 1.3.2
library(WGCNA) # 1.70.3
library(nnet) # 7.3.18
library(broom) # 1.0.3
library(ggplot2) # 3.4.1

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
indivCoeffPlotDir <- file.path(modTraitPlotsDir, "indivCoeff")

# Make output directories
system(paste(
  "mkdir", modTraitDir, plotsDir,
  modTraitPlotsDir, indivRelationPlotDir, indivCoeffPlotDir
))

# Load data -----------------------------------------
# Coldata ======
coldata <- read.csv(coldataPath, row.names = 1, stringsAsFactors = TRUE)
coldata <- coldata %>%
  mutate(
    "phenotype" = factor(coldata$phenotype, levels = unique(coldata$phenotype))
  )
coldataD <- coldata %>%
  subset(timepoint == "D") %>%
  mutate("phenotype" = relevel(phenotype, ref = "Rhythmic"))
coldataN <- coldata %>%
  subset(timepoint == "N") %>%
  mutate("phenotype" = relevel(phenotype, ref = "Rhythmic"))
rm(coldata)

# WGCNA data ========
geneInfoList <- ReadDfs(geneInfoPaths)
mergedMEs <- readRDS(MEsPath)
MEsD <- mergedMEs$D$mergedMEs
MEsN <- mergedMEs$N$mergedMEs

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

# Plot coefficients -----------------------------------------
PlotCoeff <- function(fitTabular, outcome = "phenotype", time = "D") {
  #' Plot module eigengene counts across the variable of interest as violin plot
  #'
  #' @param fitTabular results from an association test as dataframe
  #' @param outcome The variable of interest to plot
  #' @param time The variable by which data is stratified, impacts output filename
  #' @return Returns nothing, writes a pdf file.
  #' @example PlotCoeff(fitTabular, "phenotype", "D")
  # Make folder for output files if folder doesn't already exist
  outDir <- file.path(indivCoeffPlotDir, outcome)
  if (file.exists(outDir) == FALSE) {
    system(paste("mkdir", outDir))
  }
  modName <- unique(fitTabular$term)
  plotDf <- fitTabular %>% mutate(text = ifelse(fitTabular$p.value < 0.05, fitTabular$p.value, NA))
  gplot <- plotDf %>%
    ggplot(aes(
      y = estimate, x = y.level,
      ymin = conf.low, ymax = conf.high
    )) +
    geom_pointrange() +
    geom_text(aes(label = signif(text, 2), vjust = -0.5), size = 5, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(label = "Regression coefficients") +
    xlab(label = "Outcome") +
    ggtitle(sub("ME", "", modName)) +
    coord_flip() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(
    plot = gplot, path = outDir,
    filename = paste(time, outcome, modName, "pdf", sep = "."),
    width = 138.75, height = 92.5, units = "mm"
  )
}

# Run multinomial logistic regression -----------------------------------------
RunMultinom <- function(MEs, coldat, geneInfo, formulaStart,
                        formulaEnd = "", response = "base", time = "D") {
  #' Find association with module eigengen counts and a response variable
  #' @param MEs dataframe with module eigengenes' counts across samples
  #' @param coldat dataframe with metadata about samples
  #' @param geneInfo dataframe with information about genes in modules
  #' @param formulaStart the beginning of the glm forumula, a string
  #' @param formulaEnd the end of the glm formula, a string
  #' @param response The name of the response variable of interest
  #' @param time The variable by which data is stratified, impacts output file names
  #' @param type The type of glm used. Options are binomial or poisson
  #' @return Returns a named list with the multinom results object and
  #'        a dataframe that summarizes results
  #' @example res <- RunMultinom(MEs, coldat, geneInfo, "pheno ~", "", "pheno", "D")
  fitsList <- list()
  # Instantiate a dataframe with results
  initDf <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(initDf) <- c(
    "y.level", "term", "modSize", "estimate", "std.error", "statistic",
    "p.value", "conf.low", "conf.high", 
    "RR", "RR.conf.low", "RR.conf.high", "warning"
  )
  for (idx in seq_len(MEs)) { # for every module
    # Variables
    initDfIdx <- nrow(initDf) + 1
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
        fit <- multinom(as.formula(formula), data = modelData, trace = FALSE)
        # Convert to tabular
        fitTabular <- tidy(fit, conf.int = TRUE) # gets a two-sided p-value
        fitClean <- fitTabular %>% 
          filter(term != "(Intercept)") %>% 
          mutate(
            RR = exp(estimate),
            RR.conf.low = exp(conf.low),
            RR.conf.high = exp(conf.high)
          ) %>%
          mutate("modSize" = modSize, .before = "estimate")
      },
      warning = function(w) {
        cat("\nWarning for", modName, w$message)
        initDf[initDfIdx:(initDfIdx + 2), ]$warning <<- w$message
        invokeRestart("muffleWarning")
      }
    ) # end warning handling
    # Plot
    if (sum(fitClean$p.value < 0.05) >= 1) {
      # Lines plot
      PlotRelation(modelData, "phenotype", time)
      # Plot coefficients
      PlotCoeff(fitClean, "phenotype", time)
    }
    # Save to output objects
    fitsList[[modName]] <- fit
    initDf[initDfIdx:(initDfIdx + 2), 1:12] <- fitClean
  }
  # Finalize output dataframe
  allRes <- initDf %>%
    mutate("response" = response, .before = "term")
  # Save to csv
  write.csv(
    allRes,
    file.path(modTraitDir, paste(time, response, "csv", sep = "."))
  )
  # Return fits and dataframe
  return(list("fits" = fitsList, "allRes" = allRes))
}
resChronoD <- RunMultinom(
  MEsD,
  coldataD,
  geneInfoList$D,
  formulaStart = "phenotype ~",
  response = "phenotypeNN",
  time = "D"
)
resChronoN <- RunMultinom(
  MEsN,
  coldataN,
  geneInfoList$N,
  formulaStart = "phenotype ~",
  response = "phenotypeNN",
  time = "N"
)
# Save and load results -----------------------------------------
# Save ========
saveRDS(resChronoD, file = file.path(modTraitDir, "mlrResD.rds"))
saveRDS(resChronoN, file = file.path(modTraitDir, "mlrResN.rds"))
# Load =======
resChronoD <- readRDS(file.path(modTraitDir, "mlrResD.rds"))
resChronoN <- readRDS(file.path(modTraitDir, "mlrResN.rds"))

# Forest plotting -----------------------------------------
# Prepare dataframes ========
plotResD <- resChronoD$allRes %>%
  mutate(
    text = ifelse(p.value < 0.05, p.value, NA),
    term = sub("ME", "", term),
    y.level = factor(y.level, levels = c("ConstantLow", "Mixed", "ConstantHigh"))
  )
plotResN <- resChronoN$allRes %>%
  mutate(
    text = ifelse(p.value < 0.05, p.value, NA),
    term = sub("ME", "", term),
    y.level = factor(y.level, levels = c("ConstantLow", "Mixed", "ConstantHigh"))
  )

# Function to make plot ========
MakeForest <- function(plotDf, xVar, minVar, maxVar, xLabel, time, newFilename = "multinom.pdf") {
  # Intercept =======
  if (xVar == "estimate") {
    xIntercept <- 0
  } else {
    xIntercept <- 1
  }
  # Labels =====
  lbls <- c("constant (l)", "mixed", "constant (h)")
  names(lbls) <- c("ConstantLow", "Mixed", "ConstantHigh")
  # Colors =====
  colors <- c("#731B0D", "#934DC9", "#F09F39")
  names(colors) <- c("ConstantLow", "Mixed", "ConstantHigh")
  # Shapes ======
  shapes <- c(16, 18, 15)
  names(shapes) <- c("ConstantLow", "Mixed", "ConstantHigh")
  # Make plot ======
  allPlot <- plotDf %>%
    ggplot(aes(y = term, x = .data[[xVar]])) +
    geom_vline(xintercept = xIntercept, linetype = "dashed") +
    geom_point(aes(color = y.level, shape = y.level)) +
    geom_errorbarh(
      aes(xmin = .data[[minVar]], xmax = .data[[maxVar]],
          height = 0, color = y.level
          ),
      linewidth = 0.6
    ) +
    geom_text(aes(label = signif(text, 2), vjust = -0.5), show.legend = FALSE) +
    xlab(label = xLabel) +
    ylab(label = "Outcome") +
    scale_color_manual(name = "Phenotype", labels = lbls, values = colors) +
    scale_shape_manual(name = "Phenotype", labels = lbls, values = shapes) +
    facet_grid(~y.level, scales = "free", labeller = labeller(lbls)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, -10, 0, 10)
    )
  ggsave(
    plot = allPlot, filename = paste(time, newFilename, sep = "."),
    path = modTraitPlotsDir,
    width = 185, height = 250, units = "mm"
  )
}
# Execute for all results ==========
MakeForest(
  plotResD, "estimate", "conf.low", "conf.high",
  expression(paste("Regression coefficient (", beta, ")")), "D", "coeff.multinom.pdf"
)
MakeForest(
  plotResD, "RR", "RR.conf.low", "RR.conf.high",
  "Risk Ratio", "D", "RR.multinom.pdf"
)
MakeForest(
  plotResN, "estimate", "conf.low", "conf.high",
  expression(paste("Regression coefficient (", beta, ")")), "N", "coeff.multinom.pdf"
)
MakeForest(
  plotResN, "RR", "RR.conf.low", "RR.conf.high",
  "Risk Ratio", "N", "RR.multinom.pdf"
)
