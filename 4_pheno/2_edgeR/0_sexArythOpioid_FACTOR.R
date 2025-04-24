#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Differential expression analysis with edgeR
# Author: Amanda Zacharias
# Date: 2023-03-01
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# Perform single contrasts between phenotypes
# Refer to sections 3.3 - 3.5 in the edgeR user guide
# http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# Options -----------------------------------------
setwd("absolutePath/clbpHum")
projDir <- "4_pheno"
options(ggrepel.max.overlaps = Inf, box.padding = 0)

# Packages -----------------------------------------
library(optparse) # 1.7.3
library(dplyr) # 1.1.0
library(tibble) # 3.1.8
library(edgeR) # 3.28.1
library(ggplot2) # 3.4.1
library(ggrepel) # 0.9.3
library(ComplexHeatmap) # 2.2.0

# Source -----------------------------------------------
source(file.path("0_helpers", "GetResults.R"))
source(file.path("0_helpers", "MakeUpset.R"))
source(file.path("0_helpers", "MakeVolcano.R"))

# Run optparser -----------------------------------------
# Define arguments ==================
optionList <- list(
  make_option(c("-n", "--basename",
    type = "character", default = NA,
    help = "the base filename for this analysis' outputs", metavar = "character"
  )),
  make_option(c("-p", "--countPath"),
    type = "character", default = NA,
    help = "path to your count matrix", metavar = "character"
  ),
  make_option(c("-v", "--vstCountPath"),
    type = "character", default = NA,
    help = "path to your log transformed count matrix", metavar = "character"
  ),
  make_option(c("-w", "--workingDir", type = "character", default = getwd()),
    help = "the project working directory path", metavar = "character"
  )
)
# Get parameters ==================
optParser <- OptionParser(option_list = optionList)
opt <- parse_args(optParser)

# There is no parameter checking!!!!!
cat(
  "\nFormula:", opt$formula,
  "\nBase filename:", opt$basename, "\n"
)

# Set the working directory -----------------------------------------
setwd(opt$workingDir)

# Pathways -----------------------------------------
# Inputs =====
coldataPath <- file.path(projDir, "1_dataPrep", "cleanData", "coldata.csv")
id2namePath <- file.path(projDir, "1_dataPrep", "id2name.csv")

# Outputs =====
baseDir <- file.path(projDir, "2_edgeR")
dfsDir <- file.path(baseDir, "dataframes")
baseDfsDir <- file.path(dfsDir, opt$basename)
rDataDir <- file.path(baseDir, "rData")
plotsDir <- file.path(baseDir, "plots")
basePlotsDir <- file.path(plotsDir, opt$basename)
system(paste(
  "mkdir", dfsDir, baseDfsDir, rDataDir,
  plotsDir, basePlotsDir
))

# Load data -----------------------------------------
coldata <- read.csv(coldataPath, row.names = 1, stringsAsFactors = FALSE)
id2name <- read.csv(id2namePath, row.names = 1, stringsAsFactors = FALSE) %>%
  mutate(gene_name = coalesce(gene_name, gene_id))
counts <- read.csv(opt$countPath,
  row.names = 1,
  check.names = FALSE, stringsAsFactors = FALSE
)
vstCounts <- read.csv(opt$vstCountPath,
  row.names = 1,
  check.names = FALSE, stringsAsFactors = FALSE
)

# Execute edgeR ####################################
# Subset samples by those included in analysis -----------------------------------------
keepColdata <- coldata %>%
  subset(phenotype != "Rhythmic") %>%
  mutate(
    subjectId = as.factor(subjectId),
    batch = relevel(as.factor(batch), ref = 1),
    timepoint = relevel(as.factor(timepoint), ref = "D"),
    SMH_med_opioids_refactored = factor(
      ifelse(SMH_med_opioids == 2, 1, SMH_med_opioids), levels = c(0, 1)
    )
  )
keepCounts <- counts %>% dplyr::select(all_of(keepColdata$filename))
keepVstCounts <- vstCounts %>% dplyr::select(all_of(keepColdata$filename))

varOfInterest <- "SMH_med_opioids_refactored"
cat("\nLevels in variable of interest:\n")
levels(keepColdata[[varOfInterest]])

# Descriptive stats
cat("\nNumber of people in var of interest")
keepColdata %>%
  group_by_at(varOfInterest) %>%
  tally()

# Prepare objects -----------------------------------------
dge <- DGEList(counts = keepCounts)
# Base design
design <- model.matrix(~ batch + timepoint + sex, data = keepColdata)
# Make specific contrasts
opioidCurrent <- keepColdata$SMH_med_opioids_refactored == 1
# Add to design
design <- cbind(
  design, opioidCurrent
)
cat("\nFinal design coefficients:\n", colnames(design), "\n")
# (Intercept) batch2 timepointN sexM opioidCurrent

# Esimate dispersion -----------------------------------------
dge <- estimateDisp(dge, design, robust = TRUE)
pdf(
  file = file.path(basePlotsDir, "bcv.pdf"),
  width = 6, height = 6
)
plotBCV(dge)
dev.off()

# GLM -------
# Quasi-likelyhood tests are recommended when 3+ groups and
# when using bulk RNA seq with replicates
# Fit ========
fit <- glmQLFit(dge, design, robust = TRUE) # Fit model
pdf(
  file = file.path(basePlotsDir, "QLDisp.pdf"),
  width = 6, height = 6
)
plotQLDisp(fit)
dev.off()

cat("\nTerms in the model:\n")
colnames(fit)

# Test =========
# What genes change between any of the groups?
coefficients <- c(
  "opioidCurrent"
)
qlts <- list()
for (coeff in coefficients) {
  cat("\n", coeff)
  qlts[[coeff]] <- glmQLFTest(fit, coef = coeff)
}

# Results -----------------------------------------
# What genes are different opioid vs non opioid? ======
resList <- list()
for (coeff in names(qlts)) {
  resList[[coeff]] <-
    GetResults(
      fitObj = qlts[[coeff]],
      nGenesInput = nrow(keepCounts),
      id2nameObj = id2name
    )
  cat(
    "\nNumber of sig for", coeff, ":\n",
    nrow(resList[[coeff]]$sigTable)
  )
}

# Save results -----------------------------------------
# Each group
for (coeff in names(resList)) {
  write.csv(
    resList[[coeff]]$topgenesIdName,
    file.path(baseDfsDir, paste(coeff, "allRes.csv", sep = "_"))
  )
  write.csv(
    resList[[coeff]]$sigTable,
    file.path(baseDfsDir, paste(coeff, "sigRes.csv", sep = "_"))
  )
}

save.image(file.path(rDataDir, paste(opt$basename, "RData", sep = ".")))

# Plots -----------------------------------------
# Volcano =======
cat("\nMaking volcano plot\n")
colors <- list(
  "opioidCurrent" = "red"
)
maxLim <- max(
  sapply(
    seq_along(resList),
    function(n) max(-log2(resList[[n]]$topgenesIdName$PValue))
  )
) + 2 # increase from +1 so top gene isn't removed
for (coeff in names(resList)) {
  if (sum(grepl("logFC", colnames(resList[[coeff]]$topgenesIdName))) == 1) {
    cat("\n", coeff)
    # If there is only one logFC column
    title <- "opioid use"
    MakeVolcano(
      resList[[coeff]]$topgenesIdName, coeff,
      showName = FALSE,
      newTitle = title,
      newPath = basePlotsDir,
      sigColor = colors[[coeff]],
      ylims = c(0, maxLim),
      newHeight = 150
    )
    MakeVolcano(
      resList[[coeff]]$topgenesIdName, coeff,
      showName = TRUE,
      newTitle = title,
      newPath = basePlotsDir,
      sigColor = colors[[coeff]],
      ylims = c(0, maxLim),
      newHeight = 150
    )
  }
}

# Save image -----------------------------------------
save.image(file.path(rDataDir, paste(opt$basename, "RData", sep = ".")))
