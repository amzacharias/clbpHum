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
# Refer to section 3.4.2 in the edgeR user guide
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
library(UpSetR) # 1.4.0
library(ComplexHeatmap) # 2.2.0
library(extrafont) # 0.19

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
  mutate(
    # To factor !!!!
    phenotype = relevel(as.factor(phenotype), ref = "Rhythmic"),
    subjectId = as.factor(subjectId),
    batch = relevel(as.factor(batch), ref = 1)
  )
keepCounts <- counts %>% dplyr::select(all_of(keepColdata$filename))
keepVstCounts <- vstCounts %>% dplyr::select(all_of(keepColdata$filename))

varOfInterest <- "phenotype"
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
sex <- keepColdata$sex
subjectId <- keepColdata$subjectId
batch <- keepColdata$batch
timepoint <- as.factor(keepColdata$timepoint)
phenotype <- keepColdata$phenotype

design <- model.matrix(~ batch + timepoint + sex + phenotype)
colnames(design) <- gsub("phenotype", "", colnames(design))
cat("\nFinal design coefficients:\n", colnames(design), "\n")
#  (Intercept) batch2 timepointN sexM ConstantHigh ConstantLow Mixed
# Using "subjectId" in design creates linear dependencies between columns, making coefficients impossible to estimate.
# Repeated measures from the same person are being treated like technical replicates

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
coeffsToTest <- list(
  "RhythmicVsAll" = 5:7,
  "RhythmicVsConstant" = 5:6,
  "RhythmicVsCH" = 5,
  "RhythmicVsCL" = 6,
  "RhythmicVsMixed" = 7
)
qlts <- list()
for (coeffName in names(coeffsToTest)) {
  cat("\n", coeffName)
  qlts[[coeffName]] <- glmQLFTest(fit, coef = coeffsToTest[[coeffName]])
}

# Results -----------------------------------------
# What genes are different between rhythmic and other phenotype(s)? ======
resList <- list()
for (coeff in names(qlts)) {
  resList[[coeff]] <-
    GetResults(
      fitObj = qlts[[coeff]],
      nGenesInput = nrow(keepCounts),
      id2nameObj = id2name,
      idColumn = "isoform_id"
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
load(file.path(rDataDir, paste(opt$basename, "RData", sep = ".")))
# Volcano =======
cat("\nMaking volcano plot\n")
titles <- list(
  "RhythmicVsCH" = "rhythmic versus constant-h",
  "RhythmicVsCL" = "rhythmic versus constant-l",
  "RhythmicVsMixed" = "rhythmic versus mixed"
)
colors <- list(
  "RhythmicVsCH" = "#F09F39",
  "RhythmicVsCL" = "#731B0D",
  "RhythmicVsMixed" = "#8A51C2"
)
maxLim <- max(
  sapply(
    3:length(resList),
    function(n) max(-log2(resList[[n]]$topgenesIdName$PValue))
  )
) + 2
for (coeff in names(resList)) {
  if (sum(grepl("logFC", colnames(resList[[coeff]]$topgenesIdName))) == 1) {
    cat("\n", coeff)
    # If there is only one logFC column
    # PDF
    MakeVolcano(
      resList[[coeff]]$topgenesIdName, coeff,
      showName = FALSE,
      newTitle = titles[[coeff]],
      newPath = basePlotsDir,
      sigColor = colors[[coeff]],
      ylims = c(0, maxLim),
      idColumn = "isoform_id",
      newHeight = 150
    )
    MakeVolcano(
      resList[[coeff]]$topgenesIdName, coeff,
      showName = TRUE,
      newTitle = titles[[coeff]],
      newPath = basePlotsDir,
      sigColor = colors[[coeff]],
      ylims = c(0, maxLim),
      idColumn = "isoform_id",
      newHeight = 150
    )
  }
}

# Save image -----------------------------------------
save.image(file.path(rDataDir, paste(opt$basename, "RData", sep = ".")))
