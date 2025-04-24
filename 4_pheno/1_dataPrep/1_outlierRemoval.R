#-------------------------------------------------
# Title: Data quality control and normalization
# Author: Amanda Zacharias
# Date: 2023-01-02
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
# Options -----------------------------------------
setwd("absolutePath/clbpHum")
projDir <- "4_pheno"
options(ggrepel.max.overlaps = Inf, box.padding = 0)
source("0_helpers/Normalize.R")
source("0_helpers/OutlierDetect.R")
source("0_helpers/ExtractOutliers.R")
source("0_helpers/filtering.R")
source("0_helpers/vstAndPCA.R")

# Packages -----------------------------------------
library(data.table) # 1.14.2
library(tidyverse) # 1.3.2
library(dplyr) # 1.1.0
library(stringr) # 1.4.0

library(Biobase) # ‘2.46.0’; provides 'AnnotatedDataFrame()' function
library(arrayQualityMetrics) # 3.42.0
library(DESeq2) # 1.26.0
library(edgeR) # 3.28.1

library(reshape2) # 1.4.4
library(ggplot2) # 3.3.6
library(ggrepel) # 0.9.1
library(cividis) # 0.2.0

# Pathways -----------------------------------------
# Input ===========
coldataPath <- file.path("coldata.csv")
countsPath <- file.path("3_stringtie", "prepDEcounts", "transcripts.csv")

# Output =========
baseDir <- file.path(projDir, "1_dataPrep")
qcB4OutDir <- file.path(baseDir, "qcB4Out")
qcAftOutDir <- file.path(baseDir, "qcAftOut")
cleanDataDir <- file.path(baseDir, "cleanData")
filtDataDir <- file.path(baseDir, "filtData")

plotsDir <- file.path(baseDir, "plots")
pcaDir <- file.path(plotsDir, "pca")
filtHistDir <- file.path(plotsDir, "filtHist")

system(paste(
  "mkdir", qcB4OutDir, qcAftOutDir, cleanDataDir,
  filtDataDir, plotsDir, pcaDir, filtHistDir
))

# Load data -----------------------------------------
# Coldata =========
coldata <- read.csv(coldataPath,
  row.names = 1,
  stringsAsFactors = F,
  check.names = FALSE
) %>%
  mutate(phenotype = as.factor(phenotype))
rownames(coldata) <- coldata$filename
# Remove technical duplicates
dupSamples <- c("3", "7", "4", "8")
coldata <- coldata[!coldata$filename %in% dupSamples, ]

# Count data =========
rawCounts <- fread(countsPath, stringsAsFactors = FALSE, check.names = FALSE) %>%
  tibble::column_to_rownames(var = "transcript_id") %>%
  dplyr::select(all_of(coldata$filename)) %>% # match order
  data.matrix()

# Raw QC -----------------------------------------
rawOutliers <- OutlierDetect(
  rawCounts, coldata, file.path(qcB4OutDir, "raw"), "phenotype"
)

# Normalize counts -----------------------------------------
normMat <- Normalize(rawCounts)
normCounts <- normMat %>% as.data.frame()

# Normalized QC -----------------------------------------
normOutliers <- OutlierDetect(
  normMat, coldata, file.path(qcB4OutDir, "norm"), "phenotype"
)

# Remove outliers round 1 -----------------------------------------
round1Outliers <- list("raw" = rawOutliers, "norm" = normOutliers)
round1ToRemove <- ExtractOutliers(round1Outliers)
capture.output(round1ToRemove, file = file.path(baseDir, "round1ToRemove.txt"))
saveRDS(round1ToRemove, file = file.path(baseDir, "round1ToRemove.rds"))

round1ToRemove <- readRDS(file.path(baseDir, "round1ToRemove.rds"))
# No outliers!

# Save files -----------------------------------------
# Raw counts
write.csv(
  rawCounts,
  file.path(cleanDataDir, paste("raw", "counts", "csv", sep = "."))
)
# Normalized counts
write.csv(
  normCounts,
  file.path(cleanDataDir, paste("norm", "counts", "csv", sep = "."))
)
write.csv(
  coldata,
  file.path(cleanDataDir, paste("coldata", "csv", sep = "."))
)

# Normalized & vst transformed counts ===========
rawVst <- GetVst(rawCounts, coldata, "phenotype")
cleanVst <- rawVst # there were no outliers !!!
write.csv(assay(cleanVst), file.path(
  cleanDataDir, paste("vst", "counts", "csv", sep = ".")
))

# PCA -----------------------------------------
intGrps <- c(
  "phenotype"
)

# Raw PCA ===========
rawLabels <- rep(FALSE, ncol(rawVst))
# rawLabels[match(round1ToRemove, rawVst$name)] <- TRUE
MakePCA(
  vstData = rawVst,
  intGroups = intGrps,
  labels = rawLabels,
  newFilename = "rawPCA.pdf",
  newPath = pcaDir
)

# Nonspecific filtering -----------------------------------------
# Plot ===========
GetMadCutoff(normCounts, plotsDir)

# Remove lowly expressed features ===========
# Remove using variance quantiles
mads <- GetMads(normCounts)
quantile(mads$mad, probs = seq(0.5, 0.70, 0.05))
#  50%         55%         60%         65%         70%
# 0.002189348 0.010463898 0.022160366 0.045823104 0.084102574
filtThreshold <- 0.084102574
toKeep <- rownames(subset(mads, mad >= filtThreshold))

rawFiltCounts <- data.frame(rawCounts[rownames(rawCounts) %in% toKeep, ], check.names = FALSE)
normFiltCounts <- data.frame(normCounts[rownames(normCounts) %in% toKeep, ], check.names = FALSE)
vstFiltCounts <- data.frame(assay(cleanVst)[rownames(assay(cleanVst)) %in% toKeep, ], check.names = FALSE)

# What are the differences? =======
cat(
  cat("Before filtering: ", nrow(normCounts), "\n"),
  cat("After filtering: ", nrow(normFiltCounts), "\n")
)
# Before filtering:  336781
# After filtering:  101035

# Save files ======
write.csv(
  rawFiltCounts,
  file.path(filtDataDir, paste("raw", "counts", "csv", sep = "."))
)
write.csv(
  normFiltCounts,
  file.path(filtDataDir, paste("norm", "counts", "csv", sep = "."))
)
write.csv(
  vstFiltCounts,
  file.path(filtDataDir, paste("vst", "counts", "csv", sep = "."))
)

# Plot histogram of dfs before and after filtering ----------------------------------------
MakeHist <- function(before, after, plot_title, newFilename, xtitle) {
  # Prepare plotting data
  plotData <- rbind(before %>% mutate(dataset = "before"), 
                    after %>% mutate(dataset = "after")) %>% 
    mutate(dataset = factor(dataset, levels = c("before", "after")))
  # Plotting
  gplot <- plotData %>% 
    ggplot(aes(x = log2(value + 1), fill = dataset)) +
    geom_histogram(data = subset(plotData, dataset == "before"), alpha = 0.75) + 
    geom_histogram(data = subset(plotData, dataset == "after"), alpha = 0.75) + 
    scale_fill_discrete(name = "Dataset") + 
    xlab(xtitle) +
    ylab("Frequency") +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25))
    )
  ggsave(
    plot = gplot, filename = newFilename, path = filtHistDir,
    width = 185, height = 185, units = "mm"
  )
  # return(gplot)
}
# Mads ========
MakeHist(
  data.frame(
    "id" = rownames(mads),
    "value" = mads$mad
  ),
  data.frame(
    "id" = rownames(normFiltCounts),
    "value" = mads[rownames(mads) %in% toKeep, ]
  ),
  str_to_sentence("Non-specific filtering"),
  paste("mad", "pdf", sep = "."),
  expression("Log"[2] * "(MAD of normalized expression + 1)")
)

# Save image ------------------------------------------------------------------
save.image(file = file.path(baseDir, "outlierRemoval.RData"))
