#-------------------------------------------------
# Title: Get alignment rates and descriptive statistics
# Author: Amanda Zacharias
# Date: 2022-12-17
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------
library(dplyr) # 1.0.9
library(tidyverse) # 1.3.2
library(reshape2) # 1.4.4
library(stringr) # 1.4.0
library(Hmisc) # 4.7.1
library(ggplot2) # 3.3.6
library(cividis) # 0.2.0

# Pathways -----------------------------------------
alignDir <- file.path(getwd(), "2_align")
ratesDir <- file.path(alignDir, "rates")

# Load data -----------------------------------------
# Get sample names / filenames
coldata <- read.csv(file.path("coldata.csv"),
  row.names = 1, stringsAsFactors = FALSE
)

# Get rates -----------------------------------------
GetRates <- function(alignDir) {
  # Prepare variables
  logPaths <- list.files(file.path(alignDir, "summaries"),
    full.names = TRUE
  )
  # Get rates
  filenames <- c()
  uniqueList <- c()
  multiList <- c()
  overallList <- c()
  for (idx in seq_along(logPaths)) {
    filepath <- logPaths[idx]
    filename <- gsub(".txt", "", basename(filepath))
    filenames <- c(filenames, filename)
    # uniqueue rate info on 4th line
    unique <- system(paste("sed '4q;d'", filepath, sep = " "), intern = TRUE)
    uniqueNums <- str_extract(string = unique, pattern = "(?<=\\().*(?=\\))")
    uniqueNums <- as.numeric(substr(uniqueNums, 1, nchar(uniqueNums) - 1)) # remove %
    uniqueList <- c(uniqueList, uniqueNums)
    # overall mapping rate on 6th line
    overall <- system(paste("sed '6q;d'", filepath, sep = " "), intern = TRUE)
    overallNums <- as.numeric(substr(overall, 1, 5))
    overallList <- c(overallList, overallNums)
  } # finish looping through samples
  # add to list
  rateDfs <- data.frame(
    filename = filenames,
    overall = overallList,
    unique = uniqueList
  )
  cat("Done!", "\n")
  return(rateDfs)
} # finish function
ratesDf <- GetRates(alignDir)

# What samples have particularly low alignment rates?
lowDf <- subset(ratesDf, overall < 70)
print(paste(lowDf$filename, lowDf$overall))
# character(0)

# Calculate statistics -----------------------------------------
GetStats <- function(ratesDf) {
  overallMin <- min(ratesDf$overall)
  overallMean <- mean(ratesDf$overall)
  overallMax <- max(ratesDf$overall)
  overallMedian <- median(ratesDf$overall)

  uniqueMin <- min(ratesDf$unique)
  uniqueMean <- mean(ratesDf$unique)
  uniqueMax <- max(ratesDf$unique)
  uniqueMedian <- median(ratesDf$unique)

  statsDf <- data.frame(
    groups = c("overall", "unique"),
    min = c(overallMin, uniqueMin),
    mean = c(overallMean, uniqueMean),
    max = c(overallMax, uniqueMax),
    median = c(overallMedian, uniqueMedian)
  )
  return(statsDf)
}
statsDf <- GetStats(ratesDf)

# Write dataframes -----------------------------------------
write.csv(ratesDf, file.path(ratesDir, "rates.csv"))
write.csv(statsDf, file.path(ratesDir, "stats.csv"))

# Visualize stats -----------------------------------------
# Boxplot =========
MakeBoxplot <- function(df, newfilename) {
  boxplot <- df %>%
    dplyr::select(-"filename") %>%
    stack() %>%
    dplyr::rename("value" = "values", "variable" = "ind") %>%
    ggplot(aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    geom_point(
      alpha = 0.5, shape = 21,
      show.legend = FALSE, position = position_jitterdodge()
    ) +
    scale_y_continuous(name = "Alignment rate", breaks = seq(0, 100, 10), limits = c(50, 100)) +
    scale_x_discrete(name = "Alignment rate type", labels = c("Overall", "Unique")) +
    scale_fill_cividis(
      name = "Alignment\nrate type", labels = c("Overall", "Unique"),
      begin = 0.2, end = 0.8, alpha = 0.5, discrete = TRUE
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = rel(1.1), hjust = 0.5),
      axis.title = element_text(size = rel(1.2)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.25),
      axis.ticks = element_line(linewidth = 1),
      axis.text = element_text(colour = "black")
    )

  pdf(file.path(ratesDir, newfilename),
    width = 4, height = 4
  )
  print(boxplot)
  dev.off()
}
MakeBoxplot(ratesDf, "allBoxplot_2025.pdf")
