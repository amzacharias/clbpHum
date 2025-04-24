#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make volcano plot helper
# Author: Amanda Zacharias
# Date: 2023-04-09
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------
options(ggrepel.max.overlaps = Inf, box.padding = 0)

# Packages -----------------------------------------
library(ggplot2) # 3.4.1
library(ggrepel) # 0.9.3
# library(extrafont) # 0.19
library(stringr) # 1.5.0

# Pathways -----------------------------------------



# Load data -----------------------------------------
# font_import(paths = file.path("0_resources", "TimesNewRoman"), prompt = FALSE)
# Only need to run above once in an environment

# Function -----------------------------------------
MakeVolcano <- function(resDf, prefix, showName = TRUE,
                        newTitle = "Volcano", newPath = ".", fileEnd = "pdf",
                        sigColor = "#FF0000", ylims = c(0, 50),
                        idColumn = "isoform_id", newHeight=185) {
  #' Make a volcano plot
  #'
  # Copy results table to new variable that we can alter
  volcanoDf <- resDf
  # Mark significant rows
  volcanoDf$sig <- volcanoDf$FWER < 0.05
  # Order dataframe and get top 20 sig transcripts'
  if (showName == TRUE) {
    num2label <- 15
  } else if (showName == FALSE) {
    # transcript ids take up more space, so label less
    num2label <- 10
  }
  resOrdered <- volcanoDf %>% arrange(FWER)
  resOrdered$labels <- rownames(resOrdered) %in%
    rownames(resOrdered[1:num2label, ])

  # Plotting
  bcThreshold <- 0.05 / nrow(resOrdered)
  cat("\nBonferroni threshold is:", bcThreshold, "\n")

  volcano <- resOrdered %>%
    ggplot(aes(x = logFC, y = -log2(PValue))) +
    geom_point(aes(colour = sig), size = 2, alpha = 0.6, show.legend = FALSE) +
    ggtitle(str_to_sentence(newTitle)) +
    xlab(expression("Log"[2] * "(fold change)")) +
    ylab(expression("-Log"[2] * " (p-value)")) +
    ylim(ylims) +
    geom_hline(aes(yintercept = -log2(bcThreshold), linetype = "Adjusted"), colour = "red") +
    geom_hline(aes(yintercept = -log2(0.05), linetype = "Nominal")) +
    scale_linetype_manual(name = "Significance", values = c(2, 4)) +
    scale_color_manual(
      values = c("#000000", sigColor),
      labels = c("Not significant", "Significant")
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = rel(1.1), hjust = 0.5),
      axis.title = element_text(size = rel(1.2))
    )
  # Labels
  if (showName == TRUE) {
    volcano <- volcano + geom_label_repel(
      size = 4, nudge_y = 2, min.segment.length = 0.25,
      force = 12, show.legend = FALSE,
      fontface = "italic",
      # family = "Times",
      aes(label = ifelse(labels == TRUE, # label top t's.
        as.character(gene_name), ""
      ))
    )
  } else if (showName == FALSE) {
    volcano <- volcano + geom_label_repel(
      size = 4, nudge_y = 2, min.segment.length = 0.25,
      force = 10, show.legend = FALSE,
      # family = "Times",
      fontface = "italic",
      aes(label = ifelse(labels == TRUE, # label top t's.
        as.character(.data[[idColumn]]), ""
      ))
    )
  }
  # Save
  ggsave(
    plot = volcano, path = newPath,
    filename = paste(prefix, "showname", showName, "volcano", fileEnd, sep = "."),
    width = 185, height = newHeight, units = "mm"
  )
}
