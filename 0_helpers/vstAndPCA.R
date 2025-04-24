#!/usr/bin/env Rscript
#-------------------------------------------------
# Title:
# Author: Amanda Zacharias
# Date: 2023-07-17
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
#
#
#
# Options -----------------------------------------


# Packages -----------------------------------------


# Functions -----------------------------------------
GetVst <- function(df, coldat, varOfInterest) {
  dds <- DESeqDataSetFromMatrix(
    round(data.matrix(df)),
    coldat,
    as.formula(paste("~", varOfInterest))
  )
  dds <- estimateSizeFactors(dds)
  vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  return(vst)
}

MakePCA <- function(vstData, intGroups, labels, newFilename, newPath) {
  pcaData <- plotPCA(vstData,
    intgroup = intGroups,
    returnData = TRUE
  ) %>%
    mutate("labels" = labels) %>%
    mutate_at(intGroups, as.factor)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  for (group in c(intGroups, "labels")) {
    if (sum(labels) >= 1 & group != "labels") {
      pca_gplot <- data.frame(pcaData) %>%
        ggplot(aes(PC1, PC2, color = labels, fill = .data[[group]])) +
        geom_point(shape = 21, size = 4, alpha = 0.6, stroke = 1.5) +
        scale_color_manual(
          name = "Is outlier", labels = c("no", "yes"),
          values = c("grey95", "black")
        ) +
        guides(color = guide_legend(order = 1, byrow = TRUE), fill = guide_legend(order = 2))
    } else {
      pca_gplot <- data.frame(pcaData) %>%
        ggplot(aes(PC1, PC2, color = .data[[group]])) +
        geom_point(size = 4, alpha = 0.6)
    }
    pca_gplot <- pca_gplot +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle(str_to_sentence(group)) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(-1, "mm"),
        text = element_text(size = 15)
      )
    if (length(groups) > 6) {
      pca_gplot <- theme(legend.position = "none")
    }
    if (group == "phenotype"){
      lbls <- c("constant (h)", "constant (l)", "mixed", "rhythmic")
      colors <- c("#F09F39", "#731B0D", "#934DC9", "#D26339")
      names(colors) <- c("ConstantHigh", "ConstantLow", "Mixed", "Rhythmic")
      pca_gplot <- pca_gplot + 
        scale_colour_manual(name = "Phenotype", labels = lbls,
                            values = colors, aesthetics = "colour")
    }
    ggsave(
      plot = pca_gplot, path = newPath,
      filename = paste(group, newFilename, sep = "."),
      width = 185, height = 185, units = "mm", dpi = 300
    )
  }
}
