#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Use IsoformAnalyzeR to convert transcript abundance --> gene counts
# Author: Amanda Zacharias
# Date: 2023-07-06
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load StdEnv/2020 r/4.2.1
# View section "Rescue StringTie Annotation and Extract Gene Count Matrix" in vignette
# https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#importing-data-into-r
# The newer the version of IsoformSwitchAnalyzeR, the better !
# To run in the terminal: 2 options
# R_LIBS_USER=absolutePath/R/x86_64-redhat-linux-gnu-library/4.2.1 R
# OR
# R_LIBS_USER=absolutePath/R/x86_64-redhat-linux-gnu-library/4.2.1 Rscript 7_isoformAnalyzeR.R
#
# Options -----------------------------------------

# Packages -----------------------------------------
library(dplyr) # 1.1.0x
# BiocManager::install("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR) # 1.17.4

# Pathways -----------------------------------------
# Warning: using an absolute path!
setwd("absolutePath/clbpHum")
stringtieDir <- file.path(getwd(), "3_stringtie")
# Input ===========
coldataPath <- file.path("coldata.csv")
mergeGtfPath <- file.path(stringtieDir, "4_merge", "merged.gtf")
inputDir <- file.path(stringtieDir, "ballgown")

# Output ===========
isoDir <- file.path(stringtieDir, "isoformAnalyzeR")
system(paste("mkdir", isoDir))

# Load coldata --------------------------------------------------------------
# Coldata
coldata <- read.csv(coldataPath,
  row.names = 1, stringsAsFactors = F
)

# Format for isoformSwitchAnalyzeR
isoColdata <- data.frame(
  sampleID = coldata$filename,
  condition = coldata$phenotype
)

# Get counts function -----------------------------------------
#' Use IsoformSwitchAnalyzeR functions to get gene counts from transcript abundances
#'
#' @param inputDir A string path to the folder with ballgown files from StringTie.
#' @param readlength The numeric length of sequencing reads.
#' @param coldat Metadata dataframe explaining what each sample is.
#' @param gtfPath A string path to the reference gtf used for assembly.
#' @returns A named list.
#' @examples
#' countObj <- getCounts("./ballgown", 51, df, "./merge/merged.gtf")
getCounts <- function(inputDir, readlength, coldat, gtfPath) {
  quant <- importIsoformExpression(
    parentDir = inputDir,
    readLength = readlength,
    addIsofomIdAsColumn = TRUE
  )
  switchQuant <- importRdata(
    isoformCountMatrix = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix = coldat,
    isoformExonAnnoation = gtfPath,
    showProgress = TRUE
  )
  geneCounts <- extractGeneExpression(
    switchQuant,
    extractCounts = TRUE # set to FALSE for abundances
  )
  return(list("switchQuant" = switchQuant, "geneCounts" = geneCounts))
}

# Execute function -----------------------------------------
# Get counts
countObj <- getCounts(
  inputDir = inputDir,
  readlength = 150,
  coldat = isoColdata,
  gtfPath = mergeGtfPath
)
# Save files
write.csv(
  countObj$geneCounts,
  file.path(isoDir, "geneCounts.csv")
)
saveRDS(countObj$switchQuant,
  file = file.path(isoDir, "switchQuant.rds")
)

# Transcript id to gene id -----------------------------------------
GetT2G <- function(switchQuant) {
  #' Extract isoform_id to gene_id/gene_name mapping from switchQuant object
  #'
  #' @param switchQuant A switchQuant object from isoformSwitchAnalyzeR
  #' @returns A dataframe.
  #' @note Need a separate t2g for each dataframe because, for example,
  #' MSTRG.1 != MSTRG.1 in a different dataset.
  #' @examples
  #' t2g <- GetT2G(countsObj$switchQuant)
  t2g <- switchQuant$isoformFeatures %>%
    dplyr::select(c("isoform_id", "gene_id", "gene_name")) %>%
    distinct(isoform_id, gene_id, gene_name)
  return(t2g)
}
t2g <- GetT2G(countObj$switchQuant)
write.csv(t2g, file.path(isoDir, "t2g.csv"))

# Output messages -----------------------------------------
# Step 1 of 3: Identifying which algorithm was used...
# The quantification algorithm used was: StringTie
# Found 120 quantification file(s) of interest
# Step 2 of 3: Reading data...
# reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120
# Step 3 of 3: Normalizing abundance values (not counts) via edgeR...
# Done
#
# Step 1 of 7: Checking data...
# Step 2 of 7: Obtaining annotation...
# importing GTF (this may take a while)...
# 39339 ( 11.95%) isoforms were removed since they were not expressed in any samples.
# Step 3 of 7: Fixing StringTie gene annoation problems...
# 38165 isoforms were assigned the ref_gene_id and gene_name of their associated gene_id.
# This was only done when the parent gene_id were associated with a single ref_gene_id/gene_name.
# 32281 isoforms were assigned the ref_gene_id and gene_name of the most similar
# annotated isoform (defined via overlap in genomic exon coordinates).
# This was only done if the overlap met the requriements
# indicated by the three fixStringTieViaOverlap* arguments.
# We were unable to assign 2196 isoforms (located within annotated genes) to a known ref_gene_id/gene_name.
# These were removed to enable analysis of the rest of the isoform from within the merged genes.
# 3511 gene_ids which were associated with multiple ref_gene_id/gene_names
# were split into mutliple genes via their ref_gene_id/gene_names.
# 52370 genes_id were assigned their original gene_id instead of the StringTie gene_id.
# This was only done when it could be done unambiguous.
# Step 4 of 7: Calculating gene expression and isoform fractions...
# Step 5 of 7: Merging gene and isoform expression...
# |======================================================================| 100%
# Step 6 of 7: Making comparisons...
# |======================================================================| 100%
# Step 7 of 7: Making switchAnalyzeRlist object...
# The GUESSTIMATED number of genes with differential isoform usage are:
#   comparison estimated_genes_with_dtu
# 1 ConstantHigh vs Mixed                  37 - 61
# 2  ConstantLow vs Mixed                   7 - 12
# 3     Mixed vs Rhythmic              1406 - 2343
# Done
#
# Warning messages:
#   1: In importRdata(isoformCountMatrix = quant$counts, isoformRepExpression = quant$abundance,  :
#                       We found 7649 (2.27%) unstranded transcripts.
#                     These were removed as unstranded transcripts cannot be analysed
#                     2: In importRdata(isoformCountMatrix = quant$counts, isoformRepExpression = quant$abundance,  :
#                                         No CDS annotation was found in the GTF files meaning ORFs could not be annotated.
#                                       (But ORFs can still be predicted with the analyzeORF() function)
