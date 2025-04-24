#!/bin/bash
#SBATCH --job-name=mq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=mq
#SBTACH --error=mq

# Title: MultiQC to summarize FastQC results
# Author: Amanda Zacharias
# Date: 2023-07-10
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Code
# Dependencies
module load StdEnv/2020 python/3.9.6
#pip install --user multiqc
#pip install --user --upgrade multiqc

# Variables
OUTDIR=
FILENAME=
FQPATHS=

# Begin MultiQC
multiqc \
  --outdir $OUTDIR \
  --filename $FILENAME \
  --force \
  --interactive \
  --cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
  $FQPATHS

echo Job ended at $(date +'%T')
