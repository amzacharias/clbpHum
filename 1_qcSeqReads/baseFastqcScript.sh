#!/bin/bash
#SBATCH --job-name=fq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB  # Job memory request
#SBATCH --time=0-12:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=fq
#SBTACH --error=fq

# Title: This script uses fastqc to check qualtiy of sequencing reads
# Author: Amanda Zacharias
# Date: 2023-07-10
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 

echo Job started at $(date +'%T')

module load StdEnv/2020
module load nixpkgs/16.09
module load fastqc/0.11.9
module list

INDATAPATH=
OUTDIR=
fastqc -f fastq -o $OUTDIR $INDATAPATH

echo Job ended at $(date +'%T')
