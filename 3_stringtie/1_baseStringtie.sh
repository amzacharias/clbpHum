#!/bin/bash
#SBATCH --job-name=st
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=st
#SBTACH --error=st

# Title: StringTie to assemble transcripts, pass 1
# Author: Amanda Zacharias
# Date: 2023-07-12
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Load dependencies
module load StdEnv/2020 stringtie/2.1.5

# Set variables
INPUT=
OUT_GTF=
REF_GTF=

echo Stringtie Started at $(date +'%T')
stringtie $INPUT -p 5 -G $REF_GTF -o $OUT_GTF

echo Job ended at $(date +'%T')
