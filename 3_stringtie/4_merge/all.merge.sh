#!/bin/bash
#SBATCH --job-name=mergeall
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=5
#SBATCH --mem=3GB  # Job memory request
#SBATCH --time=0-2:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=mergeall.out
#SBTACH --error=mergeall.err

# Title: This script uses StringTie to merge assembled transcripts.
# Author: Amanda Zacharias
# Date: 2023-07-12
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +%T)

# Load dependencies
module load StdEnv/2020 stringtie/2.1.5
   
# Set variables
GTFS_LIST=absolutePath/clbpHum/3_stringtie/gtfLists/pass1List.txt
OUTPUT=absolutePath/clbpHum/3_stringtie/4_merge/merged.gtf
REF_GTF=absolutePath/clbpHum/0_resources/gencode.v41.annotation.gtf

echo Stringtie Started at $(date +%T)
stringtie --merge -p 20 -o $OUTPUT -G $REF_GTF $GTFS_LIST

echo Job ended at $(date +%T)
