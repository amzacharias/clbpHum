#!/bin/bash
#SBATCH --job-name=runPrepDE
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-00:30:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=runPrepDE.out
#SBTACH --error=runPrepDE.err

# Title: Use `prepDE.py` to get transcript/gene counts
# Author: Amanda Zacharias
# Date: 2023-07-12
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

module load StdEnv/2020 python/3.9.6

STRING_PATH=absolutePath/clbpHum/3_stringtie
PREPDE_PATH=${STRING_PATH}/prepDE.py

python $PREPDE_PATH -v \
  -i ${STRING_PATH}/gtfLists/pass2List.txt \
  -g ${STRING_PATH}/prepDEcounts/genes.csv \
  -t ${STRING_PATH}/prepDEcounts/transcripts.csv \
  -l 150

echo Job ended at $(date +'%T')
