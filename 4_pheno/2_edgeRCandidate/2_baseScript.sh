#!/bin/bash
#SBATCH --job-name=
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=3GB  # Job memory request
#SBATCH --time=0-00:30:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=
#SBATCH --error=

# Run edgeR 
# Author: Amanda Zacharias
# Date: 2023-07-18
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Load dependencies
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

RSCRIPTNAME=
CANDPATH=
BASENAME=

# Execute R script
Rscript ../${RSCRIPTNAME} \
  --candPath $CANDPATH \
  --basename $BASENAME \
  --countPath 4_pheno/1_dataPrep/${DATATYPE}/filtData/norm.counts.csv \
  --vstCountPath 4_pheno/1_dataPrep/${DATATYPE}/filtData/vst.counts.csv \
  --workingDir 

echo Job ended at $(date +'%T')
