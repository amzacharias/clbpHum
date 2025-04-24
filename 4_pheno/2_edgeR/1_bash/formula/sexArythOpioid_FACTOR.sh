#!/bin/bash
#SBATCH --job-name=sexArythOpioid_FACTOR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=3GB  # Job memory request
#SBATCH --time=0-00:30:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=sexArythOpioid_FACTOR.out
#SBTACH --error=sexArythOpioid_FACTOR.err

# What transcripts are different between opioid users and non-users?
# Amanda Zacharias
# March 13th, 2023

echo Job started at $(date +'%T')

# Load dependencies
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

# Execute R script
Rscript ../../0_sexArythOpioid_FACTOR.R \
  --basename "sexArythOpioid_FACTOR" \
  --countPath "4_pheno/1_dataPrep/filtData/norm.counts.csv" \
  --vstCountPath "4_pheno/1_dataPrep/filtData/vst.counts.csv" \
  --workingDir "absolutePath/clbpHum"

echo Job ended at $(date +'%T')
