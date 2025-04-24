#!/bin/bash
#SBATCH --job-name=wgcnaSTT
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB  # Job memory request
#SBATCH --time=1-04:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=wgcnaSTT.out
#SBTACH --error=wgcnaSTT.err

# Run WGCNA script
# Author: Amanda Zacharias
# Date: 2023-07-18
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Load dependencies
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

# Run R script
Rscript ./1_makeNetworkStepwise.R # 1 cpu, 800GB memory, 1-8:00:00
Rscript ./2_modTraitRelationGLM.R
Rscript ./2_modTraitRelationMLR.R
Rscript ./3_gprofilerWGCNA.R
Rscript ./3_export2cytoscape.R # 500GB, 0-4:00:00

echo Job ended at $(date +'%T')
