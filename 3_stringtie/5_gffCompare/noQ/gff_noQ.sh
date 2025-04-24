#!/bin/bash
#SBATCH --job-name=gff_noQ
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB  # Job memory request
#SBATCH --time=0-02:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=gff_noQ.out
#SBTACH --error=gff_noQ.err

# This script is used to evaluate the assembly of transcripts
# by StringTie. 

# Amanda Zacharias
# July 14th, 2022

echo "job started"

# Load dependencies
module load nixpkgs/16.09
module load gffcompare/0.11.6

# Set variables
BASEPATH=absolutePath/clbpHum

ORIGINAL_GTF=${BASEPATH}/0_resources/gencode.v41.annotation.gtf
CUSTOM_GTF=${BASEPATH}/3_stringtie/4_merge/merged.gtf
OUT_PREFIX=/${BASEPATH}/3_stringtie/5_gffCompare/noQ/stats/strtcmp_noQ

# Actually do the comparison
gffcompare -T -R -r $ORIGINAL_GTF -o $OUT_PREFIX $CUSTOM_GTF

echo "job ended"
