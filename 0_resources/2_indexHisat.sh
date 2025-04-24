#!/bin/bash
#SBATCH --job-name=indexHisat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=5
#SBATCH --mem=200GB  # Job memory request
#SBATCH --time=10-00:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=indexHisat.out
#SBTACH --error=indexHisat.err

# Building HGFM index with transcript
# Amanda Zacharias, April 25th, 2022
# Reference link used to make this script: 
# http://daehwankimlab.github.io/hisat2/howto/ 

echo "job started"

# Load HISAT2 and dependencies
module load StdEnv/2020
module load hisat2/2.2.1

# Assign variables
BASEPATH=absolutePath/clbpHum/0_resources
FNA_HUMAN_REF=${BASEPATH}/GRCh38.primary_assembly.genome.fa
GTF_HUMAN_REF=${BASEPATH}/gencode.v41.annotation.gtf
SPLICESITE_PATH=${BASEPATH}/gencode.v41.annotation.ss
EXON_PATH=${BASEPATH}/gencode.v41.annotation.exon
INDEX_NAME=${BASEPATH}/index/idx

# Make folder for index files
mkdir ${BASEPATH}/index

# Extract splice sites and exons
echo "Started extracting splice sites"
hisat2_extract_splice_sites.py -v $GTF_HUMAN_REF > $SPLICESITE_PATH
echo "Started extracting exons"
hisat2_extract_exons.py -v $GTF_HUMAN_REF > $EXON_PATH

# Build index
echo "Started building the index"
hisat2-build -p 20 -f $FNA_HUMAN_REF --ss $SPLICESITE_PATH --exon $EXON_PATH $INDEX_NAME

echo "job ended"
