#!/bin/bash
#SBATCH --job-name=hi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=10
#SBATCH --mem=55GB  # Job memory request
#SBATCH --time=1-8:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=hi
#SBTACH --error=hi

# Title: Aligning reads with Hisat2
# Author: Amanda Zacharias
# Date: 2023-07-11
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#assembly

echo Job started at $(date +'%T')

# Load dependencies
module load StdEnv/2020 samtools/1.10 hisat2/2.2.1 

# Set variables
PAIRED_END_1=
PAIRED_END_2=
INDEX=
ALIGN_PATH=
SUMMARY_PATH=
BATCH=
LANE=
SAMPLEID=

echo ALignment Started at $(date +'%T')
hisat2 -p 10 -x $INDEX -1 $PAIRED_END_1 -2 $PAIRED_END_2 \
  --dta --sensitive --no-discordant --no-mixed \
  --summary-file $SUMMARY_PATH \
  --time --verbose \
  --rg-id=${BATCH}.${LANE} \
  --rg PU:${BATCH}.${LANE} \
  --rg SM:${SAMPLEID} \
  --rg LB:LIB_${SAMPLEID} \
  --rg PL:DNBSEQ \
  -S ${OUT_PATH}.sam

echo Samtools processing started at $(date +'%T')
samtools view -b -@ 10 ${OUT_PATH}.sam > ${OUT_PATH}.bam
rm ${OUT_PATH}.sam
echo collate started at $(date +'%T')
samtools collate -@ 10 -o ${OUT_PATH}.col.bam  ${OUT_PATH}.bam 
rm ${OUT_PATH}.bam
echo fixmate started at $(date +'%T')
samtools fixmate -m -@ 10 ${OUT_PATH}.col.bam  ${OUT_PATH}.fix.bam 
rm ${OUT_PATH}.col.bam
echo sort started at $(date +'%T')
samtools sort -m 2G -@ 10 -o ${OUT_PATH}.sort.bam ${OUT_PATH}.fix.bam 
rm ${OUT_PATH}.fix.bam
echo markdup started at $(date +'%T')
samtools markdup -@ 10 -s ${OUT_PATH}.sort.bam ${OUT_PATH}.sort.mrkdup.bam 
rm ${OUT_PATH}.sort.bam
echo index started at $(date +'%T')
samtools index -b -@ 10 ${OUT_PATH}.sort.mrkdup.bam  ${OUT_PATH}.sort.mrkdup.bam.bai

echo Job ended at $(date +'%T')
