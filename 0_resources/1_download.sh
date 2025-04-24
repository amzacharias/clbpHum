#!/bin/bash
#SBATCH --job-name=download
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=150GB  # Job memory request
#SBATCH --time=00-6:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=download.out
#SBTACH --error=download.err

# Download reference from gencode
cd absolutePath/clbpHum/0_resources

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz

gunzip -v *.gtf.gz
gunzip -v *.fa.gz
