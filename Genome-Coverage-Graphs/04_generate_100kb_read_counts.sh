#!/bin/bash -l
#SBATCH --job-name=coverage_100kb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-23              # Adjust according to the number of samples
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

############################################################
# Generate 100-kb Genome Coverage
#
# Input:
#   <sample>.sam
#
# Output:
#   <sample>_100kb.txt
#
# Pipeline:
#   1. Remove SAM header lines
#   2. Retain concordant paired-end reads (YT:Z:CP)
#   3. Retain uniquely mapped reads (NH:i:1)
#   4. Extract chromosome and alignment position
#   5. Assign reads to non-overlapping 100-kb bins
#   6. Count reads per genomic bin
#
# Output format:
#   Column 1: Number of reads
#   Column 2: Chromosome
#   Column 3: 100-kb bin start position
############################################################

module load samtools

# Directory containing SAM files
SAM_DIR="/path/to/sam_files"

# Get sample corresponding to this SLURM array task
sample=$(ls ${SAM_DIR}/*.sam | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "${sample}" ".sam")

echo "========================================="
echo "Processing sample: ${base}"
echo "========================================="

grep -v "^@" "${sample}" \
    | grep "YT:Z:CP" \
    | grep "NH:i:1" \
    | cut -f3,4 \
    | sort -k1,1 -k2,2n \
    | awk '{print $1 "\t" int($2/100000)*100000}' \
    | uniq -c \
    | awk '{print $1"\t"$2"\t"$3}' \
    > "${base}_100kb.txt"

echo "Finished processing ${base}"
echo "Output: ${base}_100kb.txt"
