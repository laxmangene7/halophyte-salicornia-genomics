#!/bin/bash -l
#SBATCH --job-name=hisat2_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=06-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-23              # Adjust according to the number of samples
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

############################################################
# Align paired-end Illumina reads using HISAT2
#
# Input:
#   file_fp.<sample>_1.fq.gz
#   file_fp.<sample>_2.fq.gz
#
# Reference:
#   Salicornia_France.S.bigelovii.combined.ref.noChrun.ref
#
# Output:
#   <sample>.sam
#   <sample>.log
#
# Parameters:
#   --no-spliced-alignment   Disable splice-aware alignment
#   --no-unal                Do not report unaligned reads
#   -X 1000                  Maximum paired-end fragment length
############################################################

module load hisat2

# Directory containing trimmed FASTQ files
FASTQ_DIR="/path/to/trimmed_fastq"

# HISAT2 reference index prefix
REFERENCE="/path/to/reference/Salicornia_France.S.bigelovii.combined.ref.noChrun.ref"

# Get sample for this SLURM array task
sample=$(ls ${FASTQ_DIR}/*_1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
base=$(basename "${sample}" "_1.fq.gz")

echo "Processing sample: ${base}"

hisat2 \
    -p ${SLURM_CPUS_PER_TASK} \
    -x ${REFERENCE} \
    -1 "${FASTQ_DIR}/${base}_1.fq.gz" \
    -2 "${FASTQ_DIR}/${base}_2.fq.gz" \
    -X 1000 \
    --no-spliced-alignment \
    --no-unal \
    -S "${base}.sam" \
    > "${base}.log" 2>&1

echo "Alignment completed for ${base}"
