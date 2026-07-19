#!/bin/bash -l
#SBATCH --job-name=hisat2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=1-00:00:00
#SBATCH --mem=60G
#SBATCH --array=1-317
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

module load hisat2

############################################################
# Input files
############################################################

FASTQ_DIR="/path/to/trimmed_fastq"

REFERENCE="/path/to/Salicornia_ramosissima.ref"

############################################################
# Select sample
############################################################

FILES=(${FASTQ_DIR}/*_1.trimmed.fq.gz)
INDEX=$((SLURM_ARRAY_TASK_ID-1))
SAMPLE=${FILES[$INDEX]}
PREFIX=$(basename "$SAMPLE" _1.trimmed.fq.gz)

############################################################
# Align reads
############################################################

hisat2 \
    -p ${SLURM_CPUS_PER_TASK} \
    -x ${REFERENCE} \
    -1 ${FASTQ_DIR}/${PREFIX}_1.trimmed.fq.gz \
    -2 ${FASTQ_DIR}/${PREFIX}_2.trimmed.fq.gz \
    -X 1000 \
    --no-spliced-alignment \
    --no-unal \
    -S ${PREFIX}.sam \
    > ${PREFIX}.log 2>&1
