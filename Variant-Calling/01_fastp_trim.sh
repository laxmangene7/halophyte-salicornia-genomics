#!/bin/bash -l
#SBATCH --job-name=fastp_trim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-317
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

module load fastp

############################################################
# Directory containing paired-end FASTQ files
############################################################

DATA_DIR="/path/to/FASTQ"

############################################################
# Select sample
############################################################

FILES=(${DATA_DIR}/*_1.fq.gz)
INDEX=$((SLURM_ARRAY_TASK_ID-1))
SAMPLE=${FILES[$INDEX]}
PREFIX=$(basename "$SAMPLE" _1.fq.gz)

############################################################
# Run fastp
############################################################

fastp \
    -i ${DATA_DIR}/${PREFIX}_1.fq.gz \
    -I ${DATA_DIR}/${PREFIX}_2.fq.gz \
    -o ${PREFIX}_1.trimmed.fq.gz \
    -O ${PREFIX}_2.trimmed.fq.gz \
    --thread ${SLURM_CPUS_PER_TASK} \
    --detect_adapter_for_pe \
    --qualified_quality_phred 15 \
    --length_required 150 \
    --dedup \
    --html ${PREFIX}.fastp.html \
    --json ${PREFIX}.fastp.json
