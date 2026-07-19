#!/bin/bash -l
#SBATCH --job-name=fastp_trim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-23              # Adjust according to the number of samples
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

############################################################
# Trim paired-end Illumina reads using fastp
#
# Input:
#   *_1.fq.gz
#   *_2.fq.gz
#
# Output:
#   file_fp.<sample>_1.fq.gz
#   file_fp.<sample>_2.fq.gz
#   <sample>.html
#   <sample>.json
#
# Performs:
#   - Adapter trimming
#   - Quality filtering
#   - Length filtering
#   - Duplicate read removal
############################################################

module load fastp

# Directory containing raw FASTQ files
DATA_DIR="/ibex/scratch/projects/c2159/Genotyping_Salicornia/all_fastqs"

# Get sample corresponding to this SLURM array task
FILES=(${DATA_DIR}/*_1.fq.gz)
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${FILES[$INDEX]}

# Extract sample name
PREFIX=$(basename "${SAMPLE}" "_1.fq.gz")

echo "Processing sample: ${PREFIX}"

fastp \
    -i "${DATA_DIR}/${PREFIX}_1.fq.gz" \
    -I "${DATA_DIR}/${PREFIX}_2.fq.gz" \
    -o "file_fp.${PREFIX}_1.fq.gz" \
    -O "file_fp.${PREFIX}_2.fq.gz" \
    --thread ${SLURM_CPUS_PER_TASK} \
    --detect_adapter_for_pe \
    --qualified_quality_phred 15 \
    --length_required 150 \
    --dedup \
    --html "${PREFIX}.html" \
    --json "${PREFIX}.json"

echo "Finished processing ${PREFIX}"
