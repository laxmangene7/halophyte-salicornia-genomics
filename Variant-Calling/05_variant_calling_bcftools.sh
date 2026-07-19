#!/bin/bash -l
#SBATCH --job-name=bcftools_call
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6-00:00:00
#SBATCH --mem=200G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

module load bcftools

############################################################
# Input files
############################################################

REFERENCE="/path/to/Salicornia_ramosissima.ref.fa"

BAM_LIST="bamFile.list.txt"

############################################################
# Variant calling
############################################################

bcftools mpileup \
    -f ${REFERENCE} \
    -b ${BAM_LIST} \
    --annotate AD,DP,INFO/AD \
    --skip-indels \
    -B \
    -r chr8,chr9 \
| bcftools call \
    -m \
    --variants-only \
    --skip-variants indels \
    --group-samples \
    --output-type v \
    -o chr8_9.vcf
