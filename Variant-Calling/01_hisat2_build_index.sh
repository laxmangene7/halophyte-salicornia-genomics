#!/bin/bash -l
#SBATCH --job-name=hisat2_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=4-00:00:00
#SBATCH --mem=200G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

############################################################
# Build HISAT2 reference genome index
#
# Input:
#   Salicornia_ramosissima.ref.fa
#
# Output:
#   Salicornia_ramosissima.ref.*.ht2
#
# HISAT2 index files are required before aligning
# Illumina paired-end reads to the reference genome.
############################################################

module load hisat2

REFERENCE="Salicornia_ramosissima.ref.fa"
PREFIX="Salicornia_ramosissima.ref"

hisat2-build \
    ${REFERENCE} \
    ${PREFIX}
