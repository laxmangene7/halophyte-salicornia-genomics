#!/bin/bash -l
#SBATCH --job-name=hisat2_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=04-00:00:00
#SBATCH --mem=200G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

############################################################
# Build HISAT2 index
#
# Input:
#   Salicornia_France.S.bigelovii.combined.ref.noChrun.ref.fa
#
# Output:
#   Salicornia_France.S.bigelovii.combined.ref.noChrun.ref.*.ht2
#
# This index is used for aligning Illumina whole-genome
# resequencing reads prior to genome coverage analysis.
############################################################

module load hisat2

# Reference genome FASTA
REFERENCE="Salicornia_France.S.bigelovii.combined.ref.noChrun.ref.fa"

# Prefix for HISAT2 index files
PREFIX="Salicornia_France.S.bigelovii.combined.ref.noChrun.ref"

echo "Building HISAT2 index..."
echo "Reference: ${REFERENCE}"

hisat2-build \
    ${REFERENCE} \
    ${PREFIX}

echo "Index construction completed."
