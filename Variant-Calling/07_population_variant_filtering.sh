#!/bin/bash -l
#SBATCH --job-name=population_filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=2-10:00:00
#SBATCH --mem=300G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

############################################################
# Population-level SNP filtering
#
# This script filters variants based on:
#   - Minor allele frequency (MAF)
#   - Missing genotype rate
#   - Heterozygosity
#
# Filtering thresholds:
#   MAF            >= 0.02
#   Missing data   <= 20%
#   Heterozygosity <= 20%
############################################################

INPUT_VCF="chr3_4.filtered.vcf"

MAF=0.02
MISSING=20
HET=20

perl vcf_MAF_Missing_Het_filter.pl \
    ${INPUT_VCF} \
    ${MAF} \
    ${MISSING} \
    ${HET}
