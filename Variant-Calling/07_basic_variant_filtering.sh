#!/bin/bash -l
#SBATCH --job-name=variant_filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=60G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

module load bcftools
module load tabix

############################################################
# Input VCF
############################################################

INPUT_VCF="chr3.chr4.vcf"

############################################################
# Output VCF
############################################################

OUTPUT_VCF="chr3_4.filtered.vcf.gz"

############################################################
# Basic variant filtering
#
# Retain variants satisfying:
#   QUAL >= 30
#   INFO/DP >= 20
#   0.01 < AF < 0.99
############################################################

bcftools filter \
    ${INPUT_VCF} \
    -i 'QUAL>=30 && INFO/DP>=20 && AF>0.01 && AF<0.99' \
    -Oz \
    -o ${OUTPUT_VCF}

############################################################
# Index filtered VCF
############################################################

tabix -p vcf ${OUTPUT_VCF}
