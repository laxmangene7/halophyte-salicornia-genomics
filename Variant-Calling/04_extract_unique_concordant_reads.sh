#!/bin/bash -l
#SBATCH --job-name=unique_reads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3:00:00
#SBATCH --mem=35G
#SBATCH --array=1-317
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

module load samtools

############################################################
# Select SAM file
############################################################

SAMPLE=$(ls *.sam | sed -n ${SLURM_ARRAY_TASK_ID}p)
PREFIX=$(basename ${SAMPLE} .sam)

############################################################
# Retain uniquely mapped concordant reads
############################################################

cat \
<(samtools view -H ${PREFIX}.sam) \
<(awk '/YT:Z:CP/ && /NH:i:1/' ${PREFIX}.sam) \
> ${PREFIX}.bam
