#!/bin/bash -l
#SBATCH --job-name=sort_bam
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH --mem=35G
#SBATCH --array=1-317
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

module load samtools

############################################################
# BAM file list
############################################################

LIST=all.bamfile.list2.txt

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST})
PREFIX=$(basename ${SAMPLE} .bam)

############################################################
# Sort BAM
############################################################

samtools sort \
    ${SAMPLE} \
    -o ${PREFIX}.sorted.bam

############################################################
# Index BAM
############################################################

samtools index -c \
    ${PREFIX}.sorted.bam
