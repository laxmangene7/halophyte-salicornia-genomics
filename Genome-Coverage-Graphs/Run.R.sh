#!/bin/bash -l
#SBATCH --job-name=plot_genome_coverage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

############################################################
# Run genome coverage plotting script
#
# Input:
#   *.norm.txt
#
# Output:
#   *_Normalized.reads.pdf
############################################################

module load R

Rscript R_script.plotting.R
