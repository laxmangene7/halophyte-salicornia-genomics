# halophyte-salicornia-genomics
We assembled six chromosome-scale Salicornia genomes and analyzed 318 accessions from a resequencing panel. Our results reveal taxonomic inconsistencies and substantial genetic diversity, providing a genomic foundation for breeding and neo-domestication of this halophyte for saline agriculture.

## Phylogenetic tree
Run Salicornia.phylogenetic.tree.R 

Input files: 
- line.info.salicornia2025.txt
- kmer_matrix.319.salicornia_accessions.filtered.5.percent.maf.txt

## Genome Coverage
HISAT2 Alignment Script 

```#!/bin/bash -l
#SBATCH --job-name=run-hisat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

module load hisat2

for sample in `ls /ibex/scratch/projects/c2159/Genotyping_Salicornia/genome/updated_frozen.genome/all_trimmmed_files/*_1.fq.gz`
do
  dir="/ibex/scratch/projects/c2159/Genotyping_Salicornia/genome/updated_frozen.genome/all_trimmmed_files"
  base=$(basename $sample "_1.fq.gz")
  hisat2 -p 12 -x /ibex/scratch/projects/c2159/Genotyping_Salicornia/genome/updated_frozen.genome/new/combined_europaea.fr.S.bigelovii.ref \
         -1 ${dir}/${base}_1.fq.gz -2 ${dir}/${base}_2.fq.gz \
         -S ${base}.octa.ref.sam \
         --no-spliced-alignment --no-unal &> ${base}.log
done```



## Population Structure

