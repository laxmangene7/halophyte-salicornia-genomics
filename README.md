# halophyte-salicornia-genomics
We assembled six chromosome-scale Salicornia genomes and analyzed 318 accessions from a resequencing panel. Our results reveal taxonomic inconsistencies and substantial genetic diversity, providing a genomic foundation for breeding and neo-domestication of this halophyte for saline agriculture.

## Phylogenetic tree
Run Salicornia.phylogenetic.tree.R 

Input files: 
- line.info.salicornia2025.txt
- kmer_matrix.319.salicornia_accessions.filtered.5.percent.maf.txt

## Genome Coverage

Run genome.coverage.sh 

Input files: 
- fastp trimmed paired-end reads (the raw reads are available at NCBI)
- bed files provided in Dryad
```
