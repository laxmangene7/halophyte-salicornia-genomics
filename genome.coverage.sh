## Genome Coverage
HISAT2 Trimmed Reads Alignment

```
#!/bin/bash -l
#SBATCH --job-name=run-hisat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

module load hisat2

dir="/ibex/scratch/projects/c2159/Genotyping_Salicornia/genome/updated_frozen.genome/all_trimmmed_files"

for sample in ${dir}/*_1.fq.gz
do
  base=$(basename "$sample" "_1.fq.gz")
  hisat2 -p 12 -x /ibex/scratch/projects/c2159/Genotyping_Salicornia/genome/updated_frozen.genome/new/combined_europaea.fr.S.bigelovii.ref \
         -1 ${dir}/${base}_1.fq.gz -2 ${dir}/${base}_2.fq.gz \
         -X 1000 \
         -S ${base}.octa.ref.sam \
         --no-spliced-alignment --no-unal &> ${base}.log
done
```

SAM to 100kb Bin Counts (sam_to_100kb_bins.sh) for bed/text files

#
