## 100 kb Genome Coverage Calculation from SAM Files

# Pipeline:
# 1. Remove header lines
# 2. Keep concordant paired reads (YT:Z:CP)
# 3. Keep uniquely mapped reads (NH:i:1)
# 4. Extract chromosome and position
# 5. Assign reads to 100 kb bins
# 6. Count reads per bin


#!/bin/bash -l
#SBATCH --job-name=samstotxt1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=02-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-1
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

module load samtools

# Directory containing SAM files
dir="/path/to/sam_files"
sample=$(ls ${dir}/*.sam | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$sample" ".sam")

echo "Processing sample: $base"

grep -v "^@" "$sample" \
    | grep YT:Z:CP \
    | grep NH:i:1 \
    | cut -f 3,4 \
    | sort -k1,1 -k2,2n \
    | awk '{print $1 "\t" int($2 / 100000) * 100000}' \
    | uniq -c \
    | awk '{print $1"\t"$2"\t"$3}' \
    > "${base}_100kb.txt"
