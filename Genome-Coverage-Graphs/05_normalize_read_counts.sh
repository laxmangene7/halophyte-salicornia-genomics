#!/bin/bash -l
#SBATCH --job-name=normalize_reads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01-04:00:00
#SBATCH --mem=2G
#SBATCH --array=0-317              # Adjust according to the number of samples
#SBATCH --output="%x_%A_%a.out"
#SBATCH --error="%x_%A_%a.err"

############################################################
# Normalize 100-kb Genome Coverage
#
# This script normalizes read counts to account for
# differences in sequencing depth among samples.
#
# Normalization procedure:
#
# 1. Calculate the average read count of the top 10%
#    highest-coverage genomic bins.
#
# 2. Define a threshold equal to 10% of this average.
#
# 3. Calculate the average read count using only bins
#    with coverage above the threshold.
#
# 4. Normalize every genomic bin by this average.
#
# Input:
#   Sample list:
#       all.sample.txt
#
#   Coverage files:
#       *_100kb.txt
#
# Output:
#   <sample>.norm.txt
#
# Output columns:
#   1. Raw read count
#   2. Chromosome
#   3. 100-kb bin start position
#   4. Normalized read count
############################################################

# Read list of input files
readarray -t SAMPLE_LIST < all.sample.txt

# Select file corresponding to this SLURM array task
FILE="${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}"

# Output file prefix
PREFIX=$(basename "$FILE" "_100kb.txt")

echo "=========================================="
echo "Processing: ${FILE}"
echo "=========================================="

############################################################
# Step 1. Calculate average of top 10% highest-coverage bins
############################################################

TOTAL_LINES=$(wc -l < "$FILE")
TOP_LINES=$((TOTAL_LINES / 10))

TOP10_AVG=$(awk '{print $1}' "$FILE" \
    | sort -nr \
    | head -n "$TOP_LINES" \
    | awk '{sum+=$1} END {if(NR>0) print sum/NR}')

############################################################
# Step 2. Calculate threshold (10% of top10 average)
############################################################

THRESHOLD=$(echo "$TOP10_AVG * 0.1" | bc -l)

############################################################
# Step 3. Calculate usable average
############################################################

USABLE_AVG=$(awk -v T="$THRESHOLD" '
$1>=T{
    sum+=$1
    n++
}
END{
    if(n>0)
        print sum/n
    else
        print 1
}' "$FILE")

############################################################
# Step 4. Normalize read counts
############################################################

awk -v AVG="$USABLE_AVG" '
BEGIN{
    OFS="\t"
}
{
    print $1,$2,$3,$1/AVG
}' "$FILE" > "${PREFIX}.norm.txt"

echo "Finished: ${PREFIX}.norm.txt"
