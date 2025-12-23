## computing polymorphism (Ref, Alt, Het) sites per individuals ##
#!/bin/bash -l
#SBATCH --job-name=ref.alt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

set -euo pipefail

module purge
module load bcftools/1.16

VCF="all.diploides.sinus.persica.ref.filt.vcf.gz"
PREFIX="per_sample_stats_$(date +%Y%m%d_%H%M%S)"

echo "[$(date)] Running bcftools stats"
bcftools stats \
    --threads 6 \
    -s - \
    "$VCF" > "${PREFIX}.stats"

echo "[$(date)] Extracting PSC lines"
grep "^PSC" "${PREFIX}.stats" > "${PREFIX}.psc"

echo "[$(date)] Creating TSV summary"
echo -e "Sample\tRef\tHet\tAlt\tMissing\tPercent_Poly" > "${PREFIX}.tsv"

awk '{
    sample = $3;   # sample name
    ref    = $4;   # homozygous reference
    het    = $5;   # heterozygous
    alt    = $6;   # homozygous alternate
    miss   = $14;  # MISSING genotypes (correct column!)

    nonmiss = ref + het + alt;
    percent = (nonmiss > 0 ? (het / nonmiss) * 100 : 0);

    printf "%s\t%d\t%d\t%d\t%d\t%.5f\n",
           sample, ref, het, alt, miss, percent
}' "${PREFIX}.psc" >> "${PREFIX}.tsv"

echo "[$(date)] DONE"
echo "Output files:"
echo "  ${PREFIX}.stats"
echo "  ${PREFIX}.psc"
echo "  ${PREFIX}.tsv"

echo "Preview:"
head -n 10 "${PREFIX}.tsv"
