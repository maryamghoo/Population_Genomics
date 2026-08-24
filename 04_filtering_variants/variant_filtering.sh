#!/bin/bash -l
#SBATCH -J variantFiltering
#SBATCH -o /variant_filtering_results/logs/filtered_variants.out
#SBATCH -e /variant_filtering_results/logs/filtered_variants.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 48:00:00

cd /path/to/output_directory

# Filtering Variants
/path/to/bcftools \
filter \
--threads 16 \
-sFAIL -e'QUAL < 30 || INFO/MQ <= 30 || MEAN(FMT/DP) > 21 || MEAN(FMT/DP) < 2 || INFO/RPBZ < -5 || INFO/RPBZ > 5 || INFO/BQBZ < -5 || INFO/DP4[2]+INFO/DP4[3] <= 2' \
-g10 \
-G10 \
-o filteredvariants.vcf \
/path/to/variants.vcf

# Keep only PASS variants
/path/to/bcftools \
view -f 'PASS' \
-O v filteredvariants.vcf > passedvariants.tmp.vcf

# Rename SNP IDs
awk 'BEGIN {x=1} {OFS="\t"} !/^#/ {$3="SNP_"x++} {print}' \
passedvariants.tmp.vcf > passedvariants.vcf

# Remove temporary file
rm passedvariants.tmp.vcf

# Generate stats for passed variants
/path/to/bcftools \
stats passedvariants.vcf > passedvariants.stats
