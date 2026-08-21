#!/bin/bash
#SBATCH -J callable_vcf
#SBATCH -o /callable_vcf/logs/callable_vcf.out
#SBATCH -e /callable_vcf/logs/callable_vcf.err
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00

# Stop on errors, undefined variables, and failed pipes
set -euo pipefail

##### Settings #####

VCF="/path/to/input_snps.vcf.gz"
BED="/path/to/callable_common.merged.bed"

OUT="/path/to/snps_callable.vcf.gz"

##########

echo "Filtering SNPs to callable genome..."

bcftools view -R $BED $VCF -Oz -o $OUT

echo "Indexing VCF..."

bcftools index $OUT

echo "Finished."
